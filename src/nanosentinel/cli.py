import argparse
import sys
import time
from pathlib import Path

from . import __version__
from .concat import concatenate
from .console import Console
from .errors import RunNotReady, SentinelError
from .references import prepare_references
from .report import hms
from .rundir import DEFAULT_MINKNOW_DATA, find_fastq_pass
from .samplesheet import read_samplesheet, select_run


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="nanosentinel",
        description="Monitor a live ONT run from a samplesheet: map reads as MinKNOW writes them, show coverage on a local "
        "dashboard, record how long each sample takes to reach a target depth, and (optionally) concatenate reads.",
    )
    p.add_argument("--version", action="version", version=f"nanosentinel {__version__}")
    g = p.add_argument_group("run")
    g.add_argument("--input", required=True, help="Samplesheet (CSV/TSV): needs 'sample' and 'barcode'; optional 'run_name', 'neg_control'")
    g.add_argument("--run_name", help="Run to monitor: a folder name under --minknow_data and the run_name column. Inferred if the sheet has one run")
    g.add_argument("--run_dir", help="Exact MinKNOW run folder (or its fastq_pass folder); overrides the --run_name search")
    g.add_argument("--minknow_data", default=DEFAULT_MINKNOW_DATA, help="MinKNOW data directory (default: %(default)s)")
    g.add_argument("--outdir", help="Output directory (default: nanosentinel_<run_name>)")
    g = p.add_argument_group("references (ordinary FASTA / BED / GFF files)")
    g.add_argument("--possible_references", help="Multi-FASTA reference panel reads are mapped to (several similar references are fine)")
    g.add_argument("--reference_fasta", help="Reference the primer scheme and genes are drawn against (default: first panel entry)")
    g.add_argument("--primer_bed", help="Primer scheme BED (ARTIC-style <amplicon>_LEFT/_RIGHT names) on the coordinate reference")
    g.add_argument("--gff", help="GFF3 with gene/CDS features on the coordinate reference (labels the coverage plot)")
    g = p.add_argument_group("monitoring")
    g.add_argument("--target_depth", type=int, default=20, help="Depth that counts as covered (default: %(default)s)")
    g.add_argument("--target_fraction", type=float, default=0.9, help="Fraction of the target region that must reach --target_depth (default: %(default)s)")
    g.add_argument("--denominator", choices=["genome", "amplicons"], default="genome", help="Region the fraction is measured over; 'amplicons' uses the --primer_bed span (default: %(default)s)")
    g.add_argument("--poll_interval", type=int, default=5, help="Seconds between updates (default: %(default)s)")
    g.add_argument("--port", type=int, default=6010, help="Dashboard port (default: %(default)s)")
    g.add_argument("--host", default="127.0.0.1", help="Address the dashboard listens on; use 0.0.0.0 to view it from another machine (default: %(default)s)")
    g.add_argument("--min_identity", type=float, default=0.8, help="Minimum alignment identity for a read to count as mapped (default: %(default)s)")
    g.add_argument("--threads", type=int, default=4, help="Threads for minimap2 (default: %(default)s)")
    g.add_argument("--webhook", help="URL to POST a JSON {'text': ...} message to when every sample has reached target")
    g.add_argument("--no_browser", action="store_true", help="Do not open the dashboard in a browser")
    g.add_argument("--clear_annotations", action="store_true", help="Discard earlier annotations and map every read again")
    g.add_argument("--skip_monitor", action="store_true", help="Skip monitoring (use with --concatenate_reads for a finished run)")
    g.add_argument("--dry_run", action="store_true", help="Check inputs and prepare the reference index, then stop without monitoring")
    g = p.add_argument_group("read concatenation")
    g.add_argument("--concatenate_reads", action="store_true", help="After monitoring, concatenate each sample's reads and write samplesheet_with_reads.csv")
    g.add_argument("--gzip", action="store_true", help="Gzip concatenated output when the input FASTQs are plain text")
    g.add_argument("--overwrite", action="store_true", help="Replace existing concatenated files")
    args = p.parse_args(argv)
    if args.skip_monitor and not args.concatenate_reads:
        p.error("--skip_monitor does nothing unless --concatenate_reads is also given")
    if not args.skip_monitor and not (args.possible_references or args.reference_fasta):
        p.error("monitoring needs --possible_references and/or --reference_fasta")
    if args.denominator == "amplicons" and not args.primer_bed:
        p.error("--denominator amplicons needs --primer_bed")
    if not 1 <= args.target_depth <= 255:
        p.error("--target_depth must be between 1 and 255")
    if not 0 < args.target_fraction <= 1:
        p.error("--target_fraction must be in (0, 1]")
    return args


def locate_reads(args, run_name, console, wait):
    """Find fastq_pass; when monitoring, wait for MinKNOW to create it if the run has not started yet."""
    started, announced = time.time(), False
    while True:
        try:
            return find_fastq_pass(args.run_dir, args.minknow_data, run_name)
        except RunNotReady as exc:
            if not wait:
                raise
            if not announced:
                console.say(f"Waiting for MinKNOW to create the run's fastq_pass folder ({exc}). Ctrl+C to cancel.")
                announced = True
            console.status(f"waiting for fastq_pass ... {hms(time.time() - started)}")
            time.sleep(5)


def run(args):
    header, rows = read_samplesheet(args.input)
    run_name, rows = select_run(rows, args.run_name)
    outdir = Path(args.outdir or f"nanosentinel_{run_name or 'run'}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    console = Console(outdir / "nanosentinel.log")

    n_ctrl = sum(r.neg_control for r in rows)
    console.say(f"nanosentinel {__version__} | run {run_name or '(unnamed)'} | {len(rows) - n_ctrl} sample(s), {n_ctrl} control(s)")
    fastq_pass = locate_reads(args, run_name, console, wait=not (args.skip_monitor or args.dry_run))
    console.end_status()
    console.say(f"Reads: {fastq_pass}")
    console.say(f"Output: {outdir}")

    if not args.skip_monitor:
        from .monitor import monitor  # imported late so --help and --skip_monitor stay fast

        refinfo = prepare_references(outdir, args.possible_references, args.reference_fasta, args.primer_bed, args.gff)
        console.say(f"Reference panel: {len(refinfo.lengths)} sequence(s); primer scheme and genes are drawn against {refinfo.coordinate_ref}")
        if args.dry_run:
            console.say("Dry run: inputs are valid and the reference index was built; nothing was monitored.")
            return 0
        monitor(
            console, rows, refinfo, fastq_pass, outdir, run_name, args.host, args.port, args.target_depth,
            args.target_fraction, args.denominator, args.poll_interval, args.webhook, not args.no_browser,
            args.clear_annotations, args.min_identity, args.threads,
        )  # fmt: skip
    elif args.dry_run:
        console.say("Dry run: inputs are valid; nothing was done.")
        return 0

    if args.concatenate_reads:
        console.say("Concatenating reads ...")
        sheet, n = concatenate(rows, header, fastq_pass, outdir, run_name, args.gzip, args.overwrite, log=console.warn)
        console.say(f"Wrote {n} FASTQ file(s) to {outdir / 'reads'}")
        console.say(f"Samplesheet with read paths: {sheet}")
    return 0


def main(argv=None):
    args = parse_args(argv)
    try:
        return run(args)
    except SentinelError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
