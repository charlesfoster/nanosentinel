import argparse
import sys
import time
from pathlib import Path

from . import __version__
from .concat import concatenate
from .console import Console, error, notice
from .errors import RunNotReady, SentinelError
from .protocol import build_protocol, load_protocol, parse_region
from .references import DEFAULT_GROUP_PATTERN, check_region, prepare_references
from .report import hms
from .rundir import DEFAULT_MINKNOW_DATA, find_fastq_pass
from .samplesheet import read_samplesheet, select_run


def _monitor_arguments(p):
    g = p.add_argument_group("run")
    g.add_argument("--input", required=True, help="Samplesheet (CSV/TSV): needs 'sample' and 'barcode'; optional 'run_name', 'neg_control'")
    g.add_argument("--run_name", help="Run to monitor: a folder name under --minknow_data and the run_name column. Inferred if the sheet has one run")
    g.add_argument("--run_dir", help="Exact MinKNOW run folder (or its fastq_pass folder); overrides the --run_name search")
    g.add_argument("--minknow_data", default=DEFAULT_MINKNOW_DATA, help="MinKNOW data directory (default: %(default)s)")
    g.add_argument("--outdir", help="Output directory (default: nanosentinel_<run_name>)")
    g = p.add_argument_group("references: a protocol (recommended) ...")
    g.add_argument("--protocol", help="A protocol folder made by 'prep-protocol', or the name of a bundled one (for example hcv)")
    g = p.add_argument_group("... or ordinary files (quick: only the coordinate reference is annotated)")
    g.add_argument("--possible_references", help="Multi-FASTA reference panel reads are mapped to (several similar references are fine)")
    g.add_argument("--reference_fasta", help="Reference the primer scheme and genes are drawn against (default: first panel entry)")
    g.add_argument("--primer_bed", help="Primer scheme BED (ARTIC-style <amplicon>_LEFT/_RIGHT names) on the coordinate reference")
    g.add_argument("--gff", help="GFF3 with gene/CDS features on the coordinate reference (labels the coverage plot)")
    g = p.add_argument_group("monitoring")
    g.add_argument("--target_depth", type=int, default=20, help="Depth that counts as covered (default: %(default)s)")
    g.add_argument("--target_fraction", type=float, default=0.9, help="Fraction of the target region that must reach --target_depth (default: %(default)s)")
    g.add_argument("--denominator", choices=["genome", "amplicons"], default="genome", help="Region the fraction is measured over; 'amplicons' uses the primer-scheme spans (default: %(default)s)")
    g.add_argument("--target_region", metavar="NAME", help="Measure the target over a named region of the protocol (for example E1E2) on each sample's own reference, instead of the genome (needs --protocol)")
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


def _prep_arguments(p):
    g = p.add_argument_group("inputs (ordinary FASTA / BED / GFF files)")
    g.add_argument("--possible_references", required=True, help="Multi-FASTA reference panel reads are mapped to")
    g.add_argument("--coordinate_ref", help="Id of the panel entry to use as the coordinate reference: the shared axis genes, regions and amplicons are drawn on (default: first panel entry)")
    g.add_argument("--reference_fasta", help="Same, as a FASTA file (it is added to the panel if it is not already in it)")
    g.add_argument("--gff", help="GFF3 with gene/CDS features on the coordinate reference")
    g.add_argument("--region", action="append", default=[], metavar="NAME=START-END", help="A named region on the coordinate reference, 1-based inclusive (repeatable), for example E1E2=915-2579")
    g = p.add_argument_group("primer scheme (one of)")
    g.add_argument("--primer_bed", help="Amplicons as a BED (ARTIC-style <amplicon>_LEFT/_RIGHT names) on the coordinate reference; carried onto every reference through the alignment")
    g.add_argument("--primer_fasta", help="Primer sequences (IUPAC codes allowed); their sites are searched for on every reference. Names need LEFT or RIGHT in them")
    g.add_argument("--amplicon_pattern", help="Regular expression with one group that picks the amplicon name out of a primer name (default: text before _LEFT/_RIGHT), for example '^(Seg[0-9]+)_'")
    g.add_argument("--max_mismatch", type=int, default=5, help="Most mismatches allowed between a primer and a reference site (default: %(default)s)")
    g.add_argument("--max_3prime_mismatch", type=int, default=1, help="Most of those allowed in the last 5 bases of the primer's 3' end (default: %(default)s)")
    g = p.add_argument_group("alignment and output")
    g.add_argument("--panel_alignment", help="A ready-made alignment of the panel (including the coordinate reference), instead of running MAFFT")
    g.add_argument("--group_pattern", default=DEFAULT_GROUP_PATTERN, help="Regular expression whose first group names a reference's sibling group, such as its subtype (default: %(default)s)")
    g.add_argument("--name", help="Protocol name (default: the folder name)")
    g.add_argument("--threads", type=int, default=4, help="Threads for MAFFT (default: %(default)s)")
    g.add_argument("--outdir", required=True, help="Protocol folder to write")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="nanosentinel",
        description="Monitor a live ONT run from a samplesheet: map reads as MinKNOW writes them, show coverage on a local "
        "dashboard, record how long each sample takes to reach a target depth, and (optionally) concatenate reads.",
    )
    p.add_argument("--version", action="version", version=f"nanosentinel {__version__}")
    sub = p.add_subparsers(dest="command", metavar="{prep-protocol,monitor}", required=True)
    prep = sub.add_parser("prep-protocol", help="Prepare a reference panel once: align it, annotate every reference, QC it",
                          description="Align every reference in a panel to a coordinate reference and lift genes, regions and "
                          "amplicons onto each, writing a protocol folder that 'monitor --protocol' loads in about a second.")  # fmt: skip
    _prep_arguments(prep)
    mon = sub.add_parser("monitor", help="Monitor a run: map reads, show the dashboard, record time to target, concatenate reads",
                         description="Monitor a live ONT run from a samplesheet.")  # fmt: skip
    _monitor_arguments(mon)
    args = p.parse_args(argv)
    if args.command == "prep-protocol":
        if args.primer_bed and args.primer_fasta:
            p.error("give the primer scheme as --primer_bed or --primer_fasta, not both")
        if args.coordinate_ref and args.reference_fasta:
            p.error("give the coordinate reference as --coordinate_ref or --reference_fasta, not both")
        if args.amplicon_pattern and not args.primer_fasta:
            p.error("--amplicon_pattern applies to --primer_fasta")
        return args
    if args.skip_monitor and not args.concatenate_reads:
        p.error("--skip_monitor does nothing unless --concatenate_reads is also given")
    files = args.possible_references or args.reference_fasta
    if args.protocol and (files or args.primer_bed or args.gff):
        p.error("--protocol replaces --possible_references, --reference_fasta, --primer_bed and --gff; give one or the other")
    if not args.skip_monitor and not (args.protocol or files):
        p.error("monitoring needs --protocol, or --possible_references and/or --reference_fasta")
    if args.target_region:
        if args.denominator != "genome":
            p.error("--target_region replaces --denominator; give one or the other")
        if not args.protocol:
            p.error("--target_region needs --protocol (regions are defined when a protocol is built)")
        args.denominator = f"region:{args.target_region}"
    if args.denominator == "amplicons" and not (args.protocol or args.primer_bed):
        p.error("--denominator amplicons needs a primer scheme: --protocol, or --primer_bed")
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
                console.say(f"Waiting for MinKNOW to create the run's fastq_pass folder ({exc}). Ctrl+C to cancel.", style="warn")
                announced = True
            console.status([("waiting for fastq_pass", hms(time.time() - started), "warn")])
            time.sleep(5)


def run(args):
    header, rows = read_samplesheet(args.input)
    run_name, rows = select_run(rows, args.run_name)
    outdir = Path(args.outdir or f"nanosentinel_{run_name or 'run'}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    console = Console(outdir / "nanosentinel.log")

    n_ctrl = sum(r.neg_control for r in rows)
    console.say(f"nanosentinel {__version__} | run {run_name or '(unnamed)'} | {len(rows) - n_ctrl} sample(s), {n_ctrl} control(s)", style="title")
    fastq_pass = locate_reads(args, run_name, console, wait=not (args.skip_monitor or args.dry_run))
    console.end_status()
    console.kv("Reads", fastq_pass)
    console.kv("Output", outdir)

    if not args.skip_monitor:
        from .monitor import monitor  # imported late so --help and --skip_monitor stay fast

        if args.protocol:
            refinfo = load_protocol(args.protocol, outdir)
            if args.denominator.startswith("region:"):
                check_region(refinfo, args.denominator[len("region:"):])  # fail now, not after a wait for MinKNOW
            annotated = len(refinfo.maps)
            console.say(f"Protocol {refinfo.protocol}: {len(refinfo.lengths)} reference(s), {annotated} aligned to {refinfo.coordinate_ref}")
        else:
            refinfo = prepare_references(outdir, args.possible_references, args.reference_fasta, args.primer_bed, args.gff)
            console.say(f"Reference panel: {len(refinfo.lengths)} sequence(s); primer scheme and genes are drawn against {refinfo.coordinate_ref}")
            if len(refinfo.lengths) > 1:
                console.say("Quick mode annotates the coordinate reference only; 'prep-protocol' annotates every reference.", style="muted")
        if args.dry_run:
            console.say("Dry run: inputs are valid and the reference index was built; nothing was monitored.", style="muted")
            return 0
        monitor(
            console, rows, refinfo, fastq_pass, outdir, run_name, args.host, args.port, args.target_depth,
            args.target_fraction, args.denominator, args.poll_interval, args.webhook, not args.no_browser,
            args.clear_annotations, args.min_identity, args.threads,
        )  # fmt: skip
    elif args.dry_run:
        console.say("Dry run: inputs are valid; nothing was done.", style="muted")
        return 0

    if args.concatenate_reads:
        console.say("Concatenating reads ...", style="muted")
        sheet, n = concatenate(rows, header, fastq_pass, outdir, run_name, args.gzip, args.overwrite, log=console.warn)
        console.say(f"Wrote {n} FASTQ file(s) to {outdir / 'reads'}", style="ok")
        console.kv("Samplesheet with read paths", sheet)
    return 0


def prep(args):
    outdir = Path(args.outdir).resolve()
    console = Console()
    regions = dict(parse_region(r) for r in args.region)
    qc = build_protocol(
        outdir, args.possible_references, args.reference_fasta, args.gff, regions, args.primer_bed, args.primer_fasta,
        args.amplicon_pattern, args.panel_alignment, args.name, args.group_pattern, args.threads, args.max_mismatch,
        say=lambda line: console.say(line, style="muted"), coordinate_ref=args.coordinate_ref,
        max_tail_mismatch=args.max_3prime_mismatch,
    )  # fmt: skip
    console.kv("Protocol", outdir, style="url")
    report = (outdir / "qc.txt").read_text().rstrip("\n").splitlines()
    for line in report:
        bad = ("poorly aligned" in line or "too few close matches" in line or " absent" in line or "odd_span" in line
               or "(no site on" in line)  # fmt: skip
        console.say(line, style="warn" if bad else None)
    if qc["poorly_aligned"] or qc["unplaced_primers"]:
        console.warn(f"See {outdir / 'qc.txt'} for the references and primers that need a look.")
    console.say(f"Use it with: nanosentinel monitor --protocol {outdir} ...", style="ok")
    return 0


def main(argv=None):
    args = parse_args(argv)
    try:
        return prep(args) if args.command == "prep-protocol" else run(args)
    except SentinelError as exc:
        error(exc)
        return 1
    except KeyboardInterrupt:
        print(file=sys.stderr)
        notice("Interrupted.")
        return 130
