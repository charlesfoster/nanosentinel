"""Concatenate a run's per-barcode FASTQ chunks into one file per sample and write a samplesheet pointing at them."""

import csv
import gzip
import re
from pathlib import Path

from .errors import ConcatError


def _natural_key(path):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", Path(path).name)]


def barcode_files(fastq_pass, barcode):
    folder = Path(fastq_pass) / barcode
    if not folder.is_dir():
        return []
    files = [p for p in folder.iterdir() if p.name.endswith((".fastq.gz", ".fastq", ".fq.gz", ".fq"))]
    return sorted(files, key=_natural_key)


def _copy_binary(files, out, plain_text):
    last = b"\n"
    for f in files:
        with open(f, "rb") as src:
            while chunk := src.read(1 << 20):
                out.write(chunk)
                last = chunk[-1:]
        if plain_text and last != b"\n":  # keep records separate if a file lacks a trailing newline
            out.write(b"\n")
            last = b"\n"


def concatenate(rows, header, fastq_pass, outdir, run_name, gzip_output=False, overwrite=False, log=print):
    """Write <outdir>/reads/<sample>_<run>_<barcode>.fastq[.gz] for each row plus samplesheet_with_reads.csv."""
    outdir = Path(outdir)
    reads_dir = outdir / "reads"
    plans = []
    for r in rows:
        files = barcode_files(fastq_pass, r.barcode)
        if not files:
            log(f"No reads found for {r.sample} ({r.barcode}); skipping")
            continue
        gz = [f.name.endswith(".gz") for f in files]
        if any(gz) and not all(gz):
            raise ConcatError(f"{r.sample} ({r.barcode}) has a mix of gzipped and plain FASTQ files; cannot concatenate safely")
        compress = all(gz) or gzip_output
        stem = "_".join(x for x in (r.sample, run_name, r.barcode) if x)
        out = reads_dir / f"{stem}.fastq{'.gz' if compress else ''}"
        plans.append((r, files, all(gz), compress, out))
    if not plans:
        raise ConcatError(f"No FASTQ files found for any sample under {fastq_pass}")
    exists = [p[4] for p in plans if p[4].exists()]
    if exists and not overwrite:
        raise ConcatError(f"{len(exists)} output file(s) already exist (e.g. {exists[0]}); use --overwrite to replace them")

    reads_dir.mkdir(parents=True, exist_ok=True)
    sheet_rows = []
    for r, files, inputs_gz, compress, out in plans:
        tmp = out.with_name(out.name + ".tmp")
        if inputs_gz or not compress:
            with open(tmp, "wb") as fh:  # gzip members can be concatenated byte-for-byte
                _copy_binary(files, fh, plain_text=not inputs_gz)
        else:
            with gzip.open(tmp, "wb") as fh:
                _copy_binary(files, fh, plain_text=not inputs_gz)
        tmp.replace(out)
        sheet_rows.append({**r.fields, "long_reads": str(out.resolve())})

    columns = list(header)
    if "long_reads" not in columns:
        columns.insert(columns.index("sample") + 1, "long_reads")
    sheet = outdir / "samplesheet_with_reads.csv"
    with open(sheet, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sheet_rows)
    return sheet, len(plans)
