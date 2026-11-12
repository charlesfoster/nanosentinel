"""Replay a finished run's fastq_pass folder into a new folder as if it were being sequenced right now.

    python dev/replay_run.py SRC_fastq_pass DEST_fastq_pass [--speed N] [--dry_run]

A monitor watching DEST sees the chunks arrive in the order they were produced. Start the monitor first (it waits for the
folder to appear), then the replay:

    nanosentinel monitor --run_dir DEST_fastq_pass --protocol hcv --input samplesheet.csv ...
    python dev/replay_run.py SRC_fastq_pass DEST_fastq_pass --speed 60      # one hour of sequencing in a minute

How it behaves
  * Order and pacing come from the reads' own start times (`start_time=` or the `st:Z:` tag), never from file modification times (which are
    lost on copy or rsync). A chunk is released when its last read was produced, because that is when MinKNOW would have
    finished writing it. The first read breaks ties.
  * --speed N replays at N times real time (default 1). --speed 0 copies everything at once.
  * Each chunk is copied under a name the monitor ignores, then renamed, so a half-copied chunk is never visible.
  * Every barcode folder (and unclassified) is replayed; files that are not FASTQ are ignored.
  * DEST must not already contain FASTQ files.

Results that depend only on read timestamps (such as time to target) do not depend on the replay speed.
Standard library only, so any python3 will do.
"""

import argparse
import gzip
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

FASTQ_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
_START_TIME = re.compile(r"(?:start_time=|(?<=\s)st:Z:)(\S+)")  # older "start_time=", or the SAM-style "st:Z:" tag


@dataclass
class Chunk:
    path: Path
    rel: Path  # barcode folder / file name, as it will appear under DEST
    first: datetime = None  # start time of the first read
    last: datetime = None  # start time of the last read
    reads: int = 0
    offset: float = 0.0  # seconds from the start of the run until this chunk is complete


def natural_key(name):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(name))]


def parse_time(text):
    text = text.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        when = datetime.fromisoformat(text)
    except ValueError:
        return None
    return when if when.tzinfo else when.replace(tzinfo=timezone.utc)


def read_times(path):
    """(first start time, last start time, number of reads) of a FASTQ; times are None if it has none."""
    opener = gzip.open if str(path).endswith(".gz") else open
    first = last = None
    reads = 0
    with opener(path, "rt") as fh:
        for number, line in enumerate(fh):
            if number % 4:
                continue
            reads += 1
            match = _START_TIME.search(line)
            when = parse_time(match.group(1)) if match else None
            if when is not None:
                first = first or when
                last = when if last is None or when > last else last
    return first, last, reads


def find_chunks(src):
    """Every FASTQ chunk under the barcode folders of ``src``, unordered."""
    src = Path(src)
    if not src.is_dir():
        raise SystemExit(f"{src} is not a folder")
    chunks = []
    for folder in sorted(p for p in src.iterdir() if p.is_dir()):
        for path in sorted(folder.iterdir(), key=lambda p: natural_key(p.name)):
            if path.name.endswith(FASTQ_SUFFIXES) and not path.name.startswith("."):
                chunks.append(Chunk(path, Path(folder.name) / path.name))
    if not chunks:
        raise SystemExit(f"no FASTQ files in the barcode folders under {src}")
    return chunks


def build_schedule(chunks, warn=print):
    """Fill in each chunk's times and release offset, and return the chunks in release order."""
    for chunk in chunks:
        try:
            chunk.first, chunk.last, chunk.reads = read_times(chunk.path)
        except (OSError, EOFError) as exc:
            raise SystemExit(f"cannot read {chunk.path}: {exc}")
    # A chunk without timestamps inherits the previous chunk's time in the same barcode folder (name order).
    previous = {}
    for chunk in sorted(chunks, key=lambda c: (c.rel.parent.name, natural_key(c.rel.name))):
        folder = chunk.rel.parent.name
        if chunk.last is None:
            if folder in previous:
                chunk.first = chunk.last = previous[folder]
                warn(f"warning: {chunk.rel} has no start_time in its reads; placing it after the previous chunk of {folder}")
            else:
                continue
        previous[folder] = chunk.last
    timed = [c for c in chunks if c.last is not None]
    if not timed:
        raise SystemExit("none of the reads has a start_time= or st:Z: field, so the run cannot be ordered or paced")
    for chunk in chunks:
        if chunk.last is None:  # no timestamp anywhere before it in its folder: release it first
            chunk.first = chunk.last = min(c.first for c in timed)
            warn(f"warning: {chunk.rel} has no start_time in its reads; releasing it first")
    t0 = min(c.first for c in chunks)
    for chunk in chunks:
        chunk.offset = (chunk.last - t0).total_seconds()
    return sorted(chunks, key=lambda c: (c.last, c.first, natural_key(c.rel)))


def copy_atomic(src, dest):
    """Copy under a hidden name the monitor does not treat as a FASTQ, then rename into place."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.parent / f".{dest.name}.partial"
    try:
        shutil.copyfile(src, partial)
        os.replace(partial, dest)
    finally:
        if partial.exists():
            partial.unlink()


def hms(seconds):
    seconds = int(max(0, seconds))
    return f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def replay(chunks, dest, speed=1.0, clock=time.monotonic, sleep=time.sleep, say=print):
    """Copy the chunks (already in release order) into ``dest``, pacing them by their offsets / speed."""
    started = clock()
    for number, chunk in enumerate(chunks, 1):
        if speed > 0:
            wait = chunk.offset / speed - (clock() - started)
            if wait > 0:
                sleep(wait)
        copy_atomic(chunk.path, Path(dest) / chunk.rel)
        say(f"[run time {hms(chunk.offset)}] {chunk.rel}  ({number}/{len(chunks)}, {chunk.reads:,} reads)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="fastq_pass folder of a finished run")
    ap.add_argument("dest", help="folder to replay into (created; must not already hold FASTQ files)")
    ap.add_argument("--speed", type=float, default=1.0, help="times real time; 0 copies everything at once (default: 1)")
    ap.add_argument("--dry_run", action="store_true", help="show the schedule without copying anything")
    args = ap.parse_args(argv)
    if args.speed < 0:
        ap.error("--speed cannot be negative")
    dest = Path(args.dest)
    if dest.resolve() == Path(args.src).resolve():
        ap.error("DEST must be a different folder from SRC")
    if dest.is_dir() and any(p.name.endswith(FASTQ_SUFFIXES) for p in dest.rglob("*") if p.is_file()):
        ap.error(f"{dest} already contains FASTQ files; replay into a fresh folder")

    chunks = build_schedule(find_chunks(args.src), warn=lambda m: print(m, file=sys.stderr))
    span = chunks[-1].offset
    wall = f"{hms(span / args.speed)} at {args.speed:g}x" if args.speed > 0 else "all at once"
    print(f"{len(chunks)} chunk(s), {sum(c.reads for c in chunks):,} reads, run time {hms(span)}; replaying {wall}")
    if args.dry_run:
        for c in chunks:
            print(f"  {hms(c.offset)}  {c.rel}  ({c.reads:,} reads)")
        return 0
    dest.mkdir(parents=True, exist_ok=True)
    try:
        replay(chunks, dest, args.speed)
    except KeyboardInterrupt:
        print("\nInterrupted; the chunks copied so far are complete.", file=sys.stderr)
        return 130
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
