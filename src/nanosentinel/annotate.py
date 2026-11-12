"""Read annotation: map each finished FASTQ chunk against the reference panel.

Every finished FASTQ chunk MinKNOW writes under fastq_pass/<barcode>/ is mapped with minimap2 against the reference
panel and turned into one small annotation file, annotations/<barcode>/<chunk>.tsv. That file is the only interface
between this module and the coverage tracking: one row per alignment, with read length and start time taken from the
FASTQ itself. Files record the source size and the panel they were made against, so a restarted monitor reuses
valid annotations and redoes stale ones.

Near-identical references: minimap2 also reports secondary alignments that score within ~10% of the primary. They are
kept (primary=0) and credited to coverage, so a read that fits three sibling genotypes equally well counts towards all
three instead of being split between them; read counts and composition use the primary alignment only.
"""

import gzip
import os
import re
import subprocess
import threading
import time
from pathlib import Path

from . import tools

FASTQ_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
COLUMNS = ["read_name", "read_len", "start_time", "ref", "ref_len", "start", "end", "mapq", "identity", "primary"]
# The read's start time: "start_time=..." (older MinKNOW and guppy) or the SAM-style tag "st:Z:..." (current MinKNOW/dorado).
_START_TIME = re.compile(r"(?:start_time=|(?<=\s)st:Z:)(\S+)")
# Bump when what an annotation file holds changes, so annotations made the old way are redone instead of reused.
ANNOTATION_FORMAT = 2


def natural_key(path):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", Path(path).name)]


def fastq_files(folder):
    folder = Path(folder)
    if not folder.is_dir():
        return []
    return sorted((p for p in folder.iterdir() if p.name.endswith(FASTQ_SUFFIXES)), key=natural_key)


def read_fastq_info(path):
    """(read_name, length, start_time) for every read, in file order. Raises on truncated / unreadable files."""
    opener = gzip.open if str(path).endswith(".gz") else open
    reads = []
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            fh.readline()
            qual = fh.readline()
            if not header.startswith("@") or not qual:
                raise ValueError(f"truncated or malformed FASTQ: {path}")
            match = _START_TIME.search(header)
            reads.append((header[1:].split()[0], len(seq.strip()), match.group(1) if match else ""))
    return reads


def parse_paf(lines, min_identity=0.8):
    """PAF lines -> {read_name: [row, ...]} where the first row is the read's primary alignment.

    A row is (ref, ref_len, start, end, mapq, identity, primary). Reads whose primary alignment is below
    ``min_identity`` get no rows (they are unmapped).
    """
    grouped = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or f[5] == "*":
            continue
        block = int(f[10])
        identity = int(f[9]) / block if block else 0.0
        tags = {t[:2]: t[5:] for t in f[12:]}
        grouped.setdefault(f[0], []).append(
            {
                "row": (f[5], int(f[6]), int(f[7]), int(f[8]), int(f[11]), identity),
                "tp": tags.get("tp", "P"),
                "score": int(tags.get("AS", 0)),
            }
        )
    result = {}
    for name, recs in grouped.items():
        primaries = [r for r in recs if r["tp"] in ("P", "p")]
        if not primaries:
            continue
        best = max(primaries, key=lambda r: r["score"])
        if best["row"][5] < min_identity:
            continue
        rows = [(*best["row"], 1)]
        rows += [(*r["row"], 0) for r in recs if r is not best and r["row"][5] >= min_identity]
        result[name] = rows
    return result


def write_annotation(path, source_size, panel_hash, reads, alignments):
    """Atomically write one chunk's annotation file."""
    tmp = Path(f"{path}.tmp")
    with open(tmp, "w") as fh:
        fh.write(f"#source_size={source_size} panel={panel_hash} format={ANNOTATION_FORMAT}\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for name, length, start_time in reads:
            for ref, ref_len, start, end, mapq, identity, primary in alignments.get(name) or [("*", 0, 0, 0, 0, 0.0, 1)]:
                fh.write(f"{name}\t{length}\t{start_time}\t{ref}\t{ref_len}\t{start}\t{end}\t{mapq}\t{identity:.3f}\t{primary}\n")
    os.replace(tmp, path)


def annotation_is_current(path, source_size, panel_hash):
    try:
        with open(path) as fh:
            first = fh.readline()
    except OSError:
        return False
    return first.strip() == f"#source_size={source_size} panel={panel_hash} format={ANNOTATION_FORMAT}"


class Annotator(threading.Thread):
    """Background thread: watch the barcode folders and annotate each finished chunk."""

    def __init__(self, fastq_pass, barcodes, annotation_dir, refinfo, min_identity=0.8, threads=4, log=None,
                 poll=2.0, stable_seconds=4.0):  # fmt: skip
        super().__init__(daemon=True)
        self.fastq_pass, self.barcodes = Path(fastq_pass), list(barcodes)
        self.annotation_dir, self.refinfo = Path(annotation_dir), refinfo
        self.min_identity, self.threads = min_identity, threads
        self.log, self.poll, self.stable = log, poll, stable_seconds
        self.minimap2 = tools.find("minimap2")
        self._stop_event = threading.Event()
        self._proc = None
        self._sizes = {}  # path -> size seen on the previous pass (a chunk is finished once it stops changing)
        self._done = set()
        self._failures = {}
        self.errors = []  # human-readable problems, surfaced by the monitor
        self.files_annotated = 0

    def stop(self):
        self._stop_event.set()
        proc = self._proc
        if proc is not None and proc.poll() is None:
            proc.terminate()

    def run(self):
        while not self._stop_event.is_set():
            worked = False
            for barcode in self.barcodes:
                for path in fastq_files(self.fastq_pass / barcode):
                    if self._stop_event.is_set():
                        return
                    if path in self._done:
                        continue
                    if self._handle(barcode, path):
                        worked = True
                        break  # one chunk per barcode per pass, so a busy barcode cannot starve the others
            if not worked:
                self._stop_event.wait(self.poll)

    def _handle(self, barcode, path):
        """Annotate ``path`` if it is finished. Returns True if work was done (or attempted)."""
        try:
            stat = path.stat()
        except OSError:
            return False
        out = self.annotation_dir / barcode / (re.sub(r"\.(fastq|fq)(\.gz)?$", "", path.name) + ".tsv")
        if annotation_is_current(out, stat.st_size, self.refinfo.panel_hash):
            self._done.add(path)
            return False
        previous, self._sizes[path] = self._sizes.get(path), stat.st_size
        if previous != stat.st_size or time.time() - stat.st_mtime < self.stable:
            return False  # still being written
        try:
            reads = read_fastq_info(path)
            alignments = self._map(path)
            if alignments is None:  # stopped mid-way
                return True
            out.parent.mkdir(parents=True, exist_ok=True)
            write_annotation(out, stat.st_size, self.refinfo.panel_hash, reads, alignments)
            self._done.add(path)
            self.files_annotated += 1
        except (OSError, EOFError, ValueError, subprocess.SubprocessError) as exc:
            self._failures[path] = self._failures.get(path, 0) + 1
            self._sizes.pop(path, None)  # retry after it has been seen stable again
            if self._failures[path] >= 3:
                self._done.add(path)
                self.errors.append(f"could not annotate {path.name}: {exc}")
                if self.log:
                    self.log.warning(self.errors[-1])
        return True

    def _map(self, path):
        cmd = [self.minimap2, "-v", "1", "-x", "map-ont", "-c", "--secondary=yes", "-N", "20", "-p", "0.9",
               "-t", str(self.threads), str(self.refinfo.index), str(path)]  # fmt: skip
        self._proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out, err = self._proc.communicate()
        code, self._proc = self._proc.returncode, None
        if self._stop_event.is_set():
            return None
        if code != 0:
            raise subprocess.SubprocessError(f"minimap2 exited with {code}: {err.strip()[-200:]}")
        return parse_paf(out.splitlines(), self.min_identity)
