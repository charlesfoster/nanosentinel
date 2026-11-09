"""Per-sample depth tracking and 'time to target' (e.g. >=20x over >=90% of the genome)."""

import csv
import re
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np


class SampleTracker:
    """Tracks, per reference, how much of the (target region of the) genome is >= target_depth.

    Depth is capped at ``target_depth`` (uint8), which is all that is needed to know when a position
    first reaches the target. Reads are expected in (roughly) time order; the time recorded for the
    moment of reaching target is the latest read start time seen so far, which keeps it monotonic.
    """

    def __init__(self, ref_lengths, target_depth=20, target_fraction=0.9, masks=None):
        if not 1 <= target_depth <= 255:
            raise ValueError("target_depth must be between 1 and 255")
        if not 0 < target_fraction <= 1:
            raise ValueError("target_fraction must be in (0, 1]")
        self.ref_lengths = ref_lengths
        self.target = int(target_depth)
        self.fraction = target_fraction
        self.masks = masks or {}
        self.depth, self.covered = {}, {}
        self.reads = self.mapped = self.bases = 0
        self.ref_reads, self.ref_bases = Counter(), Counter()  # per-reference mapped reads / aligned bases
        self.len_counts = Counter()  # read length -> number of reads (all reads, mapped or not)
        self.series = []  # (time, reads, mapped, bases, best_reference, best_fraction), one per annotation file
        self.first = self.latest = None
        self.reached = False
        self.reached_at = self.reached_ref = None
        self.reads_at_target = self.bases_at_target = 0

    def denominator(self, ref):
        mask = self.masks.get(ref)
        return int(mask.sum()) if mask is not None else self.ref_lengths[ref]

    def add(self, ref, start, end, when=None, read_len=0):
        """Add one annotated read. Returns True the first time the target is reached."""
        self.reads += 1
        if read_len > 0:
            self.len_counts[read_len] += 1
        if when is not None:
            if self.latest is None or when > self.latest:
                self.latest = when
            if self.first is None or when < self.first:
                self.first = when
        length = self.ref_lengths.get(ref)
        if length is None:
            return False
        start, end = max(0, start), min(end, length)
        if end <= start:
            return False
        self.mapped += 1
        self.bases += end - start
        self.ref_reads[ref] += 1
        self.ref_bases[ref] += end - start
        depth = self.depth.get(ref)
        if depth is None:
            depth = self.depth[ref] = np.zeros(length, dtype=np.uint8)
            self.covered[ref] = 0
        seg = depth[start:end]
        arriving = seg == self.target - 1
        mask = self.masks.get(ref)
        if mask is not None:
            arriving &= mask[start:end]
        self.covered[ref] += int(np.count_nonzero(arriving))
        seg += seg < self.target  # in-place on the view: cap at target
        if not self.reached and self.fraction_of(ref) >= self.fraction - 1e-12:
            self.reached = True
            self.reached_at, self.reached_ref = self.latest, ref
            self.reads_at_target, self.bases_at_target = self.reads, self.bases
            return True
        return False

    def snapshot(self):
        """Record progress (called after each annotation file, so the series follows the data, not the wall clock)."""
        if self.reads == 0 or (self.series and self.series[-1][1] == self.reads):
            return
        ref, frac = self.best()
        self.series.append((self.latest, self.reads, self.mapped, self.bases, ref, frac))

    def dominant(self):
        """Reference with the most mapped reads, or None."""
        return max(self.ref_reads, key=self.ref_reads.get) if self.ref_reads else None

    def length_stats(self):
        """Read length summary (n, total, mean, median, n50, min, max) or None if no reads."""
        n = sum(self.len_counts.values())
        if not n:
            return None
        lengths = sorted(self.len_counts)
        total = sum(length * count for length, count in self.len_counts.items())

        def kth(k):  # k-th smallest length, 0-based
            seen = 0
            for length in lengths:
                seen += self.len_counts[length]
                if seen > k:
                    return length

        median = (kth((n - 1) // 2) + kth(n // 2)) / 2
        seen, n50 = 0, lengths[-1]
        for length in reversed(lengths):
            seen += length * self.len_counts[length]
            if seen >= total / 2:
                n50 = length
                break
        return {"n": n, "total": total, "mean": total / n, "median": median, "n50": n50, "min": lengths[0], "max": lengths[-1]}

    def fraction_of(self, ref):
        denom = self.denominator(ref)
        return self.covered.get(ref, 0) / denom if denom else 0.0

    def best(self):
        """(reference, fraction) with the highest fraction at target depth, or (None, 0.0)."""
        if not self.covered:
            return None, 0.0
        ref = max(self.covered, key=self.fraction_of)
        return ref, self.fraction_of(ref)


def _natural_key(path):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", Path(path).name)]


def parse_time(text):
    try:
        return datetime.fromisoformat(text.strip())
    except (ValueError, AttributeError):
        return None


class AnnotationReader:
    """Incrementally feeds RAMPART's per-file annotation CSVs (annotations/<barcode>/*.csv) to trackers."""

    def __init__(self, annotation_dir, trackers, min_age=2.0):
        self.dir = Path(annotation_dir)
        self.trackers = trackers  # barcode -> SampleTracker
        self.min_age = min_age
        self.done = set()
        self.run_start = None
        self.changed = False  # new data since the reports were last written

    def poll(self, final=False, now=None):
        """Read any new, finished CSVs. Returns the barcodes that reached target during this poll."""
        now = time.time() if now is None else now
        newly = []
        for barcode, tracker in self.trackers.items():
            folder = self.dir / barcode
            if not folder.is_dir():
                continue
            for path in sorted(folder.glob("*.csv"), key=_natural_key):
                if path in self.done:
                    continue
                if not final and now - path.stat().st_mtime < self.min_age:
                    break  # still being written; later files in order are even newer
                self.done.add(path)
                self.changed = True
                if self._read_file(path, tracker) and barcode not in newly:
                    newly.append(barcode)
            if tracker.first and (self.run_start is None or tracker.first < self.run_start):
                self.run_start = tracker.first
        return newly

    @staticmethod
    def _read_file(path, tracker):
        reached = False
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    ref = row["best_reference"]
                    start, end = int(row["start_coords"]), int(row["end_coords"])
                    read_len = int(row["read_len"])
                except (KeyError, ValueError, TypeError):
                    continue  # truncated / malformed line
                if tracker.add(ref, start, end, parse_time(row.get("start_time", "")), read_len):
                    reached = True
        tracker.snapshot()
        return reached
