"""Per-sample depth tracking and 'time to target' (e.g. >=20x over >=90% of the genome)."""

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np

from .annotate import ANNOTATION_FORMAT, natural_key


LOCK_MIN_READS = 200  # mapped reads to "taste" before a sibling group is chosen
LOCK_SWITCH_RATIO = 1.5  # another group must lead the locked one by this factor ...
LOCK_SWITCH_MARGIN = 100  # ... and by this many reads before the lock moves


class SampleTracker:
    """Tracks, per reference, depth along the genome and how much of it is >= target_depth.

    Reads are expected in (roughly) time order; the time recorded for the moment the target is reached is the
    latest read start time seen so far, which keeps it monotonic. Every alignment adds coverage to its reference,
    but only a read's primary alignment counts as a read (read totals, length statistics, composition).
    """

    def __init__(self, ref_lengths, target_depth=20, target_fraction=0.9, masks=None, group_of=None, denominators=None):
        if not 1 <= target_depth <= 255:
            raise ValueError("target_depth must be between 1 and 255")
        if not 0 < target_fraction <= 1:
            raise ValueError("target_fraction must be in (0, 1]")
        self.ref_lengths = ref_lengths
        self.target = int(target_depth)
        self.fraction = target_fraction
        self.masks = masks or {}
        self.denominators = denominators or {}  # reference -> size of the target, when it is not the mask's own size
        self.group_of = group_of or (lambda ref: ref)
        self.locked_group = None
        self.depth, self.covered = {}, {}
        self.reads = self.mapped = self.bases = 0
        self.ref_reads, self.ref_bases = Counter(), Counter()  # per-reference primary reads / aligned bases
        self.len_counts = Counter()  # read length -> number of reads (all reads, mapped or not)
        self.series = []  # (time, reads, mapped, bases, best_reference, best_fraction), one per annotation file
        self.first = self.latest = None
        self.reached = False
        self.reached_at = self.reached_ref = None
        self.reads_at_target = self.bases_at_target = 0

    def denominator(self, ref):
        mask = self.masks.get(ref)
        if mask is not None:
            # A region target is measured over the region's length on the coordinate axis, so a reference that lacks
            # part of the region (a partial sequence, a deletion) is held to the whole of it.
            return self.denominators.get(ref, int(mask.sum()))
        return 0 if self.masks else self.ref_lengths[ref]  # amplicon mode: a reference without amplicons cannot count

    def add(self, ref, start, end, when=None, read_len=0, primary=True):
        """Add one alignment row. Returns True the first time the target is reached."""
        if primary:
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
        if primary:
            self.mapped += 1
            self.bases += end - start
            self.ref_reads[ref] += 1
            self.ref_bases[ref] += end - start
        depth = self.depth.get(ref)
        if depth is None:
            depth = self.depth[ref] = np.zeros(length, dtype=np.uint32)
            self.covered[ref] = 0
        seg = depth[start:end]
        arriving = seg == self.target - 1
        mask = self.masks.get(ref)
        if mask is not None:
            arriving &= mask[start:end]
        self.covered[ref] += int(np.count_nonzero(arriving))
        seg += 1  # in place on the view
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
        self.update_lock()
        ref, frac = self.best()
        self.series.append((self.latest, self.reads, self.mapped, self.bases, ref, frac))

    def group_reads(self):
        """Primary reads per sibling group (for example per subtype)."""
        totals = Counter()
        for ref, n in self.ref_reads.items():
            totals[self.group_of(ref)] += n
        return totals

    def update_lock(self):
        """Taste, then lock: choose the sibling group after the first reads, and move only on a clear, large lead.

        Near-identical references share reads, so the single most-read reference flips between siblings; the group
        (for example the subtype) does not. The lock keeps the dashboard from jumping between them.
        """
        totals = self.group_reads()
        if sum(totals.values()) < LOCK_MIN_READS:
            return
        top = max(totals, key=totals.get)
        if self.locked_group is None:
            self.locked_group = top
        elif top != self.locked_group:
            now, lead = totals[self.locked_group], totals[top]
            if lead >= LOCK_SWITCH_RATIO * now and lead - now >= LOCK_SWITCH_MARGIN:
                self.locked_group = top

    def shown(self):
        """The reference to display: the best-covered sibling in the locked group (before the lock: the best overall)."""
        if self.locked_group is not None:
            members = [ref for ref in self.covered if self.group_of(ref) == self.locked_group]
            if members:
                return max(members, key=lambda ref: (self.fraction_of(ref), self.ref_reads[ref]))
        return self.best()[0]

    def dominant(self):
        """Reference with the most primary reads, or None."""
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

    def mean_depth(self, ref):
        depth = self.depth.get(ref)
        return float(depth.mean()) if depth is not None else 0.0

    def binned_depth(self, ref, nbins=300):
        """Mean depth in ``nbins`` equal windows along ``ref`` (None if the reference has no coverage)."""
        depth = self.depth.get(ref)
        if depth is None:
            return None
        nbins = min(nbins, len(depth))
        edges = np.linspace(0, len(depth), nbins + 1).astype(int)
        sums = np.add.reduceat(depth.astype(np.float64), edges[:-1])
        return (sums / np.diff(edges)).round(1).tolist()



def parse_time(text):
    try:
        return datetime.fromisoformat(text.strip())
    except (ValueError, AttributeError):
        return None


class AnnotationReader:
    """Feeds the annotation files (annotations/<barcode>/*.tsv) to the trackers as they appear."""

    def __init__(self, annotation_dir, trackers, panel_hash):
        self.dir = Path(annotation_dir)
        self.trackers = trackers  # barcode -> SampleTracker
        self.panel_hash = panel_hash
        self.done = set()
        self.run_start = None
        self.changed = False  # new data since the reports were last written

    def poll(self):
        """Read any new annotation files. Returns the barcodes that reached target during this poll."""
        newly = []
        for barcode, tracker in self.trackers.items():
            folder = self.dir / barcode
            if not folder.is_dir():
                continue
            for path in sorted(folder.glob("*.tsv"), key=natural_key):
                if path in self.done:
                    continue
                result = self._read_file(path, tracker)
                if result is None:
                    continue  # made against a different panel: the annotator will redo it
                self.done.add(path)
                self.changed = True
                if result and barcode not in newly:
                    newly.append(barcode)
            if tracker.first and (self.run_start is None or tracker.first < self.run_start):
                self.run_start = tracker.first
        return newly

    def _read_file(self, path, tracker):
        reached = False
        with open(path, newline="") as fh:
            if f"panel={self.panel_hash} format={ANNOTATION_FORMAT}" not in fh.readline():
                return None
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    ref = row["ref"]
                    start, end = int(row["start"]), int(row["end"])
                    read_len = int(row["read_len"])
                    primary = row["primary"] == "1"
                except (KeyError, ValueError, TypeError):
                    continue  # malformed line
                if tracker.add(ref, start, end, parse_time(row.get("start_time", "")), read_len, primary):
                    reached = True
        tracker.snapshot()
        return reached
