"""Tab-separated reports written during and at the end of monitoring (all written atomically)."""

import os
from datetime import timezone


def hms(seconds):
    seconds = int(max(0, seconds))
    return f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def _write(path, header, records):
    tmp = f"{path}.tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for rec in records:
            fh.write("\t".join(str(c) for c in rec) + "\n")
    os.replace(tmp, path)


def _utc(when):
    return when.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") if when else "NA"


def _minutes(when, run_start):
    return f"{(when - run_start).total_seconds() / 60:.2f}" if when and run_start else "NA"


def _bool(value):
    return "true" if value else "false"


TIMING_HEADER = [
    "sample", "barcode", "neg_control", "status", "time_to_target_min", "time_to_target_hms", "reached_at_utc",
    "reads_at_target", "mapped_bases_at_target", "best_reference", "best_fraction", "reads_total",
    "target_depth", "target_fraction", "denominator",
]  # fmt: skip


def write_timing(path, rows, trackers, run_start, target_depth, target_fraction, denominator):
    """When each sample first reached the target (sequencing time, from the reads' own timestamps)."""
    records = []
    for r in rows:
        t = trackers[r.barcode]
        ref, frac = t.best()
        status = "control" if r.neg_control else ("reached" if t.reached else "not_reached")
        got = t.reached and t.reached_at and run_start
        secs = (t.reached_at - run_start).total_seconds() if got else None
        records.append([
            r.sample, r.barcode, _bool(r.neg_control), status,
            f"{secs / 60:.1f}" if got else "NA", hms(secs) if got else "NA", _utc(t.reached_at) if got else "NA",
            t.reads_at_target if t.reached else "NA", t.bases_at_target if t.reached else "NA",
            (t.reached_ref if t.reached else ref) or "NA", f"{frac:.3f}", t.reads,
            target_depth, target_fraction, denominator,
        ])  # fmt: skip
    _write(path, TIMING_HEADER, records)


def write_timeseries(path, rows, trackers, run_start):
    """Progress over sequencing time: one line per sample per annotation file (best fraction >= target depth)."""
    header = ["sample", "barcode", "neg_control", "time_utc", "elapsed_min", "reads", "mapped_reads",
              "mapped_bases", "best_reference", "best_fraction"]  # fmt: skip
    records = []
    for r in rows:
        for when, reads, mapped, bases, ref, frac in trackers[r.barcode].series:
            records.append([r.sample, r.barcode, _bool(r.neg_control), _utc(when), _minutes(when, run_start),
                            reads, mapped, bases, ref or "NA", f"{frac:.4f}"])  # fmt: skip
    _write(path, header, records)


def write_sample_summary(path, rows, trackers):
    """One line per sample: read counts, dominant reference and read length statistics."""
    header = ["sample", "barcode", "neg_control", "reads", "mapped_reads", "mapped_pct", "mapped_bases",
              "dominant_reference", "dominant_read_pct", "dominant_mean_depth", "best_fraction",
              "total_read_bases", "mean_read_len", "median_read_len", "read_n50", "min_read_len", "max_read_len"]  # fmt: skip
    records = []
    for r in rows:
        t = trackers[r.barcode]
        dom = t.dominant()
        _, frac = t.best()
        ls = t.length_stats()
        records.append([
            r.sample, r.barcode, _bool(r.neg_control), t.reads, t.mapped,
            f"{100 * t.mapped / t.reads:.1f}" if t.reads else "NA", t.bases,
            dom or "NA", f"{100 * t.ref_reads[dom] / t.mapped:.1f}" if dom else "NA",
            f"{t.mean_depth(dom):.1f}" if dom else "NA", f"{frac:.3f}",
            *([ls["total"], f"{ls['mean']:.0f}", f"{ls['median']:.0f}", ls["n50"], ls["min"], ls["max"]] if ls else ["NA"] * 6),
        ])  # fmt: skip
    _write(path, header, records)


def write_reference_summary(path, rows, trackers):
    """One line per sample per reference with coverage: primary read share, mean depth and fraction at target depth."""
    header = ["sample", "barcode", "neg_control", "reference", "reads", "read_pct", "mapped_bases",
              "mean_depth", "fraction_at_target", "dominant"]  # fmt: skip
    records = []
    for r in rows:
        t = trackers[r.barcode]
        dom = t.dominant()
        refs = sorted(t.depth, key=lambda ref: (-t.ref_reads[ref], -t.fraction_of(ref)))
        for ref in refs:
            n = t.ref_reads[ref]
            records.append([
                r.sample, r.barcode, _bool(r.neg_control), ref, n, f"{100 * n / t.mapped:.1f}" if t.mapped else "NA",
                t.ref_bases[ref], f"{t.mean_depth(ref):.1f}", f"{t.fraction_of(ref):.4f}", _bool(ref == dom),
            ])  # fmt: skip
    _write(path, header, records)


def write_all(outdir, rows, trackers, run_start, target_depth, target_fraction, denominator):
    write_timing(outdir / "time_to_target.tsv", rows, trackers, run_start, target_depth, target_fraction, denominator)
    write_timeseries(outdir / "coverage_timeseries.tsv", rows, trackers, run_start)
    write_sample_summary(outdir / "sample_summary.tsv", rows, trackers)
    write_reference_summary(outdir / "reference_summary.tsv", rows, trackers)
