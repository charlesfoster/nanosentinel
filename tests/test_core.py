"""Three small tests, for the places where a bug would be silent."""

import gzip
from datetime import datetime, timedelta, timezone

import pytest

from nanosentinel.concat import concatenate
from nanosentinel.coverage import SampleTracker
from nanosentinel.errors import SheetError
from nanosentinel.samplesheet import read_samplesheet, select_run


def test_time_to_target_is_when_the_fraction_is_first_reached():
    t0 = datetime(2026, 1, 1, 10, 0, tzinfo=timezone.utc)
    tracker = SampleTracker({"R": 100}, target_depth=2, target_fraction=0.9)
    assert tracker.add("R", 0, 100, t0) is False  # depth 1 everywhere
    assert tracker.add("R", 0, 89, t0 + timedelta(minutes=1)) is False  # 89/100 at depth 2
    assert tracker.add("R", 0, 50, t0 + timedelta(minutes=2)) is False  # already >= 2x: no change
    assert tracker.add("R", 89, 95, t0 + timedelta(minutes=3)) is True  # 95/100
    assert tracker.reached_at == t0 + timedelta(minutes=3)
    assert tracker.reads_at_target == 4
    assert tracker.add("R", 0, 100, t0 + timedelta(minutes=9)) is False  # never fires twice
    assert tracker.reached_at == t0 + timedelta(minutes=3)


def test_samplesheet_barcodes_controls_and_repeated_samples(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text(
        "sample,long_reads,barcode,run_name,neg_control\n"
        "A,,BC01,run1,\n"
        "A,,barcode2,run1,false\n"  # repeated sample on a second barcode is fine
        "NTC,,3,run1,TRUE\n"
        "B,,BC01,run2,\n"  # a different run may reuse barcodes
    )
    header, rows = read_samplesheet(sheet)
    assert header[:3] == ["sample", "long_reads", "barcode"]  # original columns are kept
    run, picked = select_run(rows, "run1")
    assert [r.barcode for r in picked] == ["barcode01", "barcode02", "barcode03"]
    assert [r.neg_control for r in picked] == [False, False, True]
    with pytest.raises(SheetError):  # sheet spans two runs and none was chosen
        select_run(rows)
    sheet.write_text("sample,barcode\nA,BC01\nB,barcode1\n")  # same barcode, two samples
    with pytest.raises(SheetError):
        select_run(read_samplesheet(sheet)[1])


def test_concatenated_gzip_reads_back_in_order_and_sheet_gets_long_reads(tmp_path):
    fastq_pass = tmp_path / "fastq_pass" / "barcode01"
    fastq_pass.mkdir(parents=True)
    for i in (2, 10, 1):  # natural order must be 1, 2, 10 - not 1, 10, 2
        with gzip.open(fastq_pass / f"x_{i}.fastq.gz", "wt") as fh:
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample,barcode\nA,BC01\n")
    header, rows = read_samplesheet(sheet)
    out, n = concatenate(rows, header, tmp_path / "fastq_pass", tmp_path / "out", "run1")
    assert n == 1
    reads = gzip.open(tmp_path / "out" / "reads" / "A_run1_barcode01.fastq.gz", "rt").read()
    assert [l for l in reads.splitlines() if l.startswith("@")] == ["@r1", "@r2", "@r10"]
    assert out.read_text().splitlines()[0] == "sample,long_reads,barcode"  # long_reads inserted after sample


def test_read_length_stats_and_dominant_reference():
    tracker = SampleTracker({"A": 1000, "B": 1000}, target_depth=1, target_fraction=0.99)
    for length in (100, 200, 300, 400, 1000):
        tracker.add("A", 0, length, None, length)
    tracker.add("B", 0, 50, None, 50)
    tracker.add("*", 0, 0, None, 70)  # unmapped reads still count towards length stats
    stats = tracker.length_stats()
    assert (stats["n"], stats["total"], stats["min"], stats["max"]) == (7, 2120, 50, 1000)
    # sorted: 50 70 100 200 300 400 1000 -> median 200; half of 2120 bases is 1060, reached by 1000 + 400 -> N50 400
    assert stats["n50"] == 400 and stats["median"] == 200
    assert tracker.dominant() == "A" and tracker.ref_reads["A"] == 5 and tracker.mapped == 6
