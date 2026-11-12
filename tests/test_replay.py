"""dev/replay_run.py: order and pacing come from the reads' timestamps, and a half-copied chunk is never visible."""

import gzip
import importlib.util
import os
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "dev" / "replay_run.py"
spec = importlib.util.spec_from_file_location("replay_run", SCRIPT)
replay_run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay_run)


def write_chunk(folder, name, seconds, header_time=True):
    """A gzipped FASTQ whose reads start at the given seconds after 10:00:00."""
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / name, "wt") as fh:
        for i, s in enumerate(seconds):
            stamp = f" start_time=2026-11-09T10:{s // 60:02d}:{s % 60:02d}Z" if header_time else ""
            fh.write(f"@{name}_{i} runid=x{stamp}\nACGT\n+\nIIII\n")


@pytest.fixture
def run(tmp_path):
    """Two barcodes whose file names and modification times both lie about the order."""
    src = tmp_path / "fastq_pass"
    write_chunk(src / "barcode01", "a_0.fastq.gz", [0, 10, 20])  # complete at 20 s
    write_chunk(src / "barcode01", "a_1.fastq.gz", [30, 100])  # complete at 100 s
    write_chunk(src / "barcode02", "b_0.fastq.gz", [5, 45])  # complete at 45 s
    write_chunk(src / "barcode02", "b_1.fastq.gz", [60, 70])  # complete at 70 s
    (src / "barcode01" / "notes.txt").write_text("not a FASTQ")
    for number, path in enumerate(sorted(src.rglob("*.fastq.gz"), reverse=True)):  # scramble the mtimes
        os.utime(path, (1_000_000 + number, 1_000_000 + number))
    return src


def test_chunks_are_ordered_and_paced_by_their_read_timestamps_not_file_times(run):
    chunks = replay_run.build_schedule(replay_run.find_chunks(run))
    assert [str(c.rel) for c in chunks] == [
        "barcode01/a_0.fastq.gz", "barcode02/b_0.fastq.gz", "barcode02/b_1.fastq.gz", "barcode01/a_1.fastq.gz",
    ]  # fmt: skip
    assert [c.offset for c in chunks] == [20, 45, 70, 100]  # seconds from the first read of the run
    assert [c.reads for c in chunks] == [3, 2, 2, 2]  # the .txt file was ignored


def test_pacing_follows_speed(run, tmp_path):
    chunks = replay_run.build_schedule(replay_run.find_chunks(run))
    now = [0.0]
    released = []

    def sleep(seconds):
        now[0] += seconds

    def watch(*args):  # record the fake time at which each chunk lands
        released.append(now[0])

    replay_run.replay(chunks, tmp_path / "dest", speed=10, clock=lambda: now[0], sleep=sleep, say=watch)
    assert released == pytest.approx([2.0, 4.5, 7.0, 10.0])  # offsets / 10
    released.clear()
    now[0] = 0.0
    slept = []
    replay_run.replay(chunks, tmp_path / "dest0", speed=0, clock=lambda: now[0], sleep=slept.append, say=watch)
    assert slept == []  # --speed 0 copies everything at once


def test_a_half_copied_chunk_is_never_visible_and_the_final_files_are_identical(run, tmp_path, monkeypatch):
    chunks = replay_run.build_schedule(replay_run.find_chunks(run))
    dest = tmp_path / "dest"
    real_replace = os.replace
    seen = []

    def checked_replace(src, dst):
        # at the moment of the rename the final name must not exist yet, and the temporary name must not look like a FASTQ
        assert not Path(dst).exists()
        assert not Path(src).name.endswith(replay_run.FASTQ_SUFFIXES)
        seen.append(Path(dst).name)
        real_replace(src, dst)

    monkeypatch.setattr(replay_run.os, "replace", checked_replace)
    replay_run.replay(chunks, dest, speed=0, say=lambda *_: None)
    assert len(seen) == 4
    assert not list(dest.rglob("*.partial"))  # nothing left behind
    for chunk in chunks:
        assert (dest / chunk.rel).read_bytes() == chunk.path.read_bytes()  # byte-for-byte, gzip untouched


def test_a_chunk_without_timestamps_is_placed_after_its_predecessor_with_a_warning(tmp_path):
    src = tmp_path / "fastq_pass"
    write_chunk(src / "barcode01", "a_0.fastq.gz", [0, 50])
    write_chunk(src / "barcode01", "a_1.fastq.gz", [1, 2], header_time=False)
    warnings = []
    chunks = replay_run.build_schedule(replay_run.find_chunks(src), warn=warnings.append)
    assert [c.offset for c in chunks] == [50, 50] and len(warnings) == 1
    write_chunk(src / "barcode02", "b_0.fastq.gz", [1, 2], header_time=False)
    for path in (src / "barcode01").glob("*.gz"):
        path.unlink()
    with pytest.raises(SystemExit, match="start_time"):  # no timestamps anywhere: cannot order or pace
        replay_run.build_schedule(replay_run.find_chunks(src))


def test_refuses_a_destination_that_already_holds_reads(run, tmp_path):
    dest = tmp_path / "dest"
    write_chunk(dest / "barcode01", "old.fastq.gz", [0])
    with pytest.raises(SystemExit):
        replay_run.main([str(run), str(dest), "--speed", "0"])
