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


def test_paf_parsing_primary_secondary_chimeric_and_identity_filter():
    from nanosentinel.annotate import parse_paf

    def line(read, ref, start, end, matches, block, tp, score, mapq=0):
        return f"{read}\t500\t0\t500\t+\t{ref}\t9000\t{start}\t{end}\t{matches}\t{block}\t{mapq}\ttp:A:{tp}\tAS:i:{score}"

    paf = [
        line("r1", "A", 100, 600, 480, 500, "P", 900, 40),  # primary, 96% identity
        line("r1", "B", 100, 600, 470, 500, "S", 860),  # near-equal sibling: kept, but not a primary
        line("r2", "A", 100, 600, 350, 500, "P", 400),  # 70% identity: below 0.8, so r2 is unmapped
        line("r3", "A", 100, 300, 190, 200, "P", 300),  # chimeric: two primaries; the higher score wins
        line("r3", "A", 2000, 2400, 380, 400, "P", 700),
    ]
    out = parse_paf(paf, min_identity=0.8)
    assert [(r[0], r[-1]) for r in out["r1"]] == [("A", 1), ("B", 0)]
    assert "r2" not in out
    assert [(r[2], r[-1]) for r in out["r3"]] == [(2000, 1), (100, 0)]


# ---- protocols: coordinate maps, primers, lifting, projection -------------------------------------------------------

import random

import numpy as np

from nanosentinel import protocol as proto
from nanosentinel.references import ReferenceInfo


def _random_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_coordinate_map_handles_insertions_and_deletions_exactly():
    # coordinate ACGTACGT; R inserts GG after position 3 (a gap column in the coordinate row) and lacks position 6
    aligned = {"C": "ACGT--ACGT", "R": "ACGTGGAC-T"}
    maps = proto.build_maps(aligned, "C", {"C": "ACGTACGT", "R": "ACGTGGACT"})
    assert maps["C"].tolist() == [0, 1, 2, 3, 4, 5, 6, 7]
    assert maps["R"].tolist() == [0, 1, 2, 3, -1, -1, 4, 5, 7]  # the inserted GG has no coordinate; deletion skipped 6
    with pytest.raises(proto.ProtocolError):  # an alignment of some other sequence is rejected, not silently used
        proto.build_maps(aligned, "C", {"C": "ACGTACGT", "R": "ACGTGGACA"})


def test_projection_and_lifting_through_a_map():
    info = ReferenceInfo({"C": 8, "R": 9}, "C", maps={"C": np.arange(8, dtype=np.int32),
                         "R": np.array([0, 1, 2, 3, -1, -1, 4, 5, 7], dtype=np.int32)})  # fmt: skip
    depth = np.arange(1, 10, dtype=np.uint32) * 10  # depth on R: 10 .. 90
    values, valid = info.project("R", depth)
    assert valid.tolist() == [True, True, True, True, True, True, False, True]  # coordinate position 6 is deleted in R
    assert values.tolist() == [10, 20, 30, 40, 70, 80, 0, 90]
    assert info.to_coordinate("R", 4, 6) is None  # an interval that is all insertion has no place on the axis
    assert info.to_coordinate("R", 3, 9) == (3, 8)
    assert info.project("X", depth) is None  # a reference without a map cannot be projected
    assert info.group("1a_AF009606.1") == "1a" and info.group("plain") == "plain"


def test_primer_matching_is_iupac_aware_and_protects_the_3prime_end():
    codes = proto._CODES[np.frombuffer(b"TTTACGTTTT", dtype=np.uint8)]
    total, tail = proto.mismatches(codes, "ACGT")
    assert total[3] == 0  # exact site
    total, tail = proto.mismatches(codes, "ACRT")
    assert total[3] == 0  # R = A or G
    total, tail = proto.mismatches(codes, "TCGT")
    assert total[3] == 1 and tail[3] == 1  # the mismatch is at the primer's 3' end (last 5 bases cover all of it)
    n_codes = proto._CODES[np.frombuffer(b"ACNT", dtype=np.uint8)]
    assert proto.mismatches(n_codes, "ACGT")[0][0] == 1  # N never counts as a match
    assert proto.reverse_complement("AACGTR") == "YACGTT"


def _panel(seed=1, n_refs=4):
    coord = _random_seq(600, seed)
    seqs = {"C": coord}
    rng = random.Random(seed + 1)
    for i in range(n_refs - 1):
        s = list(coord)
        for pos in rng.sample(range(100, 200), 6) + rng.sample(range(330, 500), 6):  # away from the primer sites
            s[pos] = rng.choice([b for b in "ACGT" if b != s[pos]])
        seqs[f"R{i}"] = "".join(s)
    return seqs


def test_primers_are_placed_near_their_home_and_far_chance_hits_are_ignored(tmp_path):
    seqs = _panel()
    c = seqs["C"]
    left, right = c[50:70], proto.reverse_complement(c[250:270])
    maps = {k: np.arange(len(v), dtype=np.int32) for k, v in seqs.items()}
    primers = [{"name": "A1_LEFT", "seq": left, "amplicon": "A1", "side": "LEFT"},
               {"name": "A1_RIGHT", "seq": right, "amplicon": "A1", "side": "RIGHT"}]  # fmt: skip
    # plant a second, near-perfect copy of the left primer far away on one reference: it must not be chosen
    far = list(seqs["R0"]); far[500:520] = left; seqs["R0"] = "".join(far)
    sites, expected, unplaced = proto.place_primers(primers, seqs, maps)
    assert not unplaced and expected["A1_LEFT"][0] in range(55, 66)
    assert all(sites[ref]["A1_LEFT"][:2] == (50, 70) for ref in seqs)
    assert all(sites[ref]["A1_RIGHT"][:2] == (250, 270) for ref in seqs)
    amplicons = proto.build_amplicons(primers, sites)
    assert [a[:3] for a in amplicons["C"]] == [("A1", 50, 270)]


def test_features_are_lifted_with_a_status():
    maps = {"C": np.arange(100, dtype=np.int32),
            "SHORT": np.arange(40, dtype=np.int32)}  # a partial sequence covering only the first 40 positions
    lifted = proto.lift_features({"g": (11, 30), "h": (31, 90)}, maps, {"C": 100, "SHORT": 40})
    assert lifted["C"]["g"] == (11, 30, 20, 20, "ok")
    assert lifted["SHORT"]["g"][4] == "ok" and lifted["SHORT"]["h"][4] == "absent"


def test_protocol_round_trip_with_a_ready_made_alignment(tmp_path):
    seqs = _panel(seed=5, n_refs=5)
    panel = tmp_path / "panel.fasta"
    panel.write_text("".join(f">{'1a' if k == 'C' else '2b'}_{k}\n{v}\n" for k, v in seqs.items()))
    ids = {k: f"{'1a' if k == 'C' else '2b'}_{k}" for k in seqs}
    alignment = tmp_path / "aln.fasta"
    alignment.write_text("".join(f">{ids[k]}\n{v}\n" for k, v in seqs.items()))  # substitutions only: already aligned
    c = seqs["C"]
    primers = tmp_path / "primers.fa"
    primers.write_text(f">Seg1_all_LEFT\n{c[50:70]}\n>Seg1_all_RIGHT\n{proto.reverse_complement(c[250:270])}\n")
    gff = tmp_path / "g.gff"
    gff.write_text(f"{ids['C']}\tx\tgene\t11\t100\t.\t+\t.\tID=G1;Name=G1\n")
    out = tmp_path / "protocol"
    qc = proto.build_protocol(out, panel, coordinate_ref=ids["C"], gff=gff, regions={"R": (11, 60)}, primer_fasta=primers,
                              amplicon_pattern=r"^(Seg[0-9]+)_", panel_alignment=alignment, name="toy", say=lambda _: None)  # fmt: skip
    assert qc["references_with_all_amplicons"] == 5 and qc["features"]["gene:G1"]["status"] == {"ok": 5}
    info = proto.load_protocol(out, tmp_path / "run")
    assert info.protocol == "toy" and info.coordinate_ref == ids["C"] and info.amplicon_names == ["Seg1"]
    assert info.amplicons[ids["R0"]][0][:3] == ("Seg1", 50, 270)
    assert info.genes["G1"]["end"] == 100 and info.regions["R"]["start"] == 11 and info.index.exists()
    assert info.mask(ids["C"]).sum() == 220  # amplicon mode counts only bases inside the amplicons
    (out / "panel.fasta").write_text((out / "panel.fasta").read_text() + ">extra\nACGT\n")
    with pytest.raises(proto.ProtocolError, match="changed since the protocol was built"):
        proto.load_protocol(out, tmp_path / "run2")


def test_sibling_group_lock_does_not_flip_between_near_identical_references():
    group = lambda ref: ref.split("_")[0]  # noqa: E731
    t = SampleTracker({"1a_x": 1000, "1a_y": 1000, "3a_z": 1000}, target_depth=1, target_fraction=0.99, group_of=group)
    for _ in range(150):
        t.add("1a_x", 0, 500, None, 500)
    for _ in range(100):
        t.add("1a_y", 0, 500, None, 500)  # siblings split the reads
    t.snapshot()
    assert t.locked_group == "1a" and t.shown() in ("1a_x", "1a_y")
    for _ in range(60):
        t.add("3a_z", 0, 500, None, 500)  # a few reads elsewhere do not move the lock
    t.snapshot()
    assert t.locked_group == "1a"
    for _ in range(600):
        t.add("3a_z", 0, 500, None, 500)  # but a clear, large lead does
    t.snapshot()
    assert t.locked_group == "3a" and t.shown() == "3a_z"


def test_region_target_is_held_to_the_whole_region():
    from types import SimpleNamespace

    from nanosentinel.errors import MonitorError
    from nanosentinel.monitor import build_trackers

    maps = {
        "C": np.arange(100, dtype=np.int32),
        "SHORT": np.arange(40, dtype=np.int32),  # a partial sequence: covers coordinates 0-39 only
        "DEL": np.concatenate([np.arange(0, 40), np.arange(45, 100)]).astype(np.int32),  # lacks coordinates 40-44
    }
    info = ReferenceInfo({"C": 100, "SHORT": 40, "DEL": 95}, "C", regions={"G": {"start": 21, "end": 60}}, maps=maps)
    rows = [SimpleNamespace(barcode="b")]

    def fraction(ref, length):
        tracker = build_trackers(rows, info, 1, 0.9, "region:G")["b"]
        tracker.add(ref, 0, length, None, length)  # every base of the reference covered once
        return tracker.fraction_of(ref), tracker.reached

    assert fraction("C", 100) == (1.0, True)
    assert fraction("DEL", 95) == (35 / 40, False)  # 5 of the region's 40 positions do not exist on this reference
    assert fraction("SHORT", 40) == (20 / 40, False)  # half of the region is missing: it can never reach 90%
    with pytest.raises(MonitorError, match="no region of that name"):
        build_trackers(rows, info, 1, 0.9, "region:NOPE")


# ---- read start times: both header styles, and old annotations are redone ----------------------------------------------

REAL_HEADER = (  # a header from current MinKNOW/dorado (SAM-style tags), as seen on a real run
    "@292ba5f2-d707-4d95-8f13-08a8889dc60f qs:f:20.157257 mx:i:1 ch:i:148 rn:i:1 st:Z:2026-11-10T14:08:56.683355+11:00 "
    "ns:i:11123 ts:i:2266 sm:f:93.6924 sd:f:23.506744 sv:Z:pa du:f:2.2246 dx:i:0 "
    "RG:Z:17430652-b625-4d30-8eda-7d1b9a588ad2_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_barcode01 "
    "DT:Z:2026-11-10T14:06:54.012955+11:00 PU:Z:FBH55487 LB:Z:unknown SM:Z:barcode01 al:Z:barcode01"
)


def test_read_start_time_is_found_in_both_header_styles(tmp_path):
    from nanosentinel.annotate import read_fastq_info
    from nanosentinel.coverage import parse_time

    fq = tmp_path / "x.fastq"
    fq.write_text(
        f"{REAL_HEADER}\nACGT\n+\nIIII\n"  # current style: st:Z: (DT:Z: is the run start, not the read start)
        "@r2 runid=x read=2 ch=3 start_time=2026-11-10T14:09:01Z barcode=barcode01\nACGT\n+\nIIII\n"  # older style
        "@r3 runid=x ch:i:5\nACGT\n+\nIIII\n"  # no start time at all
    )
    reads = read_fastq_info(fq)
    assert [r[0] for r in reads] == ["292ba5f2-d707-4d95-8f13-08a8889dc60f", "r2", "r3"]
    assert reads[0][2] == "2026-11-10T14:08:56.683355+11:00" and reads[1][2] == "2026-11-10T14:09:01Z" and reads[2][2] == ""
    assert parse_time(reads[0][2]).isoformat() == "2026-11-10T14:08:56.683355+11:00"


def test_annotations_made_before_the_format_change_are_redone_not_reused(tmp_path):
    from nanosentinel.annotate import annotation_is_current, write_annotation
    from nanosentinel.coverage import AnnotationReader, SampleTracker

    new, old = tmp_path / "new.tsv", tmp_path / "old.tsv"
    write_annotation(new, 100, "abc123", [("r1", 500, "2026-11-10T14:08:56+11:00")], {})
    assert annotation_is_current(new, 100, "abc123")
    # the same file as an earlier version wrote it: no format field (and, for real runs, blank start times)
    old.write_text("#source_size=100 panel=abc123\nread_name\tread_len\tstart_time\tref\tref_len\tstart\tend\tmapq\tidentity\tprimary\n")
    assert not annotation_is_current(old, 100, "abc123")  # so the annotator maps it again
    reader = AnnotationReader(tmp_path, {"b": SampleTracker({"R": 10})}, "abc123")
    assert reader._read_file(old, reader.trackers["b"]) is None  # and the reader does not count it in the meantime
