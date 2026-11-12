"""dev/simulate_run.py: amplicon reads stay inside their amplicons, and divergence changes the genome."""

import gzip
import importlib.util
import random
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "dev" / "simulate_run.py"
spec = importlib.util.spec_from_file_location("simulate_run", SCRIPT)
sim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sim)


def test_amplicon_reads_start_and_end_inside_their_amplicon():
    rng = random.Random(3)
    genome = "".join(random.Random(1).choice("ACGT") for _ in range(2000))
    start, end = 500, 1300
    for _ in range(300):
        read = sim.amplicon_read(genome, start, end, rng, err=0.0)
        forward = read in genome[start:end]
        reverse = read.translate(sim.COMP)[::-1] in genome[start:end]
        assert forward or reverse  # error-free reads are an exact piece of the amplicon, on either strand
        assert 100 <= len(read) <= end - start


def test_divergence_substitutes_only_and_is_repeatable():
    genome = "".join(random.Random(2).choice("ACGT") for _ in range(5000))
    a = sim.diverge(genome, random.Random("x"), 0.03)
    assert len(a) == len(genome) and a == sim.diverge(genome, random.Random("x"), 0.03)
    differences = sum(x != y for x, y in zip(genome, a))
    assert 0.02 * len(genome) < differences < 0.04 * len(genome)
    assert sim.diverge(genome, random.Random("x"), 0) == genome


def test_amplicon_barcodes_only_use_amplicons_the_reference_has(tmp_path):
    panel = tmp_path / "panel.fasta"
    seqs = {n: "".join(random.Random(i).choice("ACGT") for _ in range(3000)) for i, n in enumerate(["A", "B"])}
    panel.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    protocol = tmp_path / "protocol"
    protocol.mkdir()
    (protocol / "amplicons.tsv").write_text(
        "reference\tamplicon\tstart\tend\n" "A\tSeg1\t101\t900\n" "A\tSeg2\t801\t1600\n" "B\tSeg1\t101\t900\n"  # B lacks Seg2
    )
    out = tmp_path / "sim"
    import sys

    argv = ["simulate_run", "--panel", str(panel), "--protocol", str(protocol), "--outdir", str(out), "--error", "0",
            "--barcode", "barcode01=B:200@amplicons"]  # fmt: skip
    old, sys.argv = sys.argv, argv
    try:
        sim.main()
    finally:
        sys.argv = old
    reads = []
    for path in out.rglob("*.fastq.gz"):
        with gzip.open(path, "rt") as fh:
            reads += [line.strip() for i, line in enumerate(fh) if i % 4 == 1]
    assert len(reads) == 200
    window = seqs["B"][100:900]
    assert all(r in window or r.translate(sim.COMP)[::-1] in window for r in reads)  # nothing from outside B's one amplicon
