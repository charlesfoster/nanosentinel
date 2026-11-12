"""Simulate a MinKNOW run (fastq_pass/barcodeNN) with ONT-like reads drawn from sequences in a reference panel.

Each --barcode names the panel entry (or entries, joined with '+', picked per read) a barcode's reads come from, so
you can make a 3a sample, a 1a sample, a mixed sample and a near-empty negative control without any sequencing data.

    python dev/simulate_run.py --panel panel.fasta --outdir scratch/sim \
        --barcode barcode01=3a_X76918.1:1500 --barcode barcode02=1a_AF009606.1:1800 \
        --barcode barcode03=3a_X76918.1:25 --barcode barcode04=3a_X76918.1+1a_AF009606.1:1600

Add @START-END (1-based, inclusive, on the reference's own coordinates) to draw a barcode's reads only from that window,
for example an E1E2 amplicon sample on H77:   --barcode barcode05=1a_AF009606.1:600@915-2579

Amplicon (tiling) reads: --protocol protocols/hcv --barcode barcode01=1a_EF407457.1:2500@amplicons
draws reads that start and end at the reference's amplicon sites (from the protocol's amplicons.tsv), in uneven amounts
per amplicon, so a reference gets only the amplicons the scheme has a site for. --divergence 0.03 gives each barcode its
own genome, the reference with 3% substitutions, so reads are not exact copies of a panel sequence.

then: nanosentinel --minknow_data scratch/sim --run_name exp1 ...
Reads are interleaved in time (one read per second overall, ~5% errors, lognormal lengths), written in small chunks
so the monitor sees them as separate files.
"""

import argparse
import datetime
import gzip
import io
import math
import os
import random

COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def mutate(seq, rng, err):
    """Apply substitutions, insertions and deletions at a total per-base rate `err`."""
    out = []
    for base in seq:
        if rng.random() < err:
            kind = rng.random()
            if kind < 0.4:
                out.append(rng.choice("ACGT"))
            elif kind < 0.7:
                out.append(base + rng.choice("ACGT"))
            # else: deletion
        else:
            out.append(base)
    return "".join(out)


def simulate_read(seq, rng, median_len, err):
    length = min(len(seq), max(100, int(rng.lognormvariate(math.log(median_len), 0.6))))
    start = rng.randint(0, len(seq) - length)
    read = seq[start : start + length]
    if rng.random() < 0.5:
        read = read.translate(COMP)[::-1]
    return mutate(read, rng, err)


def load_amplicons(protocol_dir):
    """{reference: [(name, start0, end)]} from a protocol's amplicons.tsv (1-based, inclusive in the file)."""
    spans = {}
    with open(os.path.join(protocol_dir, "amplicons.tsv")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            spans.setdefault(row["reference"], []).append((row["amplicon"], int(row["start"]) - 1, int(row["end"])))
    return spans


def diverge(seq, rng, rate):
    """The sequence with substitutions at ``rate`` (no indels, so positions stay comparable to the reference)."""
    if rate <= 0:
        return seq
    out = list(seq)
    for i, base in enumerate(out):
        if base in "ACGT" and rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(out)


def amplicon_read(seq, start, end, rng, err):
    """A read from one amplicon: usually the whole of it (ends trimmed a little), sometimes a fragment."""
    if rng.random() < 0.85:
        a, b = start + rng.randint(0, 20), end - rng.randint(0, 20)
    else:
        a, b = sorted(rng.sample(range(start, end), 2))
        b = max(b, min(end, a + 100))
    read = seq[a:b]
    if rng.random() < 0.5:
        read = read.translate(COMP)[::-1]
    return mutate(read, rng, err)


def parse_barcode(text):
    """barcodeNN=REF[+REF2]:COUNT[@START-END | @amplicons] -> (name, [refs], count, (start, end) | 'amplicons' | None)"""
    try:
        name, rest = text.split("=", 1)
        window = None
        if "@" in rest:
            rest, span = rest.rsplit("@", 1)
            if span == "amplicons":
                window = "amplicons"
            else:
                start, end = span.split("-")
                window = (int(start), int(end))
                if not 1 <= window[0] < window[1]:
                    raise ValueError
        refs, count = rest.rsplit(":", 1)
        return name, refs.split("+"), int(count), window
    except ValueError:
        raise SystemExit(f"bad --barcode {text!r}; expected barcode01=REF[+REF2]:COUNT[@START-END|@amplicons]")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, help="multi-FASTA the reads are drawn from")
    ap.add_argument("--outdir", required=True, help="becomes the --minknow_data directory")
    ap.add_argument("--barcode", action="append", required=True, help="barcodeNN=REF[+REF2]:COUNT (repeatable)")
    ap.add_argument("--run_name", default="exp1")
    ap.add_argument("--median_len", type=int, default=900)
    ap.add_argument("--error", type=float, default=0.05)
    ap.add_argument("--chunk", type=int, default=200, help="reads per FASTQ file")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flat_quality", action="store_true", help="one quality value per read instead of one per base: far smaller files (nanosentinel does not read base qualities)")
    ap.add_argument("--protocol", help="protocol folder whose amplicons.tsv defines the amplicons for @amplicons barcodes")
    ap.add_argument("--divergence", type=float, default=0.0, help="substitution rate between each barcode's genome and its reference (default: 0)")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    panel = read_fasta(args.panel)
    amplicons = {}
    if any(b.endswith("@amplicons") for b in args.barcode):
        if not args.protocol:
            raise SystemExit("@amplicons needs --protocol (a folder with amplicons.tsv)")
        amplicons = load_amplicons(args.protocol)
    specs = [parse_barcode(b) for b in args.barcode]
    for name, refs, _, window in specs:
        for ref in refs:
            if ref not in panel:
                raise SystemExit(f"{ref} (for {name}) is not in {args.panel}")
            if window == "amplicons" and ref not in amplicons:
                raise SystemExit(f"{ref} (for {name}) has no amplicons in the protocol's amplicons.tsv")
            if window and window != "amplicons" and window[1] > len(panel[ref]):
                raise SystemExit(f"window {window[0]}-{window[1]} for {name} is past the end of {ref} ({len(panel[ref])})")

    # One shuffled timeline so barcodes arrive interleaved, as on a real flow cell.
    timeline = [name for name, _, count, _ in specs for _ in range(count)]
    rng.shuffle(timeline)
    refs_for = {name: refs for name, refs, _, _ in specs}
    window_for = {name: window for name, _, _, window in specs}

    run = os.path.join(args.outdir, args.run_name, "sample1", "20261109_1000_P2S-00001-A_FLOabc_deadbeef")
    t0 = datetime.datetime(2026, 11, 9, 10, 0, 0)
    seen = {name: 0 for name, _, _, _ in specs}
    genomes, weights = {}, {}  # per barcode and reference: its own genome, and how much of each amplicon it gets
    for name, refs, _, window in specs:
        for ref in refs:
            genomes[(name, ref)] = diverge(panel[ref], random.Random(f"{args.seed}/{name}/{ref}"), args.divergence)
            if window == "amplicons":
                wrng = random.Random(f"{args.seed}/weights/{name}/{ref}")
                weights[(name, ref)] = [wrng.lognormvariate(0, 0.5) for _ in amplicons[ref]]
    handles = {}
    for i, name in enumerate(timeline):
        n = seen[name]
        seen[name] += 1
        chunk = n // args.chunk
        if (name, chunk) not in handles:
            for key in [k for k in handles if k[0] == name]:
                handles.pop(key).close()
            d = os.path.join(run, "fastq_pass", name)
            os.makedirs(d, exist_ok=True)
            # mtime=0 keeps the gzip header the same on every run, so regenerated files are byte-identical
            raw = gzip.GzipFile(os.path.join(d, f"FLOabc_pass_{name}_deadbeef_{chunk}.fastq.gz"), "wb", mtime=0)
            handles[(name, chunk)] = io.TextIOWrapper(raw, encoding="ascii", newline="\n")
        ref = rng.choice(refs_for[name])
        window = window_for[name]
        genome = genomes[(name, ref)]
        if window == "amplicons":
            _, a_start, a_end = rng.choices(amplicons[ref], weights[(name, ref)])[0]
            read = amplicon_read(genome, a_start, a_end, rng, args.error)
        else:
            source = genome[window[0] - 1 : window[1]] if window else genome
            read = simulate_read(source, rng, args.median_len, args.error)
        stamp = (t0 + datetime.timedelta(seconds=i)).strftime("%Y-%m-%dT%H:%M:%SZ")
        if args.flat_quality:
            qual = chr(33 + max(3, min(40, int(rng.gauss(14, 4))))) * len(read)
        else:
            qual = "".join(chr(33 + max(3, min(40, int(rng.gauss(14, 4))))) for _ in read)
        handles[(name, chunk)].write(
            f"@{name}_r{n} runid=deadbeef read={i} ch={i % 512 + 1} start_time={stamp} barcode={name}\n{read}\n+\n{qual}\n"
        )
    for h in handles.values():
        h.close()
    print(run)
    for name, refs, count, window in specs:
        print(f"  {name}: {count} reads from {'+'.join(refs)}" + (" (amplicon reads)" if window == "amplicons" else f" (positions {window[0]}-{window[1]} only)" if window else ""))


if __name__ == "__main__":
    main()
