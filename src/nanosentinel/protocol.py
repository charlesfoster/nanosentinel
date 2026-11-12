"""Protocols: a reference panel prepared once (aligned, annotated, indexed) and reused for every run.

``prep-protocol`` aligns every reference in the panel to the coordinate reference (MAFFT), which gives each reference a
*coordinate map* (its bases -> positions on the coordinate reference). Genes and named regions are given on the
coordinate reference; a primer scheme is either a BED on it or a FASTA of primers whose sites are searched for on every
reference. Everything is written to one folder that ``monitor --protocol`` loads in about a second:

    protocol.json         what it is: inputs, tool versions, date, genes, regions, QC summary, file hashes
    panel.fasta           the reference panel (this is what reads are mapped to)
    coordinate_maps.npz   per-reference coordinate maps
    amplicons.tsv         amplicon span on every reference (human-readable; 1-based, inclusive)
    primer_sites.tsv      where each primer was found on each reference, with mismatches
    features.tsv          genes and regions lifted onto every reference, with a status per lift
    qc.txt                the QC report

``monitor`` refuses to start if panel.fasta or the maps no longer match the hashes in protocol.json.
"""

import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__, tools
from .errors import ProtocolError
from .references import (
    DEFAULT_GROUP_PATTERN,
    ReferenceInfo,
    _fasta_text,
    build_index,
    parse_bed_amplicons,
    parse_gff_genes,
    read_fasta,
    read_sequences,
)

FORMAT = 1
BUNDLED = Path(__file__).resolve().parents[2] / "protocols"
POOR_ALIGNMENT = 0.6  # a reference with fewer of its bases placed than this is left out of the shared axis

# IUPAC codes as bit masks (A=1 C=2 G=4 T=8). An unknown base or N matches nothing, so runs of N never look like a site.
_BITS = {"A": 1, "C": 2, "G": 4, "T": 8, "U": 8, "R": 5, "Y": 10, "S": 6, "W": 9, "K": 12, "M": 3,
         "B": 14, "D": 13, "H": 11, "V": 7}  # fmt: skip
_CODES = np.zeros(256, dtype=np.uint8)
for _c, _v in _BITS.items():
    _CODES[ord(_c)] = _CODES[ord(_c.lower())] = _v
_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")


def _natural(text):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


def sha1_file(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---- alignment and coordinate maps ----------------------------------------------------------------------------------


def run_mafft(coordinate_id, sequences, threads=4):
    """Align every sequence to the coordinate reference: ``mafft --add`` onto the coordinate sequence alone.

    Insertions relative to the coordinate reference are kept (as gap columns in it), so no base is lost and every
    reference base can be mapped exactly. Returns {id: aligned sequence}.
    """
    others = {k: v for k, v in sequences.items() if k != coordinate_id}
    if not others:
        return {coordinate_id: sequences[coordinate_id]}
    with tempfile.TemporaryDirectory(prefix="nanosentinel_mafft_") as tmp:
        tmp = Path(tmp)
        (tmp / "coordinate.fa").write_text(f">{coordinate_id}\n{sequences[coordinate_id]}\n")
        (tmp / "others.fa").write_text("".join(f">{k}\n{v}\n" for k, v in others.items()))
        cmd = [tools.find("mafft"), "--thread", str(threads), "--quiet", "--add", str(tmp / "others.fa"),
               str(tmp / "coordinate.fa")]  # fmt: skip
        done = subprocess.run(cmd, capture_output=True, text=True)
        if done.returncode != 0:
            raise ProtocolError(f"MAFFT failed:\n{done.stderr.strip()[-500:]}")
    return parse_alignment(done.stdout.splitlines())


def parse_alignment(lines):
    """{id: aligned sequence} (upper case, U read as T, '.' read as a gap)."""
    rows, name = {}, None
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name is not None:
            rows[name].append(line.strip())
    return {k: "".join(v).upper().replace("U", "T").replace(".", "-") for k, v in rows.items()}


def build_maps(aligned, coordinate_id, sequences):
    """Coordinate maps from an alignment: {id: int32 array, base -> coordinate position or -1}.

    Each aligned row must be its sequence with gaps added; anything else means the alignment is not of this panel.
    """
    if coordinate_id not in aligned:
        raise ProtocolError(f"The alignment does not contain the coordinate reference {coordinate_id}")
    widths = {len(v) for v in aligned.values()}
    if len(widths) != 1:
        raise ProtocolError("The sequences in the alignment are not all the same length")
    coord_row = np.frombuffer(aligned[coordinate_id].encode(), dtype=np.uint8)
    coord_residue = coord_row != ord("-")
    coord_pos = np.cumsum(coord_residue) - 1  # coordinate position at each column (meaningful where it has a base)
    maps = {}
    for ref, seq in sequences.items():
        if ref not in aligned:
            raise ProtocolError(f"{ref} is in the panel but not in the alignment")
        row = np.frombuffer(aligned[ref].encode(), dtype=np.uint8)
        residue = row != ord("-")
        if aligned[ref].replace("-", "") != seq:
            raise ProtocolError(f"{ref} in the alignment is not the same sequence as in the panel")
        maps[ref] = np.where(coord_residue[residue], coord_pos[residue], -1).astype(np.int32)
    return maps


# ---- primers --------------------------------------------------------------------------------------------------------


def read_primers(path, amplicon_pattern=None):
    """Primers from a FASTA: [{name, seq, amplicon, side}]. ``amplicon_pattern`` has one group naming the amplicon."""
    pattern = re.compile(amplicon_pattern or r"^(.*?)_(?:LEFT|RIGHT)", re.IGNORECASE)
    primers = []
    for name, seq in read_sequences(path).items():
        match = pattern.search(name)
        side = re.search(r"LEFT|RIGHT", name, re.IGNORECASE)
        if not match or not match.groups() or not side:
            raise ProtocolError(
                f"Cannot tell which amplicon and side primer '{name}' belongs to. Names need LEFT or RIGHT in them, "
                f"and --amplicon_pattern (now {pattern.pattern!r}) must capture the amplicon name."
            )
        primers.append({"name": name, "seq": seq, "amplicon": match.group(1), "side": side.group(0).upper()})
    if not primers:
        raise ProtocolError(f"No primers found in {path}")
    return primers


def reverse_complement(seq):
    return seq.upper().translate(_COMPLEMENT)[::-1]


TAIL = 5  # the primer's last bases at its 3' end, where a mismatch stops extension


def mismatches(sequence_codes, primer, tail_first=False):
    """Mismatches of ``primer`` (IUPAC) at every start position of an encoded sequence.

    Returns (total, tail): all mismatches, and those in the last ``TAIL`` bases of the primer's 3' end. ``primer`` is
    the sequence as it reads on the searched strand, so for a RIGHT primer (searched as its reverse complement) the
    3' end comes first: pass ``tail_first=True``. None if the sequence is shorter than the primer.
    """
    codes = _CODES[np.frombuffer(primer.encode(), dtype=np.uint8)]
    n = len(sequence_codes) - len(codes) + 1
    if n <= 0:
        return None
    total = np.zeros(n, dtype=np.int16)
    tail = np.zeros(n, dtype=np.int16)
    for j, code in enumerate(codes):
        miss = (sequence_codes[j : j + n] & code) == 0
        total += miss
        if (j < TAIL) if tail_first else (j >= len(codes) - TAIL):
            tail += miss
    return total, tail


def place_primers(primers, sequences, maps, max_mismatch=5, max_tail_mismatch=1, window=150, min_support=3):
    """Find every primer on every reference.

    A primer is searched for on both strands' worth of sequence (a RIGHT primer is its reverse complement on the given
    strand), allowing IUPAC codes, up to ``max_mismatch`` mismatches and at most ``max_tail_mismatch`` of them in the
    last 5 bases of the primer's 3' end. Genotype-specific primers find chance matches
    elsewhere, so a hit only counts if it lies within ``window`` bases (on the coordinate axis) of where the primer
    belongs: the position where at least ``min_support`` references agree on a close (<= 2 mismatches) match.

    Returns (sites, expected, unplaced): sites[ref][primer] = (start, end, mismatches), 0-based, end exclusive;
    expected[primer] = (coordinate position, supporting references); unplaced = primers with too little support.
    """
    encoded = {ref: _CODES[np.frombuffer(seq.encode(), dtype=np.uint8)] for ref, seq in sequences.items()}
    candidates = {}  # (primer, ref) -> (starts, mismatches)
    anchors = {}  # primer -> coordinate positions of close matches
    for p in primers:
        query = reverse_complement(p["seq"]) if p["side"] == "RIGHT" else p["seq"]
        half = len(query) // 2
        for ref, codes in encoded.items():
            found = mismatches(codes, query, tail_first=p["side"] == "RIGHT")
            if found is None:
                continue
            mm, tail = found
            starts = np.nonzero((mm <= max_mismatch) & (tail <= max_tail_mismatch))[0]
            if not len(starts):
                continue
            candidates[(p["name"], ref)] = (starts, mm[starts])
            best = starts[np.argmin(mm[starts])]
            centre = maps[ref][best + half] if ref in maps else -1
            if mm[best] <= 2 and centre >= 0:
                anchors.setdefault(p["name"], []).append(int(centre))
    expected, unplaced = {}, []
    for p in primers:
        points = np.sort(np.array(anchors.get(p["name"], []), dtype=np.int64))
        if not len(points):
            unplaced.append(p["name"])
            continue
        counts = [(np.abs(points - c) <= window).sum() for c in points]
        centre = points[int(np.argmax(counts))]
        cluster = points[np.abs(points - centre) <= window]
        if len(cluster) < min_support:
            unplaced.append(p["name"])
            continue
        expected[p["name"]] = (int(np.median(cluster)), int(len(cluster)))
    sites = {ref: {} for ref in sequences}
    for p in primers:
        if p["name"] not in expected:
            continue
        where, half = expected[p["name"]][0], len(p["seq"]) // 2
        for ref in sequences:
            found = candidates.get((p["name"], ref))
            if found is None or ref not in maps:
                continue
            starts, mm = found
            centre = maps[ref][np.minimum(starts + half, len(maps[ref]) - 1)]
            keep = (centre >= 0) & (np.abs(centre - where) <= window)
            if not keep.any():
                continue
            starts, mm, centre = starts[keep], mm[keep], centre[keep]
            pick = min(range(len(starts)), key=lambda i: (mm[i], abs(centre[i] - where)))
            sites[ref][p["name"]] = (int(starts[pick]), int(starts[pick]) + len(p["seq"]), int(mm[pick]))
    return sites, expected, unplaced


def build_amplicons(primers, sites, min_length=150, max_length=4000):
    """Per reference, one span per amplicon: its best-matching LEFT primer to its best-matching RIGHT primer.

    Returns {ref: [(name, start, end, left_primer, right_primer, left_mismatches, right_mismatches)]} with 0-based,
    end-exclusive spans, amplicons in natural order. Amplicons whose primers were not both found are left out.
    """
    by_amplicon = {}
    for p in primers:
        by_amplicon.setdefault(p["amplicon"], {"LEFT": [], "RIGHT": []})[p["side"]].append(p["name"])
    result = {}
    for ref, found in sites.items():
        spans = []
        for name in sorted(by_amplicon, key=_natural):
            pairs = []
            for left in by_amplicon[name]["LEFT"]:
                for right in by_amplicon[name]["RIGHT"]:
                    if left not in found or right not in found:
                        continue
                    (ls, _, lm), (_, re_, rm) = found[left], found[right]
                    if min_length <= re_ - ls <= max_length:
                        pairs.append((lm + rm, re_ - ls, name, ls, re_, left, right, lm, rm))
            if pairs:
                _, _, name, start, end, left, right, lm, rm = min(pairs)
                spans.append((name, start, end, left, right, lm, rm))
        result[ref] = spans
    return result


def lift_bed_amplicons(coordinate_amplicons, maps, min_fraction=0.7):
    """Carry amplicons given on the coordinate reference onto every reference through its coordinate map."""
    result = {}
    for ref, forward in maps.items():
        spans = []
        for name, start, end in coordinate_amplicons:
            inside = np.nonzero((forward >= start) & (forward < end))[0]
            if len(inside) >= min_fraction * (end - start):
                spans.append((name, int(inside.min()), int(inside.max()) + 1, "", "", 0, 0))
        result[ref] = spans
    return result


# ---- genes and regions ----------------------------------------------------------------------------------------------


def lift_features(features, maps, lengths):
    """Lift {name: (start, end)} (1-based, inclusive, on the coordinate reference) onto every reference.

    Returns {ref: {name: (start, end, mapped_bp, expected_bp, status)}} (1-based, inclusive); status is ``ok``,
    ``odd_span`` (mapped, but the lifted span is far from the expected length), ``partial`` (50-90% of the feature
    has a counterpart) or ``absent``.
    """
    lifted = {}
    for ref, forward in maps.items():
        row = {}
        for name, (start, end) in features.items():
            expected = end - start + 1
            inside = np.nonzero((forward >= start - 1) & (forward <= end - 1))[0]
            mapped = len(inside)
            if mapped < 0.5 * expected:
                row[name] = (0, 0, mapped, expected, "absent")
                continue
            a, b = int(inside.min()) + 1, int(inside.max()) + 1
            span = b - a + 1
            if mapped < 0.9 * expected:
                status = "partial"
            elif not 0.8 * expected <= span <= 1.25 * expected:
                status = "odd_span"
            else:
                status = "ok"
            row[name] = (a, b, mapped, expected, status)
        lifted[ref] = row
    return lifted


def parse_region(text):
    match = re.fullmatch(r"([^=\s]+)=(\d+)-(\d+)", text.strip())
    if not match or int(match.group(2)) > int(match.group(3)):
        raise ProtocolError(f"--region must look like NAME=START-END (1-based, on the coordinate reference), got '{text}'")
    return match.group(1), (int(match.group(2)), int(match.group(3)))


# ---- building -------------------------------------------------------------------------------------------------------


def _write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _copy_in(src, outdir, name):
    dest = Path(outdir) / name
    if Path(src).resolve() != dest.resolve():
        shutil.copyfile(src, dest)
    return dest


def build_protocol(outdir, possible_references, reference_fasta=None, gff=None, regions=None, primer_bed=None,
                   primer_fasta=None, amplicon_pattern=None, panel_alignment=None, name=None,
                   group_pattern=DEFAULT_GROUP_PATTERN, threads=4, max_mismatch=5, say=print,
                   coordinate_ref=None, max_tail_mismatch=1):  # fmt: skip
    """Build a protocol folder. Returns the QC summary (also stored in protocol.json). ``say`` receives progress lines."""
    if primer_bed and primer_fasta:
        raise ProtocolError("Give the primer scheme as --primer_bed or --primer_fasta, not both")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    regions = dict(regions or {})

    panel = read_fasta(possible_references)
    texts = [_fasta_text(possible_references)]
    ids = {rid for rid, _, _ in panel}
    chosen = coordinate_ref
    coordinate_ref = panel[0][0]
    if reference_fasta:
        extra = read_fasta(reference_fasta)
        coordinate_ref = extra[0][0]
        if [r for r in extra if r[0] not in ids]:
            texts.append(_fasta_text(reference_fasta))
    if chosen:
        if chosen not in ids:
            raise ProtocolError(f"--coordinate_ref {chosen} is not in the panel ({possible_references})")
        coordinate_ref = chosen
    panel_text = "".join(texts)
    (outdir / "panel.fasta").write_text(panel_text)
    sequences = read_sequences(outdir / "panel.fasta")
    if coordinate_ref not in sequences:
        raise ProtocolError(f"The coordinate reference {coordinate_ref} is not in the panel")
    lengths = {k: len(v) for k, v in sequences.items()}
    group = re.compile(group_pattern)

    genes = parse_gff_genes(gff, coordinate_ref) if gff else {}
    for rname, (start, end) in regions.items():
        if end > lengths[coordinate_ref]:
            raise ProtocolError(f"Region {rname} ends at {end}, past the end of {coordinate_ref} ({lengths[coordinate_ref]})")

    if panel_alignment:
        say(f"Using the ready-made alignment {panel_alignment}")
        aligned = parse_alignment(Path(panel_alignment).read_text().splitlines())
        mafft_version = "not used"
    else:
        say(f"Aligning {len(sequences)} sequence(s) to {coordinate_ref} with MAFFT ...")
        aligned = run_mafft(coordinate_ref, sequences, threads)
        mafft_version = _tool_version("mafft", "--version")
    maps = build_maps(aligned, coordinate_ref, sequences)

    placed = {ref: float((forward >= 0).mean()) for ref, forward in maps.items()}
    poor = sorted(ref for ref, frac in placed.items() if frac < POOR_ALIGNMENT)
    good = {ref: forward for ref, forward in maps.items() if ref not in poor}

    sites, expected, unplaced, amplicons = {}, {}, [], {}
    primers = []
    if primer_fasta:
        primers = read_primers(primer_fasta, amplicon_pattern)
        say(f"Finding {len(primers)} primer(s) on {len(good)} reference(s) ...")
        sites, expected, unplaced = place_primers(primers, sequences, good, max_mismatch, max_tail_mismatch)
        amplicons = build_amplicons(primers, {ref: s for ref, s in sites.items() if ref in good})
    elif primer_bed:
        coordinate_amplicons = parse_bed_amplicons(primer_bed, coordinate_ref)
        amplicons = lift_bed_amplicons(coordinate_amplicons, good)
    amplicon_names = sorted({a[0] for spans in amplicons.values() for a in spans}, key=_natural)

    features = {f"gene:{k}": (v["start"], v["end"]) for k, v in genes.items()}
    features.update({f"region:{k}": v for k, v in regions.items()})
    lifted = lift_features(features, good, lengths)

    # ---- files
    np.savez_compressed(outdir / "coordinate_maps.npz", ids=np.array(list(good)),
                        **{f"m{i}": good[ref] for i, ref in enumerate(good)})  # fmt: skip
    _write_tsv(outdir / "amplicons.tsv",
               ["reference", "amplicon", "start", "end", "length", "left_primer", "right_primer", "left_mismatches", "right_mismatches"],
               [(ref, n, s + 1, e, e - s, lp, rp, lm, rm) for ref, spans in amplicons.items() for n, s, e, lp, rp, lm, rm in spans])  # fmt: skip
    if primers:
        by_name = {p["name"]: p for p in primers}
        _write_tsv(outdir / "primer_sites.tsv",
                   ["reference", "primer", "amplicon", "side", "start", "end", "mismatches", "sequence"],
                   [(ref, pn, by_name[pn]["amplicon"], by_name[pn]["side"], s + 1, e, mm, by_name[pn]["seq"])
                    for ref, found in sites.items() for pn, (s, e, mm) in found.items()])  # fmt: skip
    if features:
        _write_tsv(outdir / "features.tsv",
                   ["reference", "feature", "start", "end", "mapped_bp", "expected_bp", "status"],
                   [(ref, fname, *row[fname][:2], *row[fname][2:]) for ref, row in lifted.items() for fname in row])  # fmt: skip
    inputs = {}
    for label, src in (("genes", gff), ("primers", primer_fasta or primer_bed)):
        if src:
            dest = _copy_in(src, outdir, {"genes": "genes.gff", "primers": "primers" + Path(src).suffix}[label])
            inputs[label] = dest.name

    qc = _qc_summary(lengths, coordinate_ref, maps, placed, poor, amplicon_names, amplicons, expected, unplaced,
                     lifted, features, group)  # fmt: skip
    (outdir / "qc.txt").write_text(_qc_text(qc))
    meta = {
        "format": FORMAT,
        "name": name or outdir.resolve().name,
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "nanosentinel": __version__,
        "mafft": mafft_version,
        "coordinate_ref": coordinate_ref,
        "coordinate_length": lengths[coordinate_ref],
        "group_pattern": group_pattern,
        "genes": genes,
        "regions": {k: {"start": s, "end": e} for k, (s, e) in regions.items()},
        "amplicon_names": amplicon_names,
        "amplicon_pattern": amplicon_pattern,
        "primer_matching": {"max_mismatch": max_mismatch, "max_3prime_mismatch": max_tail_mismatch} if primers else None,
        "inputs": inputs,
        "unaligned_references": poor,
        "qc": qc,
        "files": {f: sha1_file(outdir / f) for f in ("panel.fasta", "coordinate_maps.npz")},
    }
    (outdir / "protocol.json").write_text(json.dumps(meta, indent=2, default=int) + "\n")
    return qc


def _tool_version(name, flag):
    try:
        done = subprocess.run([tools.find(name), flag], capture_output=True, text=True)
        return (done.stdout + done.stderr).strip().splitlines()[0]
    except Exception:
        return "unknown"


def _qc_summary(lengths, coordinate_ref, maps, placed, poor, amplicon_names, amplicons, expected, unplaced,
                lifted, features, group):  # fmt: skip
    def group_of(ref):
        match = group.match(ref)
        return match.group(1) if match else ref

    refs = list(maps)
    qc = {
        "references": len(refs),
        "aligned_references": len(refs) - len(poor),
        "poorly_aligned": poor,
        "placed_fraction_min": round(min(placed.values()), 3) if placed else None,
        "unplaced_primers": unplaced,
        "primer_support": {k: {"coordinate": v[0], "references": v[1]} for k, v in expected.items()},
    }
    if amplicon_names:
        per_amplicon, missing_groups = {}, {}
        for name in amplicon_names:
            have = [ref for ref in refs if ref not in poor and any(a[0] == name for a in amplicons.get(ref, []))]
            per_amplicon[name] = len(have)
            groups = {}
            for ref in refs:
                if ref in poor:
                    continue
                groups.setdefault(group_of(ref), [0, 0])[1] += 1
                groups[group_of(ref)][0] += ref in have
            missing_groups[name] = sorted(g for g, (h, n) in groups.items() if h == 0)
        qc["amplicon_references"] = per_amplicon
        qc["amplicon_missing_groups"] = missing_groups
        qc["references_with_all_amplicons"] = sum(
            all(any(a[0] == name for a in amplicons.get(ref, [])) for name in amplicon_names) for ref in refs if ref not in poor
        )
    if features:
        stats = {}
        for fname in features:
            rows = [lifted[ref][fname] for ref in lifted if fname in lifted[ref]]
            status = {}
            for row in rows:
                status[row[4]] = status.get(row[4], 0) + 1
            spans = [r[1] - r[0] + 1 for r in rows if r[4] in ("ok", "odd_span")]
            stats[fname] = {
                "expected_bp": features[fname][1] - features[fname][0] + 1,
                "status": status,
                "span_min": min(spans) if spans else None,
                "span_max": max(spans) if spans else None,
            }
        qc["features"] = stats
    return qc


def _qc_text(qc):
    lines = [f"References: {qc['references']} ({qc['aligned_references']} placed on the shared axis)"]
    if qc["poorly_aligned"]:
        lines.append(f"  poorly aligned (left off the shared axis): {', '.join(qc['poorly_aligned'])}")
    if "features" in qc:
        lines.append("")
        lines.append("Genes and regions lifted onto every reference (status counts; span in bases):")
        for fname, st in qc["features"].items():
            counts = ", ".join(f"{k} {v}" for k, v in sorted(st["status"].items()))
            span = f"{st['span_min']:,}-{st['span_max']:,}" if st["span_min"] else "-"
            lines.append(f"  {fname:<18} expected {st['expected_bp']:>6,}  span {span:<13} {counts}")
    if "amplicon_references" in qc:
        lines.append("")
        lines.append(f"Amplicons: {qc['references_with_all_amplicons']} of {qc['aligned_references']} references have all "
                     f"{len(qc['amplicon_references'])}")  # fmt: skip
        for name, n in qc["amplicon_references"].items():
            groups = qc["amplicon_missing_groups"][name]
            note = f"  (no site on any {', '.join(groups[:12])}{' ...' if len(groups) > 12 else ''})" if groups else ""
            lines.append(f"  {name:<10} on {n} references{note}")
    if qc["unplaced_primers"]:
        lines.append("")
        lines.append("Primers with too few close matches to place (check them): " + ", ".join(qc["unplaced_primers"]))
    return "\n".join(lines) + "\n"


# ---- loading --------------------------------------------------------------------------------------------------------


def resolve_protocol(arg):
    """A protocol folder given as a path, or the name of one bundled with nanosentinel (protocols/<name>)."""
    path = Path(arg).expanduser()
    if (path / "protocol.json").is_file():
        return path
    bundled = BUNDLED / arg
    if (bundled / "protocol.json").is_file():
        return bundled
    names = sorted(p.parent.name for p in BUNDLED.glob("*/protocol.json")) if BUNDLED.is_dir() else []
    hint = f" Bundled protocols: {', '.join(names)}." if names else ""
    raise ProtocolError(f"No protocol found at '{arg}' (a folder containing protocol.json).{hint}")


def load_protocol(arg, outdir):
    """Load a protocol folder into a ReferenceInfo, indexing the panel under ``outdir``/reference."""
    path = resolve_protocol(arg)
    meta = json.loads((path / "protocol.json").read_text())
    if meta.get("format") != FORMAT:
        raise ProtocolError(f"{path} has protocol format {meta.get('format')}; this nanosentinel reads format {FORMAT}. Rebuild it with prep-protocol.")
    for fname, want in meta["files"].items():
        if not (path / fname).is_file() or sha1_file(path / fname) != want:
            raise ProtocolError(f"{path / fname} has changed since the protocol was built. Rebuild it with prep-protocol.")
    panel_text = (path / "panel.fasta").read_text()
    lengths = {rid: n for rid, _, n in read_fasta(path / "panel.fasta")}
    maps = {}
    with np.load(path / "coordinate_maps.npz") as data:
        for i, ref in enumerate(data["ids"]):
            maps[str(ref)] = data[f"m{i}"]
    amplicons = {}
    if (path / "amplicons.tsv").is_file():
        with open(path / "amplicons.tsv", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                amplicons.setdefault(row["reference"], []).append((row["amplicon"], int(row["start"]) - 1, int(row["end"])))
    index, panel_hash = build_index(panel_text, Path(outdir) / "reference")
    return ReferenceInfo(
        lengths, meta["coordinate_ref"], meta["genes"], meta["regions"], index, panel_hash, maps, amplicons,
        meta["amplicon_names"], meta["group_pattern"], meta["name"],
    )  # fmt: skip
