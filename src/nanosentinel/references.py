"""Reference panel, coordinate reference, primer scheme and genes: parsing and minimap2 index preparation."""

import hashlib
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from . import tools
from .errors import ProtocolError


@dataclass
class ReferenceInfo:
    """What the monitor needs to know about the references."""

    lengths: dict  # reference id -> length
    coordinate_ref: str  # reference the coverage plot (and primer scheme) is drawn against
    amplicons: list = field(default_factory=list)  # [[start, end], ...] on coordinate_ref (0-based, end exclusive)
    genes: dict = field(default_factory=dict)  # name -> {start, end, strand} on coordinate_ref (1-based, inclusive)
    index: Path = None  # minimap2 index of the panel
    panel_hash: str = ""  # identifies the panel, so stale annotations are never reused


def read_fasta(path):
    """Return a list of (id, header_line_without_gt, sequence_length)."""
    path = Path(path)
    if not path.is_file():
        raise ProtocolError(f"FASTA not found: {path}")
    records, header, length = [], None, 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header.split()[0], header, length))
                header, length = line[1:].strip(), 0
            elif header is not None:
                length += len(line.strip())
    if header is not None:
        records.append((header.split()[0], header, length))
    if not records:
        raise ProtocolError(f"No sequences found in {path}")
    ids = [r[0] for r in records]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        raise ProtocolError(f"{path.name} has duplicated sequence ids: {', '.join(dupes[:5])}")
    return records


def _fasta_text(path):
    return Path(path).read_text().rstrip("\n") + "\n"


def parse_bed_amplicons(path, ref_id):
    """Group ARTIC-style primers (<name>_LEFT / <name>_RIGHT[_alt]) into amplicon spans."""
    spans = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            match = re.match(r"^(.*?)_(?:LEFT|RIGHT)", cols[3])
            if not match:
                continue
            start, end = int(cols[1]), int(cols[2])
            if cols[0] != ref_id:
                continue
            lo, hi = spans.get(match.group(1), (start, end))
            spans[match.group(1)] = (min(lo, start), max(hi, end))
    if not spans:
        raise ProtocolError(
            f"No amplicons on '{ref_id}' found in {Path(path).name}. Primer names must look like "
            "<amplicon>_LEFT / <amplicon>_RIGHT and the BED chromosome must match the coordinate reference."
        )
    return [list(s) for s in sorted(spans.values())]


def parse_gff_genes(path, ref_id):
    """Genes (or CDS if there are no 'gene' features) for the coordinate reference, 1-based inclusive."""
    found = {"gene": [], "CDS": []}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[0] != ref_id or cols[2] not in found:
                continue
            attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
            name = attrs.get("Name") or attrs.get("gene") or attrs.get("product") or attrs.get("ID") or cols[2]
            found[cols[2]].append((name, int(cols[3]), int(cols[4]), -1 if cols[6] == "-" else 1))
    features = found["gene"] or found["CDS"]
    if not features:
        raise ProtocolError(f"No gene/CDS features on '{ref_id}' in {Path(path).name}")
    genes = {}
    for name, start, end, strand in features:
        key, n = name, 2
        while key in genes:
            key, n = f"{name}_{n}", n + 1
        genes[key] = {"start": start, "end": end, "strand": strand}
    return genes


def amplicon_mask_intervals(amplicons):
    """Merge amplicon spans into non-overlapping intervals."""
    merged = []
    for start, end in sorted(amplicons):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def prepare_references(outdir, possible_references=None, reference_fasta=None, primer_bed=None, gff=None):
    """Write <outdir>/reference/panel.fasta, index it with minimap2 and return a ReferenceInfo."""
    if not possible_references and not reference_fasta:
        raise ProtocolError("Give --possible_references (a panel) and/or --reference_fasta")
    refdir = Path(outdir) / "reference"
    refdir.mkdir(parents=True, exist_ok=True)

    panel = read_fasta(possible_references) if possible_references else read_fasta(reference_fasta)
    texts = [_fasta_text(possible_references or reference_fasta)]
    lengths = {rid: n for rid, _, n in panel}
    coordinate_ref = panel[0][0]
    if reference_fasta:
        extra = read_fasta(reference_fasta)
        coordinate_ref = extra[0][0]
        if possible_references and [r for r in extra if r[0] not in lengths]:
            texts.append(_fasta_text(reference_fasta))  # the coordinate reference must be in the panel
            lengths.update({rid: n for rid, _, n in extra if rid not in lengths})
    text = "".join(texts)
    panel_hash = hashlib.sha1(text.encode()).hexdigest()[:12]
    fasta, index, stamp = refdir / "panel.fasta", refdir / "panel.mmi", refdir / "panel.sha1"
    if not (index.exists() and stamp.exists() and stamp.read_text() == panel_hash):
        fasta.write_text(text)
        done = subprocess.run(
            [tools.find("minimap2"), "-v", "1", "-x", "map-ont", "-d", str(index), str(fasta)],
            capture_output=True, text=True,
        )  # fmt: skip
        if done.returncode != 0:
            raise ProtocolError(f"minimap2 could not index the reference panel:\n{done.stderr.strip()[-500:]}")
        stamp.write_text(panel_hash)

    genes = parse_gff_genes(gff, coordinate_ref) if gff else {}
    amplicons = parse_bed_amplicons(primer_bed, coordinate_ref) if primer_bed else []
    return ReferenceInfo(lengths, coordinate_ref, amplicons, genes, index, panel_hash)
