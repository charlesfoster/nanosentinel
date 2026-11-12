"""Reference panel, coordinate reference, primer scheme and genes: parsing, minimap2 index and coordinate projection.

Every reference can carry a *coordinate map*: for each of its bases, the position of the aligned base on the coordinate
reference (-1 where it has no counterpart). That lets depth measured on any reference be drawn on one shared axis, so
genes, regions and amplicons mean the same thing for samples of different genotypes.
"""

import hashlib
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import tools
from .errors import MonitorError, ProtocolError

DEFAULT_GROUP_PATTERN = r"^([^_]+)_"  # "1a_AF009606.1" -> "1a": references sharing a prefix are siblings


@dataclass
class ReferenceInfo:
    """What the monitor needs to know about the references."""

    lengths: dict  # reference id -> length
    coordinate_ref: str  # reference the shared axis (and the gene annotation) is drawn against
    genes: dict = field(default_factory=dict)  # name -> {start, end, strand} on coordinate_ref (1-based, inclusive)
    regions: dict = field(default_factory=dict)  # name -> {start, end} on coordinate_ref (1-based, inclusive)
    index: Path = None  # minimap2 index of the panel
    panel_hash: str = ""  # identifies the panel, so stale annotations are never reused
    maps: dict = field(default_factory=dict)  # reference id -> int32 array: base -> coordinate position (-1 = none)
    amplicons: dict = field(default_factory=dict)  # reference id -> [(name, start, end)], 0-based, end exclusive
    amplicon_names: list = field(default_factory=list)  # every amplicon name, in order (the heatmap columns)
    group_pattern: str = DEFAULT_GROUP_PATTERN
    protocol: str = ""  # protocol name, or "" when running from plain files
    _inverse: dict = field(default_factory=dict, repr=False)

    @property
    def coordinate_length(self):
        return self.lengths[self.coordinate_ref]

    def group(self, ref):
        """The sibling group of a reference (for example its subtype)."""
        match = re.match(self.group_pattern, ref)
        return match.group(1) if match else ref

    def annotated(self, ref):
        return ref in self.maps

    def inverse(self, ref):
        """Array over the coordinate reference: position on ``ref`` (or -1). None if ``ref`` has no map."""
        if ref not in self.maps:
            return None
        if ref not in self._inverse:
            forward = self.maps[ref]
            inv = np.full(self.coordinate_length, -1, dtype=np.int32)
            have = np.nonzero(forward >= 0)[0]
            inv[forward[have]] = have
            self._inverse[ref] = inv
        return self._inverse[ref]

    def project(self, ref, depth):
        """Depth on ``ref`` laid out on the coordinate axis: (depth, valid) or None when ``ref`` cannot be projected."""
        inv = self.inverse(ref)
        if inv is None or depth is None:
            return None
        valid = inv >= 0
        out = np.zeros(len(inv), dtype=np.float64)
        out[valid] = depth[inv[valid]]
        return out, valid

    def to_coordinate(self, ref, start, end):
        """Lift a [start, end) interval on ``ref`` to the coordinate axis; None if none of it can be placed."""
        forward = self.maps.get(ref)
        if forward is None:
            return None
        placed = forward[start:end]
        placed = placed[placed >= 0]
        return (int(placed.min()), int(placed.max()) + 1) if len(placed) else None

    def region_length(self, name):
        region = self.regions[name]
        return region["end"] - region["start"] + 1

    def region_mask(self, ref, name):
        """Bases of ``ref`` that line up with a named region of the coordinate reference (None if it has no map)."""
        forward = self.maps.get(ref)
        if forward is None or name not in self.regions:
            return None
        region = self.regions[name]
        return (forward >= region["start"] - 1) & (forward <= region["end"] - 1)

    def mask(self, ref):
        """Bases of ``ref`` inside its amplicons, or None if the reference has none."""
        spans = self.amplicons.get(ref)
        if not spans:
            return None
        mask = np.zeros(self.lengths[ref], dtype=bool)
        for start, end in amplicon_mask_intervals([[s, e] for _, s, e in spans]):
            mask[start:end] = True
        return mask



def check_region(refinfo, name):
    """The region's name, or a MonitorError saying which regions this protocol has."""
    if name not in refinfo.regions:
        have = ", ".join(refinfo.regions) or "none"
        raise MonitorError(f"--target_region {name}: the protocol has no region of that name (regions: {have}). "
                           "Regions are defined when the protocol is built, with prep-protocol --region NAME=START-END.")  # fmt: skip
    return name


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


def read_sequences(path):
    """{id: sequence} (upper case, U read as T), in file order."""
    sequences, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name is not None:
                sequences[name].append(line.strip())
    return {k: "".join(v).upper().replace("U", "T") for k, v in sequences.items()}


def _fasta_text(path):
    return Path(path).read_text().rstrip("\n") + "\n"


def parse_bed_amplicons(path, ref_id):
    """Group ARTIC-style primers (<name>_LEFT / <name>_RIGHT[_alt]) into [(name, start, end)] spans on ``ref_id``."""
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
    return [(name, lo, hi) for name, (lo, hi) in sorted(spans.items(), key=lambda kv: kv[1])]


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


def build_index(fasta_text, refdir):
    """Write <refdir>/panel.fasta and its minimap2 index (reused while the panel is unchanged). Returns (index, hash)."""
    refdir = Path(refdir)
    refdir.mkdir(parents=True, exist_ok=True)
    panel_hash = hashlib.sha1(fasta_text.encode()).hexdigest()[:12]
    fasta, index, stamp = refdir / "panel.fasta", refdir / "panel.mmi", refdir / "panel.sha1"
    if not (index.exists() and stamp.exists() and stamp.read_text() == panel_hash):
        fasta.write_text(fasta_text)
        done = subprocess.run(
            [tools.find("minimap2"), "-v", "1", "-x", "map-ont", "-d", str(index), str(fasta)],
            capture_output=True, text=True,
        )  # fmt: skip
        if done.returncode != 0:
            raise ProtocolError(f"minimap2 could not index the reference panel:\n{done.stderr.strip()[-500:]}")
        stamp.write_text(panel_hash)
    return index, panel_hash


def prepare_references(outdir, possible_references=None, reference_fasta=None, primer_bed=None, gff=None):
    """Quick mode (plain files): index the panel and annotate the coordinate reference only."""
    if not possible_references and not reference_fasta:
        raise ProtocolError("Give --possible_references (a panel) and/or --reference_fasta")
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
    index, panel_hash = build_index("".join(texts), Path(outdir) / "reference")

    genes = parse_gff_genes(gff, coordinate_ref) if gff else {}
    amplicons = parse_bed_amplicons(primer_bed, coordinate_ref) if primer_bed else []
    return ReferenceInfo(
        lengths, coordinate_ref, genes, {}, index, panel_hash,
        maps={coordinate_ref: np.arange(lengths[coordinate_ref], dtype=np.int32)},
        amplicons={coordinate_ref: amplicons} if amplicons else {},
        amplicon_names=[name for name, _, _ in amplicons],
    )  # fmt: skip
