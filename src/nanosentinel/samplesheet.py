"""Read the samplesheet (CSV/TSV) with the columns monitoring needs.

Required: ``sample``, ``barcode``.
Optional: ``run_name`` (selects rows when a sheet spans runs) and ``neg_control``.
Every other column (``long_reads``, ``short_reads_1`` ...) is carried through untouched.
"""

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

from .errors import SheetError

BARCODE_RE = re.compile(r"^(?:barcode|bc|nb)?[\s_-]*0*(\d{1,3})$", re.IGNORECASE)
SAMPLE_RE = re.compile(r"^[^\s/\\=]+$")
TRUE = {"true", "t", "yes", "y", "1"}
FALSE = {"", "false", "f", "no", "n", "0"}


@dataclass
class Row:
    sample: str
    barcode: str  # normalised, e.g. "barcode01"
    run_name: str
    neg_control: bool
    fields: dict = field(default_factory=dict)  # every original column, as read


def normalise_barcode(value: str) -> str:
    """BC01, bc1, NB01, barcode1 and 01 all become 'barcode01'."""
    match = BARCODE_RE.match(value.strip())
    if not match or int(match.group(1)) == 0:
        raise SheetError(f"Cannot interpret barcode '{value}' (expected e.g. BC01 or barcode01)")
    return f"barcode{int(match.group(1)):02d}"


def _parse_bool(value: str, where: str) -> bool:
    v = value.strip().lower()
    if v in TRUE:
        return True
    if v in FALSE:
        return False
    raise SheetError(f"{where}: neg_control must be true/false (or yes/no, 1/0), got '{value}'")


def read_samplesheet(path):
    """Return (header, rows)."""
    path = Path(path)
    if not path.is_file():
        raise SheetError(f"Samplesheet not found: {path}")
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        header = [h.strip() for h in (reader.fieldnames or [])]
        if len(set(header)) != len(header):
            raise SheetError("Samplesheet has duplicate column names")
        missing = [c for c in ("sample", "barcode") if c not in header]
        if missing:
            raise SheetError(
                f"Samplesheet is missing required column(s): {', '.join(missing)} "
                f"(found: {', '.join(header) or 'none'})"
            )
        rows = []
        for n, raw in enumerate(reader, start=2):
            fields = {k.strip(): (v or "").strip() for k, v in raw.items() if k is not None}
            if not any(fields.values()):
                continue
            where = f"{path.name} line {n}"
            sample = fields["sample"]
            if not SAMPLE_RE.match(sample):
                raise SheetError(f"{where}: sample name '{sample}' must be non-empty with no spaces, '/', '\\' or '='")
            try:
                barcode = normalise_barcode(fields["barcode"])
            except SheetError as exc:
                raise SheetError(f"{where}: {exc}") from None
            rows.append(
                Row(
                    sample=sample,
                    barcode=barcode,
                    run_name=fields.get("run_name", ""),
                    neg_control=_parse_bool(fields.get("neg_control", ""), where),
                    fields=fields,
                )
            )
    if not rows:
        raise SheetError(f"{path.name} contains no samples")
    return header, rows


def select_run(rows, run_name=None):
    """Pick the rows belonging to one run. Returns (run_name or None, rows)."""
    names = sorted({r.run_name for r in rows if r.run_name})
    if run_name is None:
        if len(names) == 1:
            run_name = names[0]
        elif len(names) > 1:
            raise SheetError(f"Samplesheet spans several runs ({', '.join(names)}); choose one with --run_name")
    if run_name:
        rows = [r for r in rows if r.run_name in ("", run_name)]
        if not rows:
            raise SheetError(f"No rows in the samplesheet for run_name '{run_name}'")
    seen = {}
    for r in rows:
        if r.barcode in seen:
            raise SheetError(
                f"{r.barcode} is listed twice for the same run ({seen[r.barcode]} and {r.sample}); "
                "each barcode can appear only once per run"
            )
        seen[r.barcode] = r.sample
    return run_name or None, rows
