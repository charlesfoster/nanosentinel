"""Locate the MinKNOW ``fastq_pass`` directory for a run."""

import os
from pathlib import Path

from .errors import RunDirError, RunNotReady

DEFAULT_MINKNOW_DATA = "/var/lib/minknow/data"
# Never descend into these: they hold the bulk of a run's files.
_PRUNE = {"fastq_pass", "fastq_fail", "pod5", "pod5_pass", "pod5_fail", "fast5", "fast5_pass", "fast5_fail"}


def find_fastq_pass(run_dir=None, minknow_data=DEFAULT_MINKNOW_DATA, run_name=None) -> Path:
    if run_dir:
        p = Path(run_dir).expanduser()
        if p.name == "fastq_pass" and p.is_dir():
            return p.resolve()
        if (p / "fastq_pass").is_dir():
            return (p / "fastq_pass").resolve()
        raise RunNotReady(f"{p if p.name == 'fastq_pass' else p / 'fastq_pass'} does not exist")
    if not run_name:
        raise RunDirError("Cannot find the run: give --run_name (or a run_name column in the samplesheet) or --run_dir")
    root = Path(minknow_data).expanduser()
    if not root.is_dir():
        raise RunDirError(f"MinKNOW data directory not found: {root} (set --minknow_data, or use --run_dir)")
    matches = []
    for dirpath, dirnames, _ in os.walk(root):
        if "fastq_pass" in dirnames:
            candidate = Path(dirpath) / "fastq_pass"
            if run_name in candidate.relative_to(root).parts[:-1]:
                matches.append(candidate)
        dirnames[:] = [d for d in dirnames if d not in _PRUNE]
    if not matches:
        raise RunNotReady(f"no fastq_pass folder found under {root} for run '{run_name}'")
    if len(matches) > 1:
        listing = "\n  ".join(str(m) for m in sorted(matches))
        raise RunDirError(f"Run '{run_name}' matches several folders; pass the exact one with --run_dir:\n  {listing}")
    return matches[0].resolve()
