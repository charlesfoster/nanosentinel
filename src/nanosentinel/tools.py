"""Locate the external programs (minimap2) that live next to the running Python (the pixi environment)."""

import os
import shutil
import sys

from .errors import MonitorError


def _path():
    return os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")


def find(name):
    exe = shutil.which(name, path=_path())
    if not exe:
        raise MonitorError(
            f"Required tool not found: {name}. Run nanosentinel through its launcher (bin/nanosentinel) "
            "or inside the pixi environment."
        )
    return exe
