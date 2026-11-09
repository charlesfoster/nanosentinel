"""Start / stop RAMPART with all of its (very chatty) output going to a log file."""

import atexit
import os
import shutil
import signal
import socket
import subprocess
import sys
import time
import webbrowser
from pathlib import Path

from .errors import MonitorError


def _env():
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + env.get("PATH", "")
    return env


def check_tools():
    path = _env()["PATH"]
    missing = [t for t in ("rampart", "snakemake", "minimap2") if not shutil.which(t, path=path)]
    if missing:
        raise MonitorError(
            f"Required tool(s) not found: {', '.join(missing)}. Run nanosentinel through its launcher (bin/nanosentinel) "
            "or inside the pixi environment."
        )


def check_path(fastq_pass):
    """RAMPART's file watcher ignores any path containing a hidden (dot-prefixed) folder, silently."""
    hidden = [part for part in Path(fastq_pass).parts if part.startswith(".") and part not in (".", "..")]
    if hidden:
        raise MonitorError(
            f"RAMPART ignores reads under hidden folders, and {fastq_pass} contains '{hidden[0]}'. "
            "Point --run_dir / --minknow_data at a path with no dot-prefixed folder."
        )


def check_ports(ports):
    for port in ports:
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                raise MonitorError(f"Port {port} is already in use (is another RAMPART running?). Use --rampart_port.")


def launch(protocol_dir, fastq_pass, barcode_names, title, ports, workdir, logfile, clear=False):
    """Start RAMPART with ``workdir`` as its working directory.

    RAMPART 1.2.0 silently ignores --annotatedPath (its code reads ``args.annotatedDir``), so annotations
    always land in <cwd>/annotations; ``workdir`` is therefore where they, and Snakemake's state, end up.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        shutil.which("rampart", path=_env()["PATH"]),
        "--protocol", str(protocol_dir),
        "--ports", str(ports[0]), str(ports[1]),
        "--basecalledPath", str(fastq_pass),
    ]  # fmt: skip
    if title:
        cmd += ["--title", title]
    if barcode_names:
        cmd += ["--barcodeNames", *barcode_names]
    if clear:
        cmd.append("--clearAnnotated")
    log = open(logfile, "ab")
    log.write(("\n$ " + " ".join(cmd) + "\n").encode())
    log.flush()
    proc = subprocess.Popen(
        cmd, cwd=workdir, env=_env(), stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
        start_new_session=True,  # own process group so we can stop RAMPART and its Snakemake children together
    )  # fmt: skip
    atexit.register(stop, proc)  # safety net: no exit path (even an exception) may leave RAMPART running
    return proc


def wait_until_listening(proc, port, timeout=45):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if proc.poll() is not None:
            return False
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return True
        time.sleep(0.5)
    return False


def open_browser(port):
    try:
        webbrowser.open(f"http://localhost:{port}", new=1)
    except Exception:
        pass


def stop(proc, grace=10):
    """Stop RAMPART and its children. Repeated Ctrl+C escalates to SIGKILL instead of leaving them running."""
    if proc.poll() is not None:
        return
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        proc.wait(timeout=grace)
        return
    except (ProcessLookupError, subprocess.TimeoutExpired, KeyboardInterrupt):
        pass
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5)
    except (subprocess.TimeoutExpired, KeyboardInterrupt):
        pass


def tail(path, n=8):
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-n:])
