"""The live monitoring loop: annotate reads, track coverage, serve the dashboard, write reports."""

import json
import shutil
import signal
import socket
import statistics
import sys
import threading
import time
import urllib.request
import webbrowser

import numpy as np

from . import dashboard, report
from .annotate import Annotator
from .coverage import AnnotationReader, SampleTracker
from .errors import MonitorError
from .references import amplicon_mask_intervals
from .report import hms


def build_trackers(rows, refinfo, target_depth, target_fraction, denominator="genome"):
    masks = {}
    if denominator == "amplicons":
        if not refinfo.amplicons:
            raise MonitorError("--denominator amplicons needs --primer_bed")
        mask = np.zeros(refinfo.lengths[refinfo.coordinate_ref], dtype=bool)
        for start, end in amplicon_mask_intervals(refinfo.amplicons):
            mask[start:end] = True
        masks[refinfo.coordinate_ref] = mask
    return {r.barcode: SampleTracker(refinfo.lengths, target_depth, target_fraction, masks) for r in rows}


def send_webhook(url, text, log):
    try:
        req = urllib.request.Request(
            url, data=json.dumps({"text": text}).encode(), headers={"Content-Type": "application/json"}
        )
        urllib.request.urlopen(req, timeout=10).close()
    except Exception as exc:  # a failed notification must never stop monitoring
        log.warning(f"Webhook notification failed: {exc}")


class StopSignal:
    """Ctrl+C, SIGTERM or Enter (interactive terminals only) all set the same flag."""

    def __init__(self):
        self.event = threading.Event()

    def _handler(self, signum, frame):
        if self.event.is_set():
            raise KeyboardInterrupt  # second Ctrl+C: stop waiting for a clean shutdown
        self.event.set()

    def install(self):
        signal.signal(signal.SIGINT, self._handler)
        signal.signal(signal.SIGTERM, self._handler)
        if sys.stdin is not None and sys.stdin.isatty():

            def wait_for_enter():
                try:
                    if sys.stdin.readline():  # '' means EOF: ignore it rather than stopping
                        self.event.set()
                except Exception:
                    pass

            threading.Thread(target=wait_for_enter, daemon=True).start()
            return True
        return False


def monitor(console, rows, refinfo, fastq_pass, outdir, run_name, host="127.0.0.1", port=6010, target_depth=20,
            target_fraction=0.9, denominator="genome", poll_interval=5, webhook=None, open_browser=True,
            clear_annotations=False, min_identity=0.8, threads=4):  # fmt: skip
    with socket.socket() as s:
        if s.connect_ex((host if host != "0.0.0.0" else "127.0.0.1", port)) == 0:
            raise MonitorError(f"Port {port} is already in use (another nanosentinel running?). Use --port.")
    trackers = build_trackers(rows, refinfo, target_depth, target_fraction, denominator)
    annotations = outdir / "annotations"
    if clear_annotations:
        shutil.rmtree(annotations, ignore_errors=True)
    reader = AnnotationReader(annotations, trackers, refinfo.panel_hash)
    samples = [r for r in rows if not r.neg_control]
    controls = [r for r in rows if r.neg_control]
    names = {r.barcode: r.sample for r in rows}

    def snapshot_data(final=False):
        return dashboard.build_data(run_name, rows, trackers, refinfo, reader.run_start, target_depth, target_fraction,
                                    denominator, annotator.errors, final)  # fmt: skip

    def write_reports():
        report.write_all(outdir, rows, trackers, reader.run_start, target_depth, target_fraction, denominator)
        reader.changed = False

    annotator = Annotator(fastq_pass, [r.barcode for r in rows], annotations, refinfo, min_identity, threads, console.log)
    try:
        server = dashboard.DashboardServer(host, port)
    except OSError as exc:
        raise MonitorError(f"Cannot serve the dashboard on {host}:{port}: {exc}") from None
    stop = StopSignal()
    try:
        server.start()
        annotator.start()
        interactive = stop.install()
        reader.poll()
        server.update(snapshot_data())
        url = f"http://{'localhost' if host in ('127.0.0.1', '0.0.0.0') else host}:{port}"
        if open_browser:
            try:
                webbrowser.open(url, new=1)
            except Exception:
                pass
        console.say(f"Dashboard: {url}")
        console.say(f"Monitoring {len(samples)} sample(s) and {len(controls)} control(s) "
                    f"for >={target_depth}x over >={target_fraction:.0%} of the {denominator}")  # fmt: skip
        console.say("Press Enter or Ctrl+C to stop." if interactive else "Press Ctrl+C to stop.")

        all_notified = warned = False
        last_poll = last_report = last_log = last_data = 0.0
        write_reports()
        while not stop.event.is_set():
            now = time.time()
            if now - last_poll >= poll_interval:
                last_poll = now
                seen = len(reader.done)
                newly = reader.poll()
                fresh = len(reader.done) != seen  # new annotation files were read in this poll
                if newly:
                    _announce(console, newly, names, trackers, reader.run_start, {r.barcode for r in controls})
                    write_reports()
                    last_report = now
                if fresh or now - last_data >= 30:
                    server.update(snapshot_data())
                    last_data = now
                if annotator.errors and not warned:
                    warned = True
                    console.warn(annotator.errors[0] + (f" (+{len(annotator.errors) - 1} more, see the log)" if len(annotator.errors) > 1 else ""))
                if not all_notified and samples and all(trackers[r.barcode].reached for r in samples):
                    all_notified = True
                    msg = f"All {len(samples)} non-control sample(s) have reached the target depth. Monitoring continues."
                    console.say("*** " + msg, bell=True)
                    if webhook:
                        send_webhook(webhook, f"nanosentinel ({run_name or 'run'}): {msg}", console.log)
            if reader.changed and now - last_report >= 60:
                write_reports()
                last_report = now
            line = _status_line(trackers, samples, reader.run_start, interactive, annotator.files_annotated)
            console.status(line)
            if not console.tty and now - last_log >= 300:
                last_log = now
                console.say(line)
            stop.event.wait(0.5)
    finally:
        # Ignore further interrupts only while shutting down, so a second Ctrl+C cannot leave anything half-stopped.
        old_int = signal.signal(signal.SIGINT, signal.SIG_IGN)
        old_term = signal.signal(signal.SIGTERM, signal.SIG_IGN)
        console.end_status()
        console.say("Stopping ...")
        annotator.stop()
        if annotator.is_alive():
            annotator.join(15)
        server.stop()
        signal.signal(signal.SIGINT, old_int if callable(old_int) else signal.default_int_handler)
        signal.signal(signal.SIGTERM, old_term if callable(old_term) else signal.SIG_DFL)
        reader.poll()
        write_reports()
        dashboard.write_snapshot(outdir / "dashboard.html", snapshot_data(final=True))
    _summary(console, samples, controls, trackers, outdir)


def _announce(console, newly, names, trackers, run_start, control_barcodes):
    """Say which samples just reached target (one line each, or a single line when many arrive at once)."""
    for b in [b for b in newly if b in control_barcodes]:
        console.say(f"  WARNING: negative control {names[b]} reached the target depth - check for contamination", bell=True)
    newly = [b for b in newly if b not in control_barcodes]
    if len(newly) > 5:
        console.say(f"{len(newly)} samples reached target depth (see time_to_target.tsv)", bell=True)
        return
    for b in newly:
        t = trackers[b]
        when = hms((t.reached_at - run_start).total_seconds()) if t.reached_at and run_start else "?"
        console.say(f"  {names[b]} reached target at run time {when}", bell=True)


def _status_line(trackers, samples, run_start, interactive, files):
    reads = sum(t.reads for t in trackers.values())
    done = sum(trackers[r.barcode].reached for r in samples)
    latest = max((t.latest for t in trackers.values() if t.latest), default=None)
    run_time = hms((latest - run_start).total_seconds()) if latest and run_start else "--:--:--"
    hint = "Enter/Ctrl+C to stop" if interactive else "Ctrl+C to stop"
    return f"run time {run_time} | reads {reads:,} | at target {done}/{len(samples)} | {hint}"


def _summary(console, samples, controls, trackers, outdir):
    done = [r for r in samples if trackers[r.barcode].reached]
    console.say(f"Finished: {len(done)}/{len(samples)} sample(s) reached the target depth. Reports and dashboard.html in {outdir}")
    for r in samples:
        t = trackers[r.barcode]
        if not t.reached:
            ref, frac = t.best()
            console.say(f"  not reached: {r.sample} ({t.reads:,} reads, best {frac:.0%} on {ref or 'no reference'})")
    counts = [trackers[r.barcode].reads for r in samples]
    median = statistics.median(counts) if counts else 0
    for r in controls:
        n = trackers[r.barcode].reads
        pct = f" ({n / median:.1%} of the median sample)" if median else ""
        console.say(f"  control {r.sample}: {n:,} reads{pct}")
