"""The live monitoring loop: RAMPART for display, our own tracking for time-to-target."""

import json
import signal
import statistics
import sys
import threading
import time
import urllib.request

import numpy as np

from . import rampart, report
from .coverage import AnnotationReader, SampleTracker
from .errors import MonitorError
from .protocol import amplicon_mask_intervals
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


def monitor(console, rows, refinfo, protocol_dir, fastq_pass, outdir, run_name, ports, target_depth, target_fraction,
            denominator="genome", poll_interval=10, webhook=None, open_browser=True, clear_annotations=False):  # fmt: skip
    rampart.check_tools()
    rampart.check_path(fastq_pass)
    rampart.check_ports(ports)
    trackers = build_trackers(rows, refinfo, target_depth, target_fraction, denominator)
    annotations = outdir / "annotations"  # RAMPART writes here because its working directory is outdir
    reader = AnnotationReader(annotations, trackers)
    rampart_log = outdir / "rampart.log"
    samples = [r for r in rows if not r.neg_control]
    controls = [r for r in rows if r.neg_control]
    names = {r.barcode: r.sample for r in rows}

    def write_reports():
        report.write_all(outdir, rows, trackers, reader.run_start, target_depth, target_fraction, denominator)
        reader.changed = False

    proc = rampart.launch(
        protocol_dir, fastq_pass, [f"{r.barcode}={r.sample}" for r in rows], run_name, ports,
        outdir, rampart_log, clear=clear_annotations,
    )  # fmt: skip
    stop = StopSignal()
    try:
        if not rampart.wait_until_listening(proc, ports[0]):
            raise MonitorError(f"RAMPART did not start. Last log lines ({rampart_log}):\n{rampart.tail(rampart_log)}")
        interactive = stop.install()
        if open_browser:
            rampart.open_browser(ports[0])
        console.say(f"RAMPART: http://localhost:{ports[0]}")
        console.say(f"Monitoring {len(samples)} sample(s) and {len(controls)} control(s) "
                    f"for >={target_depth}x over >={target_fraction:.0%} of the {denominator}")  # fmt: skip
        console.say("Press Enter or Ctrl+C to stop." if interactive else "Press Ctrl+C to stop.")

        all_notified = False
        last_poll = last_log = 0.0
        write_reports()
        last_report = time.time()
        while not stop.event.is_set():
            now = time.time()
            if proc.poll() is not None:
                raise MonitorError(f"RAMPART exited unexpectedly. Last log lines ({rampart_log}):\n{rampart.tail(rampart_log)}")
            if now - last_poll >= poll_interval:
                last_poll = now
                newly = reader.poll()
                if newly:
                    _announce(console, newly, names, trackers, reader.run_start, {r.barcode for r in controls})
                    write_reports()
                    last_report = now
                if not all_notified and samples and all(trackers[r.barcode].reached for r in samples):
                    all_notified = True
                    msg = f"All {len(samples)} non-control sample(s) have reached the target depth. Monitoring continues."
                    console.say("*** " + msg, bell=True)
                    if webhook:
                        send_webhook(webhook, f"nanosentinel ({run_name or 'run'}): {msg}", console.log)
            if reader.changed and now - last_report >= 60:
                write_reports()
                last_report = now
            line = _status_line(trackers, samples, reader.run_start, interactive)
            console.status(line)
            if not console.tty and now - last_log >= 300:
                last_log = now
                console.say(line)
            stop.event.wait(0.5)
    finally:
        # Ignore further interrupts only while RAMPART is being stopped, so a second Ctrl+C cannot orphan it.
        old_int = signal.signal(signal.SIGINT, signal.SIG_IGN)
        old_term = signal.signal(signal.SIGTERM, signal.SIG_IGN)
        console.end_status()
        console.say("Stopping RAMPART ...")
        rampart.stop(proc)
        signal.signal(signal.SIGINT, old_int if callable(old_int) else signal.default_int_handler)
        signal.signal(signal.SIGTERM, old_term if callable(old_term) else signal.SIG_DFL)
        reader.poll(final=True)
        write_reports()
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


def _status_line(trackers, samples, run_start, interactive):
    reads = sum(t.reads for t in trackers.values())
    done = sum(trackers[r.barcode].reached for r in samples)
    latest = max((t.latest for t in trackers.values() if t.latest), default=None)
    run_time = hms((latest - run_start).total_seconds()) if latest and run_start else "--:--:--"
    hint = "Enter/Ctrl+C to stop" if interactive else "Ctrl+C to stop"
    return f"run time {run_time} | reads {reads:,} | at target {done}/{len(samples)} | {hint}"


def _summary(console, samples, controls, trackers, outdir):
    done = [r for r in samples if trackers[r.barcode].reached]
    console.say(f"Finished: {len(done)}/{len(samples)} sample(s) reached the target depth. Reports in {outdir}")
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
