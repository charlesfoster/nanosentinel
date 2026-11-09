"""The live dashboard: a JSON export of the trackers, a tiny local web server, and a self-contained final snapshot."""

import json
import threading
from collections import Counter
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

PAGE = Path(__file__).with_name("dashboard.html")
_HIST_WIDTHS = (10, 20, 25, 50, 100, 200, 500, 1000, 2000, 5000)


def _clean(x):
    return round(x, 4) if isinstance(x, float) else x


def _downsample(points, limit=150):
    if len(points) <= limit:
        return points
    step = len(points) / limit
    picked = [points[int(i * step)] for i in range(limit)]
    return picked[:-1] + [points[-1]]


def _length_histogram(trackers):
    counts = Counter()
    for t in trackers.values():
        counts.update(t.len_counts)
    if not counts:
        return {"width": 50, "counts": []}
    total, seen, cap = sum(counts.values()), 0, max(counts)
    for length in sorted(counts):  # trim the extreme tail so a few very long reads do not flatten the plot
        seen += counts[length]
        if seen >= 0.995 * total:
            cap = length
            break
    width = next((w for w in _HIST_WIDTHS if cap / w <= 80), _HIST_WIDTHS[-1])
    bins = [0] * (cap // width + 1)
    for length, n in counts.items():
        bins[min(length, cap) // width] += n
    return {"width": width, "counts": bins}


def build_data(run_name, rows, trackers, refinfo, run_start, target_depth, target_fraction, denominator,
               errors=(), final=False, bins=300):  # fmt: skip
    """Everything the page needs, as plain JSON-able Python."""

    def secs(when):
        return (when - run_start).total_seconds() if when and run_start else None

    samples = []
    for r in rows:
        t = trackers[r.barcode]
        best, frac = t.best()
        ls = t.length_stats()
        series = [[round(secs(when) / 60, 2), round(f, 4)] for when, _, _, _, _, f in t.series if secs(when) is not None]
        top = t.ref_reads.most_common(6)
        composition = [[ref, n] for ref, n in top]
        rest = t.mapped - sum(n for _, n in top)
        if rest > 0:
            composition.append(["Other", rest])
        samples.append({
            "name": r.sample, "barcode": r.barcode, "control": r.neg_control,
            "reads": t.reads, "mapped_pct": round(100 * t.mapped / t.reads, 1) if t.reads else None,
            "best_ref": best, "best_fraction": _clean(frac), "dominant": t.dominant(),
            "reached": t.reached, "reached_s": secs(t.reached_at) if t.reached else None,
            "mean_depth": round(t.mean_depth(best), 1) if best else None,
            "median_len": ls and round(ls["median"]), "n50": ls and ls["n50"],
            "series": _downsample(series), "composition": composition,
            "coverage": {"ref": best, "length": refinfo.lengths[best], "bins": t.binned_depth(best, bins)} if best else None,
            "amp": t.amplicon_depths(refinfo.coordinate_ref, refinfo.amplicons) if refinfo.amplicons else None,
        })  # fmt: skip
    latest = max((t.latest for t in trackers.values() if t.latest), default=None)
    return {
        "run": run_name or "", "final": final, "updated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        "elapsed_s": secs(latest), "reads": sum(t.reads for t in trackers.values()),
        "target_depth": target_depth, "target_fraction": target_fraction, "denominator": denominator,
        "coordinate_ref": refinfo.coordinate_ref, "amplicons": refinfo.amplicons,
        "genes": [{"name": k, **v} for k, v in refinfo.genes.items()],
        "samples": samples, "length_hist": _length_histogram(trackers), "warnings": list(errors),
    }  # fmt: skip


def dumps(data):
    return json.dumps(data, allow_nan=False, separators=(",", ":"))


def write_snapshot(path, data):
    """A single self-contained HTML file (data embedded) that opens straight from disk."""
    blob = dumps(data).replace("</", "<\\/")
    html = PAGE.read_text().replace("/*__DATA__*/null", blob)
    tmp = Path(f"{path}.tmp")
    tmp.write_text(html)
    tmp.replace(path)


class DashboardServer:
    """Serves the page and the latest data.json from memory (nothing on disk is exposed)."""

    def __init__(self, host, port):
        self.page = PAGE.read_bytes()
        self.data = b"{}"
        outer = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                path = self.path.split("?")[0]
                if path in ("/", "/index.html"):
                    body, kind = outer.page, "text/html; charset=utf-8"
                elif path == "/data.json":
                    body, kind = outer.data, "application/json"
                else:
                    self.send_error(404)
                    return
                self.send_response(200)
                self.send_header("Content-Type", kind)
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *args):  # keep the terminal quiet
                pass

        self.httpd = ThreadingHTTPServer((host, port), Handler)
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)

    def start(self):
        self.thread.start()

    def update(self, data):
        self.data = dumps(data).encode()

    def stop(self):
        self.httpd.shutdown()
        self.httpd.server_close()
