"""Terminal output: a few real messages plus one in-place status line, in colour when the terminal supports it.

Everything is also written, uncoloured, to a log file. Colour follows the usual conventions: it is off when the stream
is not a terminal, when NO_COLOR is set or TERM=dumb, and can be forced on with FORCE_COLOR.
"""

import logging
import os
import shutil
import sys

_CODES = {
    "bold": "1",
    "dim": "2",
    "red": "31",
    "green": "32",
    "yellow": "33",
    "blue": "34",
    "magenta": "35",
    "cyan": "36",
}

# Semantic styles used by callers, so message sites say what a line means rather than which colour it is.
STYLES = {
    "title": ("bold", "cyan"),
    "ok": ("green",),
    "good": ("bold", "green"),
    "warn": ("yellow",),
    "alert": ("bold", "red"),
    "muted": ("dim",),
    "url": ("bold", "blue"),
    "label": ("dim",),
}


def color_enabled(stream):
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    return stream.isatty() and os.environ.get("TERM", "") != "dumb"


def paint(text, style, enabled=True):
    """Wrap `text` in the ANSI codes for a semantic style name (or a tuple of colour names)."""
    if not enabled or not style:
        return text
    names = STYLES.get(style, (style,)) if isinstance(style, str) else style
    codes = ";".join(_CODES[n] for n in names if n in _CODES)
    return f"\033[{codes}m{text}\033[0m" if codes else text


def error(message):
    """Print a fatal error to stderr (used by the CLI before exiting)."""
    on = color_enabled(sys.stderr)
    print(f"{paint('Error:', 'alert', on)} {message}", file=sys.stderr)


def notice(message, style="warn"):
    on = color_enabled(sys.stderr)
    print(paint(message, style, on), file=sys.stderr)


class Console:
    def __init__(self, logfile=None):
        self.tty = sys.stdout.isatty()
        self.color_out = color_enabled(sys.stdout)
        self.color_err = color_enabled(sys.stderr)
        self._status_shown = False
        self.log = logging.getLogger("nanosentinel")
        self.log.setLevel(logging.INFO)
        self.log.handlers.clear()
        self.log.propagate = False
        if logfile:
            handler = logging.FileHandler(logfile)
            handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            self.log.addHandler(handler)
        else:
            self.log.addHandler(logging.NullHandler())  # keep the root logger's last-resort handler from echoing

    def _clear(self):
        if self._status_shown:
            sys.stdout.write("\r\033[K")
            self._status_shown = False

    def say(self, message, bell=False, style=None):
        self._clear()
        print(("\a" if bell and self.tty else "") + paint(message, style, self.color_out), flush=True)
        self.log.info(message)

    def kv(self, label, value, style=None):
        """A "Label: value" line with a muted label (and optionally a styled value)."""
        self._clear()
        print(f"{paint(label + ':', 'label', self.color_out)} {paint(value, style, self.color_out)}", flush=True)
        self.log.info(f"{label}: {value}")

    def warn(self, message):
        self._clear()
        print(f"{paint('Warning:', 'warn', self.color_err)} {message}", file=sys.stderr, flush=True)
        self.log.warning(message)

    @staticmethod
    def plain(parts):
        """Plain-text form of status parts, for the log and for non-terminal output."""
        return " | ".join(f"{label} {value}".strip() for label, value, _ in parts)

    def status(self, parts):
        """Redraw a single status line (TTY only; elsewhere it is only logged by the caller).

        `parts` is a list of (label, value, style) or a plain string. Parts that do not fit the terminal width are
        dropped from the right rather than cut mid-colour-code.
        """
        if not self.tty:
            return
        width = shutil.get_terminal_size((100, 20)).columns - 1
        if isinstance(parts, str):
            line = parts[:width]
        else:
            sep = paint(" | ", "muted", self.color_out)
            shown, used = [], 0
            for label, value, style in parts:
                cost = len(f"{label} {value}".strip()) + (3 if shown else 0)
                if shown and used + cost > width:
                    break
                shown.append(f"{paint(label, 'label', self.color_out)} {paint(value, style, self.color_out)}".strip())
                used += cost
            line = sep.join(shown)
        sys.stdout.write("\r\033[K" + line)
        sys.stdout.flush()
        self._status_shown = True

    def end_status(self):
        self._clear()
