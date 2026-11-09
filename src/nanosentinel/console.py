"""Quiet terminal output: a few real messages plus one in-place status line. Everything also goes to a log file."""

import logging
import shutil
import sys


class Console:
    def __init__(self, logfile=None):
        self.tty = sys.stdout.isatty()
        self._status_shown = False
        self.log = logging.getLogger("nanosentinel")
        self.log.setLevel(logging.INFO)
        self.log.handlers.clear()
        self.log.propagate = False
        if logfile:
            handler = logging.FileHandler(logfile)
            handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            self.log.addHandler(handler)

    def _clear(self):
        if self._status_shown:
            sys.stdout.write("\r\033[K")
            self._status_shown = False

    def say(self, message, bell=False):
        self._clear()
        print(("\a" if bell and self.tty else "") + message, flush=True)
        self.log.info(message)

    def warn(self, message):
        self._clear()
        print(f"Warning: {message}", file=sys.stderr, flush=True)
        self.log.warning(message)

    def status(self, line):
        """Redraw a single status line (TTY only; elsewhere it is only logged by the caller)."""
        if not self.tty:
            return
        width = shutil.get_terminal_size((100, 20)).columns - 1
        sys.stdout.write("\r\033[K" + line[:width])
        sys.stdout.flush()
        self._status_shown = True

    def end_status(self):
        self._clear()
