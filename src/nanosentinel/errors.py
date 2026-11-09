"""Exceptions that carry a message meant to be shown to the user as-is."""


class SentinelError(Exception):
    """Base class: the CLI prints the message (no traceback) and exits 1."""


class SheetError(SentinelError):
    pass


class RunDirError(SentinelError):
    pass


class RunNotReady(RunDirError):
    """The run folder / fastq_pass does not exist (yet): worth waiting for."""


class ProtocolError(SentinelError):
    pass


class MonitorError(SentinelError):
    pass


class ConcatError(SentinelError):
    pass
