"""Make Greek letters printable on a Windows console."""

from __future__ import annotations

import sys


def configure_stdout() -> None:
    """Encode stdout as UTF-8.

    The default Windows console encoding is cp1252, and printing τ or π then
    raises UnicodeEncodeError before a report can be written.
    """
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, OSError, ValueError):
        return
