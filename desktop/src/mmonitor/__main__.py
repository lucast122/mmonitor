"""Allow running MMonitor as ``python -m mmonitor``."""
import sys

from mmonitor.cli.main import main


def _exception_handler(exc_type, exc_value, exc_tb):
    print(f"{exc_type.__name__}: {exc_value}")


sys.excepthook = _exception_handler

if __name__ == "__main__":
    main()
