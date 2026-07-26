"""
Module entry point for CemCT.

This file enables CemCT to be executed using:

    python -m cemct
"""

from cemct.cli import main


if __name__ == "__main__":
    raise SystemExit(main())