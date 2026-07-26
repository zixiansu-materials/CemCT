"""
Command-line interface for CemCT.

The command-line interface provides user-facing access to CemCT without
placing input(), print(), or argument-parsing logic inside the scientific
calculation modules.

Initial commands
----------------
cemct --help
cemct --version
cemct inspect INPUT_FILE
"""

from __future__ import annotations

import argparse
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import numpy as np

from cemct.io import load_labelled_tiff


def get_cemct_version() -> str:
    """
    Return the installed CemCT package version.
    """
    try:
        return version("cemct")
    except PackageNotFoundError:
        return "0+unknown"


def inspect_tiff(file_path: str | Path) -> dict[str, object]:
    """
    Inspect an integer-labelled three-dimensional TIFF volume.

    Parameters
    ----------
    file_path
        Path to the labelled 3D TIFF volume.

    Returns
    -------
    dict
        Input path, shape, data type, phase labels, voxel counts, and phase
        percentages.
    """
    path = Path(file_path)
    volume = load_labelled_tiff(path)

    labels, voxel_counts = np.unique(
        volume,
        return_counts=True,
    )

    total_voxels = int(volume.size)

    phase_information = []

    for label, voxel_count in zip(labels, voxel_counts):
        count = int(voxel_count)

        phase_information.append(
            {
                "label": int(label),
                "voxel_count": count,
                "volume_fraction_percent": (
                    100.0 * count / total_voxels
                ),
            }
        )

    return {
        "file_path": path.resolve(),
        "shape_zyx": tuple(int(value) for value in volume.shape),
        "dtype": str(volume.dtype),
        "total_voxels": total_voxels,
        "phases": phase_information,
    }


def print_inspection_result(
    inspection_result: dict[str, object],
) -> None:
    """
    Print a TIFF inspection result to the command line.
    """
    print("CemCT labelled-volume inspection")
    print("=" * 40)
    print(f"Input file: {inspection_result['file_path']}")
    print(f"Shape (Z, Y, X): {inspection_result['shape_zyx']}")
    print(f"Data type: {inspection_result['dtype']}")
    print(f"Total voxels: {inspection_result['total_voxels']}")
    print()
    print("Observed phase labels")
    print("-" * 40)
    print(
        f"{'Label':>10}"
        f"{'Voxel count':>15}"
        f"{'Fraction (%)':>15}"
    )

    phases = inspection_result["phases"]

    for phase in phases:
        print(
            f"{phase['label']:>10d}"
            f"{phase['voxel_count']:>15d}"
            f"{phase['volume_fraction_percent']:>15.6f}"
        )


def run_inspect_command(arguments: argparse.Namespace) -> int:
    """
    Execute the ``cemct inspect`` command.
    """
    result = inspect_tiff(arguments.input_file)
    print_inspection_result(result)

    return 0


def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the CemCT argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="cemct",
        description=(
            "Post-segmentation XRmicroCT data analysis "
            "for cementitious materials."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {get_cemct_version()}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="commands",
    )

    inspect_parser = subparsers.add_parser(
        "inspect",
        help="Inspect labels and dimensions in a labelled 3D TIFF.",
        description=(
            "Display the shape, data type, observed labels, voxel counts, "
            "and phase percentages of a labelled 3D TIFF volume."
        ),
    )

    inspect_parser.add_argument(
        "input_file",
        type=Path,
        help="Path to the labelled 3D TIFF volume.",
    )

    inspect_parser.set_defaults(
        command_function=run_inspect_command,
    )

    return parser


def main(
    argv: Sequence[str] | None = None,
) -> int:
    """
    Run the CemCT command-line interface.
    """
    parser = build_parser()
    arguments = parser.parse_args(argv)

    if not hasattr(arguments, "command_function"):
        parser.print_help()
        return 0

    try:
        return int(arguments.command_function(arguments))
    except (FileNotFoundError, TypeError, ValueError, OSError) as error:
        print(
            f"CemCT error: {error}",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())