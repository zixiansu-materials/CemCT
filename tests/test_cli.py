"""Tests for the CemCT command-line interface."""

import numpy as np
import pytest
import tifffile

from cemct.cli import build_parser
from cemct.cli import inspect_tiff
from cemct.cli import main


def test_cli_without_command_prints_help(capsys):
    return_code = main([])

    captured = capsys.readouterr()

    assert return_code == 0
    assert "usage:" in captured.out
    assert "inspect" in captured.out


def test_cli_help(capsys):
    with pytest.raises(SystemExit) as exception:
        main(["--help"])

    captured = capsys.readouterr()

    assert exception.value.code == 0
    assert "Post-segmentation XRmicroCT" in captured.out


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as exception:
        main(["--version"])

    captured = capsys.readouterr()

    assert exception.value.code == 0
    assert "cemct" in captured.out


def test_inspect_labelled_tiff(tmp_path):
    volume = np.zeros((4, 5, 6), dtype=np.uint8)
    volume[:, :, :3] = 1
    volume[:, :, 3:] = 3

    input_path = tmp_path / "labelled_volume.tif"

    tifffile.imwrite(
        input_path,
        volume,
        photometric="minisblack",
    )

    result = inspect_tiff(input_path)

    assert result["shape_zyx"] == (4, 5, 6)
    assert result["dtype"] == "uint8"
    assert result["total_voxels"] == 120
    assert len(result["phases"]) == 2

    observed_labels = {
        phase["label"]
        for phase in result["phases"]
    }

    assert observed_labels == {1, 3}


def test_cli_inspect_command(tmp_path, capsys):
    volume = np.ones((3, 4, 5), dtype=np.uint8)
    volume[:, :, 2:] = 2

    input_path = tmp_path / "sample.tif"

    tifffile.imwrite(
        input_path,
        volume,
        photometric="minisblack",
    )

    return_code = main(
        [
            "inspect",
            str(input_path),
        ]
    )

    captured = capsys.readouterr()

    assert return_code == 0
    assert "Shape (Z, Y, X): (3, 4, 5)" in captured.out
    assert "Observed phase labels" in captured.out


def test_cli_missing_file_returns_error(tmp_path, capsys):
    missing_path = tmp_path / "missing.tif"

    return_code = main(
        [
            "inspect",
            str(missing_path),
        ]
    )

    captured = capsys.readouterr()

    assert return_code == 1
    assert "CemCT error:" in captured.err


def test_parser_contains_inspect_command():
    parser = build_parser()

    assert parser.prog == "cemct"