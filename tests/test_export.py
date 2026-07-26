"""Tests for the unified CemCT export layer."""

import json

import numpy as np
import pandas as pd
import tifffile

from cemct.export import (
    export_analysis_bundle,
    export_array_archive,
    export_binary_mask,
    export_excel_workbook,
    export_metadata_json,
    export_scalar_volume,
)


def test_export_binary_mask(tmp_path):
    mask = np.zeros((4, 5, 6), dtype=bool)
    mask[1:3, 2:4, 3:5] = True

    output_path = tmp_path / "binary_mask.tif"

    returned_path = export_binary_mask(
        mask=mask,
        output_path=output_path,
    )

    saved_volume = tifffile.imread(returned_path)

    assert returned_path.exists()
    assert saved_volume.dtype == np.uint8
    assert set(np.unique(saved_volume)) == {0, 255}


def test_export_scalar_volume(tmp_path):
    volume = np.arange(
        4 * 5 * 6,
        dtype=np.float64,
    ).reshape((4, 5, 6))

    output_path = tmp_path / "scalar_volume.tif"

    returned_path = export_scalar_volume(
        volume=volume,
        output_path=output_path,
    )

    saved_volume = tifffile.imread(returned_path)

    assert returned_path.exists()
    assert saved_volume.dtype == np.float32
    np.testing.assert_allclose(
        saved_volume,
        volume.astype(np.float32),
    )


def test_export_excel_workbook(tmp_path):
    summary = pd.DataFrame(
        {
            "Direction": ["X", "Y", "Z"],
            "Permeability (m^2)": [1e-18, 2e-18, 3e-18],
        }
    )

    output_path = tmp_path / "results.xlsx"

    returned_path = export_excel_workbook(
        tables={"Permeability Summary": summary},
        output_path=output_path,
        metadata={
            "phase_label": 3,
            "voxel_size_um": 0.7,
        },
    )

    workbook = pd.ExcelFile(returned_path)

    assert returned_path.exists()
    assert "Permeability Summary" in workbook.sheet_names
    assert "Metadata" in workbook.sheet_names


def test_export_metadata_json(tmp_path):
    metadata = {
        "phase_label": np.int64(3),
        "voxel_size_um": np.float64(0.7),
        "directions": ("X", "Y", "Z"),
    }

    output_path = tmp_path / "metadata.json"

    returned_path = export_metadata_json(
        metadata=metadata,
        output_path=output_path,
    )

    with returned_path.open(
        mode="r",
        encoding="utf-8",
    ) as file:
        loaded_metadata = json.load(file)

    assert loaded_metadata["phase_label"] == 3
    assert loaded_metadata["voxel_size_um"] == 0.7
    assert loaded_metadata["directions"] == ["X", "Y", "Z"]


def test_export_array_archive(tmp_path):
    arrays = {
        "X pressure": np.array([1.0, 0.5, 0.0]),
        "Y pressure": np.array([1.0, 0.7, 0.0]),
    }

    output_path = tmp_path / "arrays.npz"

    returned_path = export_array_archive(
        arrays=arrays,
        output_path=output_path,
    )

    loaded_archive = np.load(returned_path)

    assert returned_path.exists()
    assert "X_pressure" in loaded_archive.files
    assert "Y_pressure" in loaded_archive.files


def test_export_complete_analysis_bundle(tmp_path):
    input_path = tmp_path / "sample.tif"
    input_path.touch()

    summary = pd.DataFrame(
        {
            "Direction": ["X"],
            "Intrinsic permeability (m^2)": [1.0e-18],
        }
    )

    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[:, 2, 2] = True

    scalar_volume = np.ones((5, 5, 5), dtype=float)

    created_paths = export_analysis_bundle(
        input_file=input_path,
        analysis_name="permeability",
        tables={"Summary": summary},
        masks={"through_connected_X": mask},
        scalar_volumes={
            "local_diameter_um": scalar_volume,
        },
        arrays={
            "X_pressure": np.array([1.0, 0.5, 0.0]),
        },
        metadata={
            "phase_label": 3,
            "voxel_size_um": 0.7,
        },
        output_directory=tmp_path / "outputs",
    )

    assert created_paths["excel"].exists()
    assert created_paths["metadata_json"].exists()
    assert created_paths["array_archive"].exists()

    assert created_paths["binary_masks"][
        "through_connected_X"
    ].exists()

    assert created_paths["scalar_volumes"][
        "local_diameter_um"
    ].exists()