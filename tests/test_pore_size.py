"""Tests for the local-thickness pore-size module."""

import numpy as np
import pytest

import cemct.pore_size as pore_size_module
from cemct.pore_size import (
    calculate_local_diameter_map,
    calculate_pore_size_distribution,
)


def fake_local_thickness(im, method="dt", sizes=25):
    """
    Return a controlled radius map for testing.

    Every selected voxel is assigned a radius of one voxel.
    """

    return np.asarray(im, dtype=np.float32)


def test_local_radius_is_converted_to_physical_diameter(monkeypatch):
    """A one-voxel radius at 0.7 um should produce a 1.4 um diameter."""

    monkeypatch.setattr(
        pore_size_module.ps.filters,
        "local_thickness",
        fake_local_thickness,
    )

    volume = np.zeros((2, 2, 2), dtype=np.uint8)
    volume[0, :, :] = 3

    selected_mask, diameter_map = calculate_local_diameter_map(
        volume=volume,
        phase_label=3,
        voxel_size_um=0.7,
    )

    assert np.count_nonzero(selected_mask) == 4
    assert np.allclose(diameter_map[selected_mask], 1.4)
    assert np.all(diameter_map[~selected_mask] == 0.0)


def test_pore_size_summary_and_classes(monkeypatch):
    """Controlled 2 um diameters should belong to the 1-10 um class."""

    monkeypatch.setattr(
        pore_size_module.ps.filters,
        "local_thickness",
        fake_local_thickness,
    )

    volume = np.zeros((2, 2, 2), dtype=np.uint8)
    volume[0, :, :] = 3

    result = calculate_pore_size_distribution(
        volume=volume,
        phase_label=3,
        voxel_size_um=1.0,
        number_of_bins=4,
    )

    summary = result["summary"]
    classes = result["size_classes"]

    assert summary["selected_voxel_count"] == 4
    assert summary["selected_phase_fraction_percent"] == pytest.approx(50.0)

    assert summary["minimum_local_diameter_um"] == pytest.approx(2.0)
    assert summary["median_local_diameter_um"] == pytest.approx(2.0)
    assert summary["mean_local_diameter_um"] == pytest.approx(2.0)
    assert summary["maximum_local_diameter_um"] == pytest.approx(2.0)

    small = classes.loc[
        classes["class_name"] == "Small capillary pores"
    ].iloc[0]

    large = classes.loc[
        classes["class_name"] == "Large capillary pores"
    ].iloc[0]

    macro = classes.loc[
        classes["class_name"] == "Macropores"
    ].iloc[0]

    assert small["voxel_count"] == 0
    assert large["voxel_count"] == 4
    assert macro["voxel_count"] == 0

    assert large["porosity_contribution_percent"] == pytest.approx(50.0)


def test_distribution_contains_all_selected_voxels(monkeypatch):
    """Histogram voxel counts should sum to the selected voxel count."""

    monkeypatch.setattr(
        pore_size_module.ps.filters,
        "local_thickness",
        fake_local_thickness,
    )

    volume = np.zeros((3, 3, 3), dtype=np.uint8)
    volume[1, :, :] = 3

    result = calculate_pore_size_distribution(
        volume=volume,
        phase_label=3,
        voxel_size_um=1.0,
        number_of_bins=5,
    )

    histogram_total = int(
        result["distribution"]["voxel_count"].sum()
    )

    assert histogram_total == 9


def test_missing_phase_label_raises_error():
    """An absent phase label must produce a clear error."""

    volume = np.zeros((3, 3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="absent"):
        calculate_local_diameter_map(
            volume=volume,
            phase_label=3,
            voxel_size_um=1.0,
        )


def test_invalid_voxel_size_raises_error():
    """Voxel size must be positive."""

    volume = np.ones((3, 3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="greater than zero"):
        calculate_local_diameter_map(
            volume=volume,
            phase_label=1,
            voxel_size_um=0.0,
        )