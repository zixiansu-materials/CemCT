"""Tests for the CemCT permeability module."""

import numpy as np
import pytest

from cemct.permeability import (
    calculate_darcy_permeability,
    get_boundary_width,
    get_domain_geometry,
    normalise_directions,
    to_physical_xyz,
    validate_labelled_volume,
)


def test_darcy_permeability_known_value():
    permeability = calculate_darcy_permeability(
        flow_rate_m3_s=2.0e-12,
        dynamic_viscosity_pa_s=1.0e-3,
        sample_length_m=1.0e-3,
        cross_sectional_area_m2=2.0e-6,
        pressure_drop_pa=2.0,
    )

    assert permeability == pytest.approx(5.0e-13)


def test_zyx_to_xyz_transposition():
    mask_zyx = np.zeros((2, 3, 4), dtype=bool)
    mask_zyx[1, 2, 3] = True

    mask_xyz = to_physical_xyz(mask_zyx)

    assert mask_xyz.shape == (4, 3, 2)
    assert mask_xyz[3, 2, 1]


def test_directional_boundary_width():
    assert get_boundary_width("X", 3) == [3, 0, 0]
    assert get_boundary_width("Y", 3) == [0, 3, 0]
    assert get_boundary_width("Z", 3) == [0, 0, 3]


def test_x_direction_domain_geometry():
    length_m, area_m2 = get_domain_geometry(
        original_shape_zyx=(30, 20, 10),
        direction="X",
        voxel_size_m=2.0e-6,
        boundary_width=3,
    )

    expected_length = (10 + 6) * 2.0e-6
    expected_area = 20 * 30 * (2.0e-6) ** 2

    assert length_m == pytest.approx(expected_length)
    assert area_m2 == pytest.approx(expected_area)


def test_missing_phase_label_raises_error():
    volume = np.zeros((4, 4, 4), dtype=np.uint8)
    volume[1:3, 1:3, 1:3] = 2

    with pytest.raises(ValueError, match="absent"):
        validate_labelled_volume(volume, phase_label=3)


def test_invalid_direction_raises_error():
    with pytest.raises(ValueError, match="Invalid directions"):
        normalise_directions(("X", "A"))


def test_negative_pressure_drop_raises_error():
    with pytest.raises(ValueError):
        calculate_darcy_permeability(
            flow_rate_m3_s=1.0,
            dynamic_viscosity_pa_s=1.0,
            sample_length_m=1.0,
            cross_sectional_area_m2=1.0,
            pressure_drop_pa=-1.0,
        )