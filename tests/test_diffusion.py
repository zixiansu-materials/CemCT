"""Tests for finite-difference diffusion analysis."""

from types import SimpleNamespace

import numpy as np
import pytest

import cemct.diffusion as diffusion_module
from cemct.diffusion import (
    calculate_diffusion_tortuosity,
    calculate_directional_diffusion,
)


def fake_tortuosity_fd(im, axis, solver=None):
    """Return controlled finite-difference results."""

    porosity = float(np.count_nonzero(im) / im.size)

    concentration = np.zeros(im.shape, dtype=np.float32)

    coordinate_values = np.linspace(
        1.0,
        0.0,
        im.shape[axis],
        dtype=np.float32,
    )

    reshape_dimensions = [1, 1, 1]
    reshape_dimensions[axis] = im.shape[axis]

    concentration_values = coordinate_values.reshape(
        reshape_dimensions
    )

    concentration[:] = np.broadcast_to(
        concentration_values,
        im.shape,
    )

    concentration[~im] = np.nan

    return SimpleNamespace(
        tortuosity=2.0,
        formation_factor=4.0,
        original_porosity=porosity,
        effective_porosity=0.5,
        im_conc=concentration,
        time=0.25,
    )


def test_single_direction_transport(monkeypatch):
    """Formation factor should be converted to relative diffusivity."""

    monkeypatch.setattr(
        diffusion_module.ps.simulations,
        "tortuosity_fd",
        fake_tortuosity_fd,
    )

    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[2, 2, :] = True

    result = calculate_diffusion_tortuosity(
        connected_mask=mask,
        direction="X",
        solver=object(),
    )

    assert result["status"] == "calculated"
    assert result["axis"] == 2
    assert result["connected_voxel_count"] == 5

    assert result["formation_factor"] == pytest.approx(4.0)
    assert result["relative_effective_diffusivity"] == pytest.approx(0.25)
    assert result["diffusion_tortuosity"] == pytest.approx(2.0)
    assert result["effective_porosity"] == pytest.approx(0.5)

    assert result["concentration_field"].shape == mask.shape


def test_directional_pipeline_skips_nonpercolating_directions(
    monkeypatch,
):
    """Only the X direction should be simulated for an X-spanning line."""

    monkeypatch.setattr(
        diffusion_module.ps.simulations,
        "tortuosity_fd",
        fake_tortuosity_fd,
    )

    volume = np.zeros((5, 5, 5), dtype=np.uint8)
    volume[2, 2, :] = 3

    result = calculate_directional_diffusion(
        volume=volume,
        phase_label=3,
        directions=("X", "Y", "Z"),
        solver=object(),
    )

    assert result["direction_results"]["X"]["status"] == "calculated"

    assert (
        result["direction_results"]["Y"]["status"]
        == "not_percolating"
    )

    assert (
        result["direction_results"]["Z"]["status"]
        == "not_percolating"
    )

    assert "X" in result["concentration_fields"]
    assert "Y" not in result["concentration_fields"]
    assert "Z" not in result["concentration_fields"]


def test_direction_axis_mapping(monkeypatch):
    """X, Y, and Z must map to NumPy axes 2, 1, and 0."""

    monkeypatch.setattr(
        diffusion_module.ps.simulations,
        "tortuosity_fd",
        fake_tortuosity_fd,
    )

    masks = {
        "X": np.zeros((4, 4, 4), dtype=bool),
        "Y": np.zeros((4, 4, 4), dtype=bool),
        "Z": np.zeros((4, 4, 4), dtype=bool),
    }

    masks["X"][2, 2, :] = True
    masks["Y"][2, :, 2] = True
    masks["Z"][:, 2, 2] = True

    expected_axes = {
        "X": 2,
        "Y": 1,
        "Z": 0,
    }

    for direction, mask in masks.items():
        result = calculate_diffusion_tortuosity(
            connected_mask=mask,
            direction=direction,
            solver=object(),
        )

        assert result["axis"] == expected_axes[direction]


def test_nonpercolating_mask_raises_error():
    """A mask not touching both opposing faces must be rejected."""

    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[2, 2, 1:4] = True

    with pytest.raises(ValueError, match="does not touch both"):
        calculate_diffusion_tortuosity(
            connected_mask=mask,
            direction="X",
            solver=object(),
        )


def test_invalid_direction_raises_error():
    """Only X, Y, and Z are valid physical directions."""

    mask = np.ones((3, 3, 3), dtype=bool)

    with pytest.raises(ValueError, match="X.*Y.*Z"):
        calculate_diffusion_tortuosity(
            connected_mask=mask,
            direction="A",
            solver=object(),
        )