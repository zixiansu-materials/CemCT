"""Tests for the PyTrax random-walk analysis module."""

import numpy as np
import pytest

import cemct.random_walk as random_walk_module
from cemct.random_walk import (
    calculate_connected_random_walk,
    run_random_walk,
)


class FakeRandomWalk:
    """Controlled replacement for PyTrax RandomWalk."""

    def __init__(self, image, seed=False):
        self.im = np.asarray(image)
        self.seed = seed

    def run(
        self,
        nt=1000,
        nw=1,
        same_start=False,
        stride=1,
        num_proc=1,
    ):
        self.nt = int(nt)
        self.nw = int(nw)
        self.stride = int(stride)

        number_of_saved_points = self.nt // self.stride

        self.real_coords = np.zeros(
            (
                number_of_saved_points,
                self.nw,
                3,
            ),
            dtype=int,
        )

    def calc_msd(self):
        number_of_saved_points = self.nt // self.stride

        time_steps = (
            np.arange(number_of_saved_points, dtype=float)
            * self.stride
        )

        # Mean MSD slope = 1/5, therefore mean tortuosity = 5.
        self.msd = time_steps / 5.0

        self.axial_msd = np.zeros(
            (number_of_saved_points, 3),
            dtype=float,
        )

        # PyTrax scales axial MSD by the number of dimensions.
        # Axis 0 = Z: scaled slope 1/2 -> tau Z = 2.
        # Axis 1 = Y: scaled slope 1/3 -> tau Y = 3.
        # Axis 2 = X: scaled slope 1/4 -> tau X = 4.
        self.axial_msd[:, 0] = time_steps / (3.0 * 2.0)
        self.axial_msd[:, 1] = time_steps / (3.0 * 3.0)
        self.axial_msd[:, 2] = time_steps / (3.0 * 4.0)


def test_random_walk_tortuosity_calculation(monkeypatch):
    """Controlled MSD slopes should produce known tortuosities."""

    monkeypatch.setattr(
        random_walk_module.pt,
        "RandomWalk",
        FakeRandomWalk,
    )

    mask = np.ones((4, 4, 4), dtype=np.uint8)

    result = run_random_walk(
        mask=mask,
        number_of_steps=100,
        number_of_walkers=10,
        stride=10,
    )

    summary = result["summary"]

    assert summary["mean_tortuosity"] == pytest.approx(5.0)
    assert summary["X_tortuosity"] == pytest.approx(4.0)
    assert summary["Y_tortuosity"] == pytest.approx(3.0)
    assert summary["Z_tortuosity"] == pytest.approx(2.0)

    assert summary["mean_r_squared"] == pytest.approx(1.0)
    assert summary["X_r_squared"] == pytest.approx(1.0)
    assert summary["Y_r_squared"] == pytest.approx(1.0)
    assert summary["Z_r_squared"] == pytest.approx(1.0)


def test_255_binary_mask_is_accepted(monkeypatch):
    """A 0/255 binary mask should be converted internally."""

    monkeypatch.setattr(
        random_walk_module.pt,
        "RandomWalk",
        FakeRandomWalk,
    )

    mask = np.zeros((4, 4, 4), dtype=np.uint8)
    mask[:, 2, 2] = 255

    result = run_random_walk(
        mask=mask,
        number_of_steps=100,
        number_of_walkers=10,
        stride=10,
    )

    assert result["summary"]["accessible_voxel_count"] == 4
    assert result["transport_mask"].dtype == bool


def test_connected_label_pipeline(monkeypatch):
    """The high-level function should extract an X-connected label."""

    monkeypatch.setattr(
        random_walk_module.pt,
        "RandomWalk",
        FakeRandomWalk,
    )

    volume = np.zeros((5, 5, 5), dtype=np.uint8)
    volume[2, 2, :] = 3

    result = calculate_connected_random_walk(
        volume=volume,
        phase_label=3,
        connectivity_mode="X",
        number_of_steps=100,
        number_of_walkers=10,
        stride=10,
    )

    assert result["summary"]["phase_label"] == 3
    assert result["summary"]["connectivity_mode"] == "X"
    assert result["summary"]["connected_voxel_count"] == 5
    assert result["summary"]["X_connectivity"] == pytest.approx(1.0)


def test_empty_mask_raises_error():
    """An empty transport mask should be rejected."""

    mask = np.zeros((3, 3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="no accessible voxels"):
        run_random_walk(
            mask=mask,
            number_of_steps=100,
            number_of_walkers=10,
            stride=10,
        )


def test_invalid_stride_raises_error():
    """Number of steps must be divisible by stride."""

    mask = np.ones((3, 3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="divisible"):
        run_random_walk(
            mask=mask,
            number_of_steps=101,
            number_of_walkers=10,
            stride=10,
        )