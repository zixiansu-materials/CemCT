"""Tests for the connectivity-analysis module."""

import numpy as np
import pytest

from cemct.connectivity import analyse_connectivity


def test_x_direction_connectivity():
    """
    A five-voxel line spanning X should be X-connected.

    One additional internal voxel is isolated, giving:

        X connectivity = 5 / 6
        Isolated fraction = 1 / 6
    """

    volume = np.zeros((5, 5, 5), dtype=np.uint8)

    # Face-connected line from X-min to X-max.
    volume[2, 2, :] = 3

    # Internal isolated pore voxel.
    volume[1, 1, 1] = 3

    result = analyse_connectivity(
        volume=volume,
        phase_label=3,
    )

    summary = result["summary"]
    masks = result["masks"]

    assert summary["selected_voxel_count"] == 6
    assert summary["connected_component_count"] == 2

    assert summary["X_connected_voxel_count"] == 5
    assert summary["Y_connected_voxel_count"] == 0
    assert summary["Z_connected_voxel_count"] == 0

    assert summary["X_connectivity"] == pytest.approx(5 / 6)
    assert summary["Y_connectivity"] == pytest.approx(0.0)
    assert summary["Z_connectivity"] == pytest.approx(0.0)

    assert summary["isolated_voxel_count"] == 1
    assert summary["isolated_fraction"] == pytest.approx(1 / 6)

    assert np.count_nonzero(masks["through_connected_X"]) == 5
    assert np.count_nonzero(masks["isolated"]) == 1


def test_edge_contact_is_not_face_connected():
    """
    Voxels sharing only an edge must remain separate components.
    """

    volume = np.zeros((3, 3, 3), dtype=np.uint8)

    volume[0, 0, 0] = 3
    volume[1, 1, 0] = 3

    result = analyse_connectivity(
        volume=volume,
        phase_label=3,
    )

    assert result["summary"]["connected_component_count"] == 2


def test_corner_contact_is_not_face_connected():
    """
    Voxels sharing only a corner must remain separate components.
    """

    volume = np.zeros((3, 3, 3), dtype=np.uint8)

    volume[0, 0, 0] = 3
    volume[1, 1, 1] = 3

    result = analyse_connectivity(
        volume=volume,
        phase_label=3,
    )

    assert result["summary"]["connected_component_count"] == 2


def test_missing_label_raises_error():
    """An absent analysis label should produce a clear error."""

    volume = np.zeros((3, 3, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="absent"):
        analyse_connectivity(
            volume=volume,
            phase_label=3,
        )