import numpy as np
import pytest

from cemct.phase_fraction import calculate_phase_fractions


def test_two_phase_fraction():
    volume = np.array(
        [
            [[1, 1], [1, 1]],
            [[2, 2], [2, 2]],
        ],
        dtype=np.uint8,
    )

    phases = {
        1: "Phase A",
        2: "Phase B",
    }

    results = calculate_phase_fractions(
        volume=volume,
        phases=phases,
        voxel_size_um=(1.0, 2.0, 3.0),
    )

    assert len(results) == 2

    phase_a = results.loc[results["label"] == 1].iloc[0]
    phase_b = results.loc[results["label"] == 2].iloc[0]

    assert phase_a["voxel_count"] == 4
    assert phase_b["voxel_count"] == 4

    assert phase_a["volume_fraction_percent"] == pytest.approx(50.0)
    assert phase_b["volume_fraction_percent"] == pytest.approx(50.0)

    # Each voxel has a physical volume of 1 x 2 x 3 = 6 um^3.
    assert phase_a["physical_volume_um3"] == pytest.approx(24.0)
    assert phase_b["physical_volume_um3"] == pytest.approx(24.0)