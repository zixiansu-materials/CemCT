"""
Synthetic structural validation for CemCT.

This script validates:

1. Multi-phase fraction calculation using a controlled two-phase volume.
2. Directional connectivity using a known X-through-connected channel.
3. Separation between a through-connected channel and an isolated pore region.

Array convention:

    axis 0 = Z
    axis 1 = Y
    axis 2 = X
"""

from pathlib import Path

import numpy as np
import pandas as pd

from cemct.connectivity import analyse_connectivity
from cemct.phase_fraction import calculate_phase_fractions


OUTPUT_DIRECTORY = Path(__file__).parent / "results"
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)


def make_validation_record(
    test_name,
    expected,
    calculated,
    tolerance=1.0e-10,
):
    """Create one numerical validation record."""

    absolute_error = abs(float(calculated) - float(expected))
    passed = absolute_error <= tolerance

    return {
        "test_name": test_name,
        "expected": float(expected),
        "calculated": float(calculated),
        "absolute_error": absolute_error,
        "tolerance": float(tolerance),
        "status": "PASS" if passed else "FAIL",
    }


def validate_phase_fractions():
    """
    Validate phase fractions using an exactly divided two-phase volume.

    The volume contains:

        Label 1 = 500 voxels
        Label 2 = 500 voxels

    Therefore, both phase fractions must equal 0.5.
    """

    volume = np.ones(
        (10, 10, 10),
        dtype=np.uint8,
    )

    volume[5:, :, :] = 2

    results = calculate_phase_fractions(
        volume=volume,
        phases={
            1: "Phase 1",
            2: "Phase 2",
        },
        voxel_size_um=1.0,
    )

    results_by_label = results.set_index("label")

    records = [
        make_validation_record(
            test_name="Phase 1 voxel count",
            expected=500,
            calculated=results_by_label.loc[1, "voxel_count"],
        ),
        make_validation_record(
            test_name="Phase 2 voxel count",
            expected=500,
            calculated=results_by_label.loc[2, "voxel_count"],
        ),
        make_validation_record(
            test_name="Phase 1 volume fraction",
            expected=0.5,
            calculated=results_by_label.loc[1, "volume_fraction"],
        ),
        make_validation_record(
            test_name="Phase 2 volume fraction",
            expected=0.5,
            calculated=results_by_label.loc[2, "volume_fraction"],
        ),
    ]

    return records


def validate_directional_connectivity():
    """
    Validate directional connectivity using an ideal synthetic geometry.

    Geometry:

        X-through-connected channel = 20 voxels
        Isolated pore cube           = 8 voxels
        Total selected pore voxels   = 28 voxels

    Expected X connectivity:

        20 / 28 = 0.7142857142857143

    The isolated cube does not touch the through-connected channel.
    """

    volume = np.ones(
        (20, 20, 20),
        dtype=np.uint8,
    )

    pore_label = 3

    # X-through-connected pore channel.
    volume[10, 10, :] = pore_label

    # Isolated 2 x 2 x 2 pore region.
    volume[2:4, 2:4, 2:4] = pore_label

    results = analyse_connectivity(
        volume=volume,
        phase_label=pore_label,
    )

    summary = results["summary"]

    expected_selected_voxels = 28
    expected_x_connected_voxels = 20
    expected_x_connectivity = 20 / 28

    records = [
        make_validation_record(
            test_name="Selected pore voxel count",
            expected=expected_selected_voxels,
            calculated=summary["selected_voxel_count"],
        ),
        make_validation_record(
            test_name="X-connected voxel count",
            expected=expected_x_connected_voxels,
            calculated=summary["X_connected_voxel_count"],
        ),
        make_validation_record(
            test_name="X-directional connectivity",
            expected=expected_x_connectivity,
            calculated=summary["X_connectivity"],
        ),
        make_validation_record(
            test_name="Y-connected voxel count",
            expected=0,
            calculated=summary["Y_connected_voxel_count"],
        ),
        make_validation_record(
            test_name="Z-connected voxel count",
            expected=0,
            calculated=summary["Z_connected_voxel_count"],
        ),
        make_validation_record(
            test_name="Isolated voxel count",
            expected=8,
            calculated=summary["isolated_voxel_count"],
        ),
    ]

    return records


def main():
    """Run all synthetic structural validation cases."""

    validation_records = []

    validation_records.extend(
        validate_phase_fractions()
    )

    validation_records.extend(
        validate_directional_connectivity()
    )

    validation_table = pd.DataFrame(
        validation_records
    )

    output_path = (
        OUTPUT_DIRECTORY
        / "synthetic_structural_validation.csv"
    )

    validation_table.to_csv(
        output_path,
        index=False,
    )

    print("\nCemCT Synthetic Structural Validation")
    print("=" * 72)
    print(validation_table.to_string(index=False))
    print("=" * 72)

    all_passed = (
        validation_table["status"]
        .eq("PASS")
        .all()
    )

    print(f"\nAll validation cases passed: {all_passed}")
    print(f"Results saved to: {output_path.resolve()}")

    if not all_passed:
        raise RuntimeError(
            "One or more synthetic structural validation cases failed."
        )


if __name__ == "__main__":
    main()