"""
Analytical permeability validation for CemCT.

This script compares the CemCT pore-network permeability result with the
analytical Hagen-Poiseuille solution for an ideal straight cylindrical
capillary.

The validation contains two checks:

1. Darcy-law consistency

   The reported intrinsic permeability is recalculated from the CemCT flow
   rate, pressure drop, simulation length, fluid viscosity, and sample
   cross-sectional area.

2. Cylindrical-capillary analytical comparison

   The permeability estimated by the complete image-to-network workflow is
   compared with the Hagen-Poiseuille analytical solution.

Array convention:

    axis 0 = Z
    axis 1 = Y
    axis 2 = X

The cylindrical capillary is aligned with the physical X direction.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemct.permeability import analyse_permeability


OUTPUT_DIRECTORY = Path(__file__).parent / "results"
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

PORE_LABEL = 3
SOLID_LABEL = 1

VOXEL_SIZE_UM = 1.0
VOXEL_SIZE_M = VOXEL_SIZE_UM * 1.0e-6

DYNAMIC_VISCOSITY_PA_S = 1.0
INLET_PRESSURE_PA = 1.0
OUTLET_PRESSURE_PA = 0.0

# Exact Darcy-law reconstruction should agree almost exactly.
DARCY_CONSISTENCY_TOLERANCE = 1.0e-8

# The complete SNOW2/OpenPNM calculation contains voxelisation and
# pore-network approximation errors. This tolerance must be declared
# before examining the result and must not be adjusted afterwards.
CAPILLARY_RELATIVE_TOLERANCE = 0.25


def create_cylindrical_x_channel(
    shape=(60, 60, 60),
    radius_voxels=10,
):
    """
    Create a straight cylindrical capillary aligned with X.

    Parameters
    ----------
    shape
        Volume shape in (Z, Y, X) order.
    radius_voxels
        Nominal cylindrical radius in voxels.

    Returns
    -------
    volume
        Integer-labelled validation volume.
    cross_section_mask
        Two-dimensional voxelised capillary cross-section.
    """

    nz, ny, nx = shape

    volume = np.full(
        shape,
        fill_value=SOLID_LABEL,
        dtype=np.uint8,
    )

    centre_z = (nz - 1) / 2.0
    centre_y = (ny - 1) / 2.0

    z_coordinates, y_coordinates = np.ogrid[
        :nz,
        :ny,
    ]

    squared_radius = (
        (z_coordinates - centre_z) ** 2
        + (y_coordinates - centre_y) ** 2
    )

    cross_section_mask = (
        squared_radius
        <= float(radius_voxels) ** 2
    )

    three_dimensional_mask = np.broadcast_to(
        cross_section_mask[:, :, np.newaxis],
        shape,
    )

    volume[three_dimensional_mask] = PORE_LABEL

    return volume, cross_section_mask


def calculate_analytical_permeability(
    cross_section_mask,
    sample_shape,
):
    """
    Calculate the Hagen-Poiseuille sample-scale permeability.

    The equivalent radius is calculated from the actual voxelised pore area:

        r_equivalent = sqrt(A_pore / pi)

    This reduces the error caused by representing a circle using square
    voxels.
    """

    nz, ny, _ = sample_shape

    pore_cross_section_voxels = int(
        np.count_nonzero(cross_section_mask)
    )

    pore_cross_section_area_m2 = (
        pore_cross_section_voxels
        * VOXEL_SIZE_M**2
    )

    sample_cross_section_area_m2 = (
        nz
        * VOXEL_SIZE_M
        * ny
        * VOXEL_SIZE_M
    )

    equivalent_radius_m = np.sqrt(
        pore_cross_section_area_m2 / np.pi
    )

    analytical_porosity = (
        pore_cross_section_area_m2
        / sample_cross_section_area_m2
    )

    analytical_permeability_m2 = (
        np.pi
        * equivalent_radius_m**4
        / (
            8.0
            * sample_cross_section_area_m2
        )
    )

    alternative_permeability_m2 = (
        analytical_porosity
        * equivalent_radius_m**2
        / 8.0
    )

    if not np.isclose(
        analytical_permeability_m2,
        alternative_permeability_m2,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise RuntimeError(
            "The two analytical permeability expressions do not agree."
        )

    analytical_values = {
        "pore_cross_section_voxels": (
            pore_cross_section_voxels
        ),
        "pore_cross_section_area_m2": (
            pore_cross_section_area_m2
        ),
        "sample_cross_section_area_m2": (
            sample_cross_section_area_m2
        ),
        "equivalent_radius_m": equivalent_radius_m,
        "equivalent_radius_um": (
            equivalent_radius_m * 1.0e6
        ),
        "analytical_porosity": analytical_porosity,
        "analytical_permeability_m2": (
            analytical_permeability_m2
        ),
    }

    return analytical_values


def calculate_relative_error(
    calculated,
    expected,
):
    """Calculate absolute relative error."""

    calculated = float(calculated)
    expected = float(expected)

    if expected == 0.0:
        return abs(calculated - expected)

    return abs(calculated - expected) / abs(expected)


def create_validation_figure(
    analytical_permeability_m2,
    cemct_permeability_m2,
    relative_error_percent,
    output_path,
):
    """Create a publication-quality permeability comparison figure."""

    figure, axis = plt.subplots(
        figsize=(5.8, 4.8),
        dpi=150,
    )

    labels = [
        "Analytical\nHagen–Poiseuille",
        "CemCT\nSNOW2 + OpenPNM",
    ]

    values = [
        analytical_permeability_m2,
        cemct_permeability_m2,
    ]

    colours = [
        "#4C78A8",
        "#72B7B2",
    ]

    bars = axis.bar(
        labels,
        values,
        color=colours,
        edgecolor="black",
        linewidth=0.7,
        width=0.62,
    )

    for bar, value in zip(
        bars,
        values,
    ):
        axis.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height(),
            f"{value:.3e}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    axis.set_ylabel(
        r"Intrinsic permeability, $k$ (m$^2$)"
    )

    axis.set_title(
        "Analytical Permeability Validation",
        fontsize=12,
        fontweight="bold",
    )

    axis.text(
        0.5,
        0.94,
        f"Relative error = {relative_error_percent:.2f}%",
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )

    axis.grid(
        True,
        axis="y",
        linestyle="--",
        linewidth=0.6,
        alpha=0.4,
    )

    axis.set_axisbelow(
        True
    )

    axis.spines["top"].set_visible(
        False
    )

    axis.spines["right"].set_visible(
        False
    )

    maximum_value = max(values)

    if maximum_value > 0.0:
        axis.set_ylim(
            0.0,
            maximum_value * 1.25,
        )

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def main():
    """Run the analytical permeability validation."""

    volume, cross_section_mask = (
        create_cylindrical_x_channel(
            shape=(60, 60, 60),
            radius_voxels=10,
        )
    )

    analytical_values = (
        calculate_analytical_permeability(
            cross_section_mask=cross_section_mask,
            sample_shape=volume.shape,
        )
    )

    print("\nCemCT Analytical Permeability Validation")
    print("=" * 78)

    print(
        "Volume shape (Z, Y, X):",
        volume.shape,
    )

    print(
        "Equivalent capillary radius (um):",
        f"{analytical_values['equivalent_radius_um']:.6f}",
    )

    print(
        "Analytical porosity:",
        f"{analytical_values['analytical_porosity']:.6f}",
    )

    print(
        "Analytical permeability (m^2):",
        f"{analytical_values['analytical_permeability_m2']:.8e}",
    )

    print("\nRunning CemCT SNOW2/OpenPNM calculation...")
    print("This calculation may take several minutes.\n")

    permeability_results = analyse_permeability(
        volume=volume,
        phase_label=PORE_LABEL,
        voxel_size_um=VOXEL_SIZE_UM,
        directions=("X",),
        boundary_width=3,
        dynamic_viscosity_pa_s=(
            DYNAMIC_VISCOSITY_PA_S
        ),
        inlet_pressure_pa=INLET_PRESSURE_PA,
        outlet_pressure_pa=OUTLET_PRESSURE_PA,
        accuracy="standard",
        sigma=0.4,
        r_max=4,
        solver_tolerance=1.0e-10,
        maximum_iterations=5000,
    )

    summary = permeability_results["summary"]

    if len(summary) != 1:
        raise RuntimeError(
            "Expected exactly one X-direction permeability result."
        )

    result_row = summary.iloc[0]

    simulation_status = str(
        result_row["Status"]
    ).lower()

    if simulation_status not in {
        "solved",
        "calculated",
        "success",
        "successful",
    }:
        raise RuntimeError(
            "The X-direction permeability simulation did not succeed. "
            f"Returned status: {result_row['Status']}"
        )

    cemct_permeability_m2 = float(
        result_row["Intrinsic permeability (m^2)"]
    )

    flow_rate_m3_s = abs(
        float(result_row["Flow rate (m^3/s)"])
    )

    pressure_drop_pa = abs(
        float(result_row["Pressure drop (Pa)"])
    )

    simulation_length_m = float(
        result_row["Simulation length (m)"]
    )

    sample_area_m2 = float(
        result_row["Cross-sectional area (m^2)"]
    )

    darcy_reconstructed_permeability_m2 = (
        flow_rate_m3_s
        * DYNAMIC_VISCOSITY_PA_S
        * simulation_length_m
        / (
            sample_area_m2
            * pressure_drop_pa
        )
    )

    darcy_relative_error = calculate_relative_error(
        calculated=cemct_permeability_m2,
        expected=darcy_reconstructed_permeability_m2,
    )

    analytical_relative_error = calculate_relative_error(
        calculated=cemct_permeability_m2,
        expected=analytical_values[
            "analytical_permeability_m2"
        ],
    )

    darcy_passed = (
        darcy_relative_error
        <= DARCY_CONSISTENCY_TOLERANCE
    )

    analytical_passed = (
        analytical_relative_error
        <= CAPILLARY_RELATIVE_TOLERANCE
    )

    validation_table = pd.DataFrame(
        [
            {
                "validation_test": (
                    "Darcy-law consistency"
                ),
                "expected_permeability_m2": (
                    darcy_reconstructed_permeability_m2
                ),
                "cemct_permeability_m2": (
                    cemct_permeability_m2
                ),
                "absolute_error_m2": abs(
                    cemct_permeability_m2
                    - darcy_reconstructed_permeability_m2
                ),
                "relative_error": darcy_relative_error,
                "relative_error_percent": (
                    darcy_relative_error * 100.0
                ),
                "acceptance_tolerance": (
                    DARCY_CONSISTENCY_TOLERANCE
                ),
                "status": (
                    "PASS" if darcy_passed else "FAIL"
                ),
            },
            {
                "validation_test": (
                    "Hagen-Poiseuille capillary"
                ),
                "expected_permeability_m2": (
                    analytical_values[
                        "analytical_permeability_m2"
                    ]
                ),
                "cemct_permeability_m2": (
                    cemct_permeability_m2
                ),
                "absolute_error_m2": abs(
                    cemct_permeability_m2
                    - analytical_values[
                        "analytical_permeability_m2"
                    ]
                ),
                "relative_error": analytical_relative_error,
                "relative_error_percent": (
                    analytical_relative_error * 100.0
                ),
                "acceptance_tolerance": (
                    CAPILLARY_RELATIVE_TOLERANCE
                ),
                "status": (
                    "PASS"
                    if analytical_passed
                    else "FAIL"
                ),
            },
        ]
    )

    csv_output_path = (
        OUTPUT_DIRECTORY
        / "analytical_permeability_validation.csv"
    )

    figure_output_path = (
        OUTPUT_DIRECTORY
        / "analytical_permeability_validation.png"
    )

    validation_table.to_csv(
        csv_output_path,
        index=False,
    )

    create_validation_figure(
        analytical_permeability_m2=(
            analytical_values[
                "analytical_permeability_m2"
            ]
        ),
        cemct_permeability_m2=(
            cemct_permeability_m2
        ),
        relative_error_percent=(
            analytical_relative_error * 100.0
        ),
        output_path=figure_output_path,
    )

    print("\nCemCT permeability summary:")
    print(summary.to_string(index=False))

    print("\nValidation results:")
    print(validation_table.to_string(index=False))

    print(
        "\nDarcy-law consistency passed:",
        darcy_passed,
    )

    print(
        "Hagen-Poiseuille comparison passed:",
        analytical_passed,
    )

    print(
        "Analytical comparison relative error (%):",
        f"{analytical_relative_error * 100.0:.6f}",
    )

    print(
        "CSV saved to:",
        csv_output_path.resolve(),
    )

    print(
        "Figure saved to:",
        figure_output_path.resolve(),
    )

    if not darcy_passed:
        raise RuntimeError(
            "The Darcy-law consistency validation failed."
        )

    if not analytical_passed:
        raise RuntimeError(
            "The Hagen-Poiseuille analytical comparison failed. "
            "Do not modify the tolerance before investigating the "
            "voxelisation and pore-network assumptions."
        )

    print("\nAll analytical permeability validation cases passed: True")


if __name__ == "__main__":
    main()