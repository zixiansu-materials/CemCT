"""
Analytical diffusion validation for CemCT.

This script compares CemCT finite-difference diffusion results against
analytical solutions for ideal straight transport domains.

Validation cases
----------------
1. Fully open domain

       Effective porosity = 1.0
       Deff / D0          = 1.0
       Formation factor   = 1.0
       Diffusion tortuosity = 1.0

2. Straight X-direction channel occupying 25% of the sample cross-section

       Effective porosity = 0.25
       Deff / D0          = 0.25
       Formation factor   = 4.0
       Diffusion tortuosity = 1.0

Array convention
----------------
    axis 0 = Z
    axis 1 = Y
    axis 2 = X

The validation calculations are intentionally maintained outside the regular
unit-test suite because finite-difference simulations are more computationally
expensive than ordinary unit tests.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemct.diffusion import calculate_directional_diffusion


OUTPUT_DIRECTORY = Path(__file__).parent / "results"
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

PORE_LABEL = 3

RELATIVE_TOLERANCE = 0.02
ABSOLUTE_TOLERANCE = 1.0e-6


def create_fully_open_domain(
    shape=(20, 20, 20),
):
    """
    Create a completely open diffusion domain.

    Every voxel belongs to the transported phase.
    """

    return np.full(
        shape,
        fill_value=PORE_LABEL,
        dtype=np.uint8,
    )


def create_straight_x_channel(
    shape=(20, 20, 20),
):
    """
    Create a straight X-direction channel.

    The channel cross-section is 10 x 10 voxels inside a 20 x 20
    sample cross-section. Therefore, its exact effective porosity is:

        epsilon = (10 x 10) / (20 x 20) = 0.25

    The channel is straight and has no geometrical tortuosity.
    """

    volume = np.ones(
        shape,
        dtype=np.uint8,
    )

    volume[5:15, 5:15, :] = PORE_LABEL

    return volume


def calculate_relative_error(
    calculated,
    expected,
):
    """Calculate relative error while safely handling zero."""

    calculated = float(calculated)
    expected = float(expected)

    if expected == 0.0:
        return abs(calculated - expected)

    return abs(calculated - expected) / abs(expected)


def create_validation_record(
    case_name,
    direction,
    metric,
    expected,
    calculated,
    run_time_seconds,
):
    """Create one analytical validation record."""

    expected = float(expected)
    calculated = float(calculated)

    absolute_error = abs(calculated - expected)

    relative_error = calculate_relative_error(
        calculated=calculated,
        expected=expected,
    )

    passed = np.isclose(
        calculated,
        expected,
        rtol=RELATIVE_TOLERANCE,
        atol=ABSOLUTE_TOLERANCE,
    )

    return {
        "validation_case": case_name,
        "direction": direction,
        "metric": metric,
        "expected": expected,
        "calculated": calculated,
        "absolute_error": absolute_error,
        "relative_error": relative_error,
        "relative_error_percent": relative_error * 100.0,
        "relative_tolerance": RELATIVE_TOLERANCE,
        "absolute_tolerance": ABSOLUTE_TOLERANCE,
        "run_time_seconds": float(run_time_seconds),
        "status": "PASS" if passed else "FAIL",
    }


def extract_validation_records(
    case_name,
    diffusion_results,
    expected_values,
):
    """Extract expected and calculated values from one CemCT result."""

    summary = diffusion_results["summary"]

    records = []

    for _, row in summary.iterrows():
        direction = str(row["direction"]).upper()
        simulation_status = str(row["status"]).lower()

        if simulation_status not in {
            "calculated",
            "solved",
            "success",
            "successful",
        }:
            raise RuntimeError(
                f"{case_name}, direction {direction}, "
                f"returned status: {row['status']}"
            )

        run_time_seconds = float(
            row.get(
                "run_time_seconds",
                np.nan,
            )
        )

        metric_columns = {
            "effective_porosity": "effective_porosity",
            "relative_effective_diffusivity": (
                "relative_effective_diffusivity"
            ),
            "formation_factor": "formation_factor",
            "diffusion_tortuosity": "diffusion_tortuosity",
        }

        for metric_name, column_name in metric_columns.items():
            records.append(
                create_validation_record(
                    case_name=case_name,
                    direction=direction,
                    metric=metric_name,
                    expected=expected_values[metric_name],
                    calculated=row[column_name],
                    run_time_seconds=run_time_seconds,
                )
            )

    return records


def run_fully_open_validation():
    """Run X-, Y-, and Z-direction validation for a fully open domain."""

    volume = create_fully_open_domain()

    results = calculate_directional_diffusion(
        volume=volume,
        phase_label=PORE_LABEL,
        directions=("X", "Y", "Z"),
        solver=None,
        solver_tolerance=1.0e-10,
        maximum_iterations=5000,
    )

    expected_values = {
        "effective_porosity": 1.0,
        "relative_effective_diffusivity": 1.0,
        "formation_factor": 1.0,
        "diffusion_tortuosity": 1.0,
    }

    records = extract_validation_records(
        case_name="Fully open domain",
        diffusion_results=results,
        expected_values=expected_values,
    )

    return records, results


def run_straight_channel_validation():
    """Run X-direction validation for an ideal straight channel."""

    volume = create_straight_x_channel()

    results = calculate_directional_diffusion(
        volume=volume,
        phase_label=PORE_LABEL,
        directions=("X",),
        solver=None,
        solver_tolerance=1.0e-10,
        maximum_iterations=5000,
    )

    expected_porosity = 0.25

    expected_values = {
        "effective_porosity": expected_porosity,
        "relative_effective_diffusivity": expected_porosity,
        "formation_factor": 1.0 / expected_porosity,
        "diffusion_tortuosity": 1.0,
    }

    records = extract_validation_records(
        case_name="Straight X channel",
        diffusion_results=results,
        expected_values=expected_values,
    )

    return records, results


def create_validation_figure(
    validation_table,
    output_path,
):
    """Create a publication-quality analytical validation figure."""

    plot_metrics = [
        (
            "relative_effective_diffusivity",
            r"Relative effective diffusivity, $D_\mathrm{eff}/D_0$",
        ),
        (
            "formation_factor",
            r"Formation factor, $F$",
        ),
        (
            "diffusion_tortuosity",
            r"Diffusion tortuosity, $\tau$",
        ),
    ]

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(12.5, 4.2),
        dpi=150,
    )

    case_order = [
        "Fully open domain — X",
        "Fully open domain — Y",
        "Fully open domain — Z",
        "Straight X channel — X",
    ]

    colours = [
        "#4C78A8",
        "#72B7B2",
    ]

    for axis, (metric_name, ylabel) in zip(
        axes,
        plot_metrics,
    ):
        metric_table = validation_table.loc[
            validation_table["metric"].eq(metric_name)
        ].copy()

        metric_table["case_direction"] = (
            metric_table["validation_case"]
            + " — "
            + metric_table["direction"]
        )

        metric_table["case_direction"] = pd.Categorical(
            metric_table["case_direction"],
            categories=case_order,
            ordered=True,
        )

        metric_table = metric_table.sort_values(
            "case_direction"
        )

        x_positions = np.arange(
            len(metric_table)
        )

        bar_width = 0.36

        axis.bar(
            x_positions - bar_width / 2,
            metric_table["expected"],
            width=bar_width,
            label="Analytical",
            color=colours[0],
            edgecolor="black",
            linewidth=0.6,
        )

        axis.bar(
            x_positions + bar_width / 2,
            metric_table["calculated"],
            width=bar_width,
            label="CemCT",
            color=colours[1],
            edgecolor="black",
            linewidth=0.6,
        )

        axis.set_xticks(
            x_positions
        )

        axis.set_xticklabels(
            [
                "Open\nX",
                "Open\nY",
                "Open\nZ",
                "Channel\nX",
            ],
            rotation=0,
        )

        axis.set_ylabel(
            ylabel
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

    axes[0].legend(
        frameon=False,
        loc="upper right",
    )

    figure.suptitle(
        "Analytical Validation of Directional Diffusion",
        fontsize=13,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=(0, 0, 1, 0.94)
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


def main():
    """Run all analytical diffusion validation cases."""

    validation_records = []

    print("\nRunning fully open-domain validation...")

    open_records, open_results = (
        run_fully_open_validation()
    )

    validation_records.extend(
        open_records
    )

    print("Running straight X-channel validation...")

    channel_records, channel_results = (
        run_straight_channel_validation()
    )

    validation_records.extend(
        channel_records
    )

    validation_table = pd.DataFrame(
        validation_records
    )

    csv_output_path = (
        OUTPUT_DIRECTORY
        / "analytical_diffusion_validation.csv"
    )

    figure_output_path = (
        OUTPUT_DIRECTORY
        / "analytical_diffusion_validation.png"
    )

    validation_table.to_csv(
        csv_output_path,
        index=False,
    )

    create_validation_figure(
        validation_table=validation_table,
        output_path=figure_output_path,
    )

    print("\nCemCT Analytical Diffusion Validation")
    print("=" * 100)

    display_columns = [
        "validation_case",
        "direction",
        "metric",
        "expected",
        "calculated",
        "relative_error_percent",
        "status",
    ]

    print(
        validation_table[
            display_columns
        ].to_string(
            index=False
        )
    )

    print("=" * 100)

    all_passed = (
        validation_table["status"]
        .eq("PASS")
        .all()
    )

    maximum_relative_error_percent = (
        validation_table[
            "relative_error_percent"
        ].max()
    )

    print(
        "\nAll analytical diffusion validation cases passed:",
        all_passed,
    )

    print(
        "Maximum relative error (%):",
        f"{maximum_relative_error_percent:.6f}",
    )

    print(
        "CSV saved to:",
        csv_output_path.resolve(),
    )

    print(
        "Figure saved to:",
        figure_output_path.resolve(),
    )

    if not all_passed:
        failed_rows = validation_table.loc[
            validation_table["status"].eq("FAIL")
        ]

        print("\nFailed validation records:")
        print(
            failed_rows[
                display_columns
            ].to_string(
                index=False
            )
        )

        raise RuntimeError(
            "One or more analytical diffusion validation cases failed."
        )


if __name__ == "__main__":
    main()