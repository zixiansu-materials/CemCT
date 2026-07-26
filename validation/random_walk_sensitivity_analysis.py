"""
Random-Walk Numerical Sensitivity and Convergence Analysis for CemCT

Purpose
-------
This script evaluates the sensitivity of CemCT random-walk transport results
to two principal numerical parameters:

    1. Number of walkers
    2. Number of random-walk steps

A fully open three-dimensional domain is used. Its expected diffusion
tortuosity is approximately 1 because no solid obstacles obstruct transport.

The analysis investigates:

    - Mean diffusion tortuosity
    - Directional X, Y, and Z tortuosity
    - Mean MSD fitting R-squared
    - Runtime
    - Difference relative to the highest-resolution reference calculation

The calculations use a deterministic random seed to improve reproducibility.

Outputs
-------
    1. CSV table of all sensitivity-analysis results
    2. PNG summary figure
    3. Markdown summary report

Notes
-----
This script performs numerical sensitivity analysis. It does not replace the
analytical random-walk validation. Lower-resolution calculations are allowed
to show larger statistical variation.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Mapping

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemct.random_walk import calculate_connected_random_walk


# =============================================================================
# Configuration
# =============================================================================

PHASE_LABEL = 3

DOMAIN_SHAPE = (64, 64, 64)

WALKER_COUNTS = (
    500,
    1000,
    3000,
    5000,
)

FIXED_STEPS_FOR_WALKER_STUDY = 600

STEP_COUNTS = (
    200,
    400,
    600,
    1000,
)

FIXED_WALKERS_FOR_STEP_STUDY = 3000

STRIDE = 10
NUMBER_OF_PROCESSES = 1

FIT_START_FRACTION = 0.05
FIT_END_FRACTION = 0.35

EXPECTED_TORTUOSITY = 1.0

MAXIMUM_REFERENCE_TORTUOSITY_ERROR_PERCENT = 15.0
MINIMUM_REFERENCE_R_SQUARED = 0.97
MAXIMUM_FINAL_CONVERGENCE_CHANGE_PERCENT = 10.0


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUT_DIRECTORY = (
    PROJECT_ROOT
    / "validation"
    / "results"
)

CSV_OUTPUT_PATH = (
    OUTPUT_DIRECTORY
    / "random_walk_sensitivity_results.csv"
)

FIGURE_OUTPUT_PATH = (
    OUTPUT_DIRECTORY
    / "random_walk_sensitivity_analysis.png"
)

REPORT_OUTPUT_PATH = (
    OUTPUT_DIRECTORY
    / "random_walk_sensitivity_summary.md"
)


# =============================================================================
# Utility functions
# =============================================================================

def normalise_name(value: str) -> str:
    """
    Normalise a mapping key for robust comparison.
    """
    return "".join(
        character.lower()
        for character in str(value)
        if character.isalnum()
    )


def get_mapping_value(
    mapping: Mapping[str, Any],
    candidate_names: tuple[str, ...],
    default: float = np.nan,
) -> float:
    """
    Extract a numeric value using several possible key names.
    """
    normalised_keys = {
        normalise_name(key): key
        for key in mapping
    }

    for candidate_name in candidate_names:
        normalised_candidate = normalise_name(candidate_name)

        if normalised_candidate in normalised_keys:
            original_key = normalised_keys[normalised_candidate]
            value = mapping[original_key]

            try:
                return float(value)
            except (TypeError, ValueError):
                continue

    return float(default)


def get_directional_value(
    summary: Mapping[str, Any],
    direction: str,
    property_name: str,
) -> float:
    """
    Extract a directional result from the random-walk summary.
    """
    return get_mapping_value(
        summary,
        (
            f"{direction}_{property_name}",
            f"{direction.lower()}_{property_name}",
            f"{direction} {property_name}",
            f"{direction.lower()} {property_name}",
        ),
    )


def calculate_directional_mean(
    summary: Mapping[str, Any],
    property_name: str,
) -> float:
    """
    Calculate the mean of valid X, Y, and Z directional values.
    """
    values = np.asarray(
        [
            get_directional_value(
                summary,
                direction,
                property_name,
            )
            for direction in ("X", "Y", "Z")
        ],
        dtype=float,
    )

    valid_values = values[np.isfinite(values)]

    if valid_values.size == 0:
        return np.nan

    return float(np.mean(valid_values))


def extract_random_walk_metrics(
    result: Mapping[str, Any],
) -> dict[str, float]:
    """
    Extract the principal random-walk metrics.
    """
    summary = result["summary"]

    x_tortuosity = get_directional_value(
        summary,
        "X",
        "tortuosity",
    )

    y_tortuosity = get_directional_value(
        summary,
        "Y",
        "tortuosity",
    )

    z_tortuosity = get_directional_value(
        summary,
        "Z",
        "tortuosity",
    )

    mean_tortuosity = get_mapping_value(
        summary,
        (
            "mean_tortuosity",
            "mean tortuosity",
        ),
    )

    if not np.isfinite(mean_tortuosity):
        mean_tortuosity = calculate_directional_mean(
            summary,
            "tortuosity",
        )

    x_r_squared = get_directional_value(
        summary,
        "X",
        "r_squared",
    )

    y_r_squared = get_directional_value(
        summary,
        "Y",
        "r_squared",
    )

    z_r_squared = get_directional_value(
        summary,
        "Z",
        "r_squared",
    )

    mean_r_squared = get_mapping_value(
        summary,
        (
            "mean_r_squared",
            "mean r squared",
            "mean_r2",
        ),
    )

    if not np.isfinite(mean_r_squared):
        mean_r_squared = calculate_directional_mean(
            summary,
            "r_squared",
        )

    return {
        "X_tortuosity": x_tortuosity,
        "Y_tortuosity": y_tortuosity,
        "Z_tortuosity": z_tortuosity,
        "mean_tortuosity": mean_tortuosity,
        "X_r_squared": x_r_squared,
        "Y_r_squared": y_r_squared,
        "Z_r_squared": z_r_squared,
        "mean_r_squared": mean_r_squared,
    }


# =============================================================================
# Random-walk calculation
# =============================================================================

def run_random_walk_case(
    volume: np.ndarray,
    number_of_walkers: int,
    number_of_steps: int,
) -> dict[str, Any]:
    """
    Run one deterministic random-walk calculation.
    """
    print(
        "\nRunning random-walk calculation:"
        f" walkers={number_of_walkers},"
        f" steps={number_of_steps}"
    )

    start_time = time.perf_counter()

    result = calculate_connected_random_walk(
        volume=volume,
        phase_label=PHASE_LABEL,
        connectivity_mode="any",
        number_of_steps=number_of_steps,
        number_of_walkers=number_of_walkers,
        stride=STRIDE,
        same_start=False,
        number_of_processes=NUMBER_OF_PROCESSES,
        deterministic_seed=True,
        fit_start_fraction=FIT_START_FRACTION,
        fit_end_fraction=FIT_END_FRACTION,
        return_coordinates=False,
    )

    runtime_seconds = time.perf_counter() - start_time

    metrics = extract_random_walk_metrics(result)

    metrics["runtime_seconds"] = float(runtime_seconds)

    print(
        "Completed:"
        f" mean tortuosity={metrics['mean_tortuosity']:.6f},"
        f" mean R-squared={metrics['mean_r_squared']:.6f},"
        f" runtime={runtime_seconds:.2f} s"
    )

    return {
        "result": result,
        "metrics": metrics,
    }


# =============================================================================
# Sensitivity study
# =============================================================================

def run_sensitivity_analysis(
    volume: np.ndarray,
) -> pd.DataFrame:
    """
    Run walker-count and step-count sensitivity studies.
    """
    calculation_cache: dict[
        tuple[int, int],
        dict[str, Any],
    ] = {}

    result_rows: list[dict[str, Any]] = []

    def get_or_run_case(
        number_of_walkers: int,
        number_of_steps: int,
    ) -> dict[str, Any]:
        cache_key = (
            int(number_of_walkers),
            int(number_of_steps),
        )

        if cache_key not in calculation_cache:
            calculation_cache[cache_key] = (
                run_random_walk_case(
                    volume=volume,
                    number_of_walkers=number_of_walkers,
                    number_of_steps=number_of_steps,
                )
            )
        else:
            print(
                "\nReusing previously calculated case:"
                f" walkers={number_of_walkers},"
                f" steps={number_of_steps}"
            )

        return calculation_cache[cache_key]

    print("\n" + "=" * 72)
    print("Walker-count sensitivity study")
    print("=" * 72)

    for number_of_walkers in WALKER_COUNTS:
        case = get_or_run_case(
            number_of_walkers=number_of_walkers,
            number_of_steps=(
                FIXED_STEPS_FOR_WALKER_STUDY
            ),
        )

        result_rows.append(
            {
                "study": "walker_count",
                "varied_parameter": (
                    "number_of_walkers"
                ),
                "varied_value": number_of_walkers,
                "number_of_walkers": number_of_walkers,
                "number_of_steps": (
                    FIXED_STEPS_FOR_WALKER_STUDY
                ),
                **case["metrics"],
            }
        )

    print("\n" + "=" * 72)
    print("Step-count sensitivity study")
    print("=" * 72)

    for number_of_steps in STEP_COUNTS:
        case = get_or_run_case(
            number_of_walkers=(
                FIXED_WALKERS_FOR_STEP_STUDY
            ),
            number_of_steps=number_of_steps,
        )

        result_rows.append(
            {
                "study": "step_count",
                "varied_parameter": (
                    "number_of_steps"
                ),
                "varied_value": number_of_steps,
                "number_of_walkers": (
                    FIXED_WALKERS_FOR_STEP_STUDY
                ),
                "number_of_steps": number_of_steps,
                **case["metrics"],
            }
        )

    results = pd.DataFrame(result_rows)

    results[
        "analytical_tortuosity_error_percent"
    ] = (
        np.abs(
            results["mean_tortuosity"]
            - EXPECTED_TORTUOSITY
        )
        / EXPECTED_TORTUOSITY
        * 100.0
    )

    results[
        "difference_from_reference_percent"
    ] = np.nan

    results[
        "change_from_previous_case_percent"
    ] = np.nan

    for study_name in (
        "walker_count",
        "step_count",
    ):
        study_mask = results["study"].eq(study_name)

        study_indices = (
            results.loc[study_mask]
            .sort_values("varied_value")
            .index
        )

        sorted_values = results.loc[
            study_indices,
            "mean_tortuosity",
        ]

        reference_value = float(
            sorted_values.iloc[-1]
        )

        results.loc[
            study_indices,
            "difference_from_reference_percent",
        ] = (
            np.abs(
                sorted_values - reference_value
            )
            / abs(reference_value)
            * 100.0
        )

        previous_change = (
            sorted_values
            .pct_change()
            .abs()
            * 100.0
        )

        results.loc[
            study_indices,
            "change_from_previous_case_percent",
        ] = previous_change

    results["analytical_acceptance"] = (
        results[
            "analytical_tortuosity_error_percent"
        ].le(
            MAXIMUM_REFERENCE_TORTUOSITY_ERROR_PERCENT
        )
        & results["mean_r_squared"].ge(
            MINIMUM_REFERENCE_R_SQUARED
        )
    )

    return results


# =============================================================================
# Figure
# =============================================================================

def create_sensitivity_figure(
    results: pd.DataFrame,
) -> None:
    """
    Create and save the numerical sensitivity figure.
    """
    walker_results = (
        results[
            results["study"].eq("walker_count")
        ]
        .sort_values("number_of_walkers")
        .copy()
    )

    step_results = (
        results[
            results["study"].eq("step_count")
        ]
        .sort_values("number_of_steps")
        .copy()
    )

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(12, 9),
    )

    axes[0, 0].plot(
        walker_results["number_of_walkers"],
        walker_results["mean_tortuosity"],
        marker="o",
        linewidth=1.8,
        color="#2166AC",
        label="Calculated",
    )

    axes[0, 0].axhline(
        EXPECTED_TORTUOSITY,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Expected value = 1",
    )

    axes[0, 0].set_xlabel(
        "Number of walkers"
    )

    axes[0, 0].set_ylabel(
        "Mean diffusion tortuosity"
    )

    axes[0, 0].set_title(
        "(a) Walker-count convergence"
    )

    axes[0, 0].legend(
        frameon=False
    )

    axes[0, 1].plot(
        walker_results["number_of_walkers"],
        walker_results["mean_r_squared"],
        marker="s",
        linewidth=1.8,
        color="#B2182B",
    )

    axes[0, 1].axhline(
        MINIMUM_REFERENCE_R_SQUARED,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label=(
            f"Acceptance = "
            f"{MINIMUM_REFERENCE_R_SQUARED:.2f}"
        ),
    )

    axes[0, 1].set_xlabel(
        "Number of walkers"
    )

    axes[0, 1].set_ylabel(
        r"Mean MSD fit $R^2$"
    )

    axes[0, 1].set_title(
        "(b) MSD fitting stability"
    )

    axes[0, 1].legend(
        frameon=False
    )

    axes[1, 0].plot(
        step_results["number_of_steps"],
        step_results["mean_tortuosity"],
        marker="o",
        linewidth=1.8,
        color="#1B7837",
        label="Calculated",
    )

    axes[1, 0].axhline(
        EXPECTED_TORTUOSITY,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Expected value = 1",
    )

    axes[1, 0].set_xlabel(
        "Number of random-walk steps"
    )

    axes[1, 0].set_ylabel(
        "Mean diffusion tortuosity"
    )

    axes[1, 0].set_title(
        "(c) Step-count convergence"
    )

    axes[1, 0].legend(
        frameon=False
    )

    runtime_labels = [
        (
            f"W={int(row.number_of_walkers)}"
            if row.study == "walker_count"
            else f"S={int(row.number_of_steps)}"
        )
        for row in results.itertuples()
    ]

    runtime_colours = [
        (
            "#4393C3"
            if study == "walker_count"
            else "#7FBF7B"
        )
        for study in results["study"]
    ]

    axes[1, 1].bar(
        runtime_labels,
        results["runtime_seconds"],
        color=runtime_colours,
        edgecolor="black",
        linewidth=0.6,
    )

    axes[1, 1].set_xlabel(
        "Sensitivity-analysis case"
    )

    axes[1, 1].set_ylabel(
        "Runtime (s)"
    )

    axes[1, 1].set_title(
        "(d) Computational runtime"
    )

    axes[1, 1].tick_params(
        axis="x",
        rotation=45,
    )

    for axis in axes.flat:
        axis.grid(
            True,
            linestyle="--",
            linewidth=0.6,
            alpha=0.35,
        )

        axis.set_axisbelow(True)

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        "CemCT Random-Walk Numerical Sensitivity Analysis",
        fontsize=15,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=(0, 0, 1, 0.96)
    )

    figure.savefig(
        FIGURE_OUTPUT_PATH,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)


# =============================================================================
# Summary report
# =============================================================================

def create_summary_report(
    results: pd.DataFrame,
) -> dict[str, Any]:
    """
    Evaluate final reference cases and write a Markdown report.
    """
    walker_results = (
        results[
            results["study"].eq("walker_count")
        ]
        .sort_values("number_of_walkers")
    )

    step_results = (
        results[
            results["study"].eq("step_count")
        ]
        .sort_values("number_of_steps")
    )

    final_walker_case = walker_results.iloc[-1]
    previous_walker_case = walker_results.iloc[-2]

    final_step_case = step_results.iloc[-1]
    previous_step_case = step_results.iloc[-2]

    walker_final_change_percent = (
        abs(
            final_walker_case["mean_tortuosity"]
            - previous_walker_case["mean_tortuosity"]
        )
        / abs(final_walker_case["mean_tortuosity"])
        * 100.0
    )

    step_final_change_percent = (
        abs(
            final_step_case["mean_tortuosity"]
            - previous_step_case["mean_tortuosity"]
        )
        / abs(final_step_case["mean_tortuosity"])
        * 100.0
    )

    walker_reference_passed = bool(
        final_walker_case[
            "analytical_acceptance"
        ]
        and walker_final_change_percent
        <= MAXIMUM_FINAL_CONVERGENCE_CHANGE_PERCENT
    )

    step_reference_passed = bool(
        final_step_case[
            "analytical_acceptance"
        ]
        and step_final_change_percent
        <= MAXIMUM_FINAL_CONVERGENCE_CHANGE_PERCENT
    )

    overall_passed = bool(
        walker_reference_passed
        and step_reference_passed
    )

    report_text = f"""# CemCT Random-Walk Sensitivity Analysis

## Configuration

- Domain shape: `{DOMAIN_SHAPE}`
- Selected phase label: `{PHASE_LABEL}`
- Expected open-domain tortuosity: `{EXPECTED_TORTUOSITY:.3f}`
- Deterministic seed: `True`
- MSD fit interval: `{FIT_START_FRACTION:.2f}` to `{FIT_END_FRACTION:.2f}`

## Walker-count reference case

- Number of walkers: `{int(final_walker_case["number_of_walkers"])}`
- Number of steps: `{int(final_walker_case["number_of_steps"])}`
- Mean tortuosity: `{final_walker_case["mean_tortuosity"]:.6f}`
- Mean R-squared: `{final_walker_case["mean_r_squared"]:.6f}`
- Analytical error: `{final_walker_case["analytical_tortuosity_error_percent"]:.3f}%`
- Change relative to preceding case: `{walker_final_change_percent:.3f}%`
- Status: `{"PASS" if walker_reference_passed else "REVIEW"}`

## Step-count reference case

- Number of walkers: `{int(final_step_case["number_of_walkers"])}`
- Number of steps: `{int(final_step_case["number_of_steps"])}`
- Mean tortuosity: `{final_step_case["mean_tortuosity"]:.6f}`
- Mean R-squared: `{final_step_case["mean_r_squared"]:.6f}`
- Analytical error: `{final_step_case["analytical_tortuosity_error_percent"]:.3f}%`
- Change relative to preceding case: `{step_final_change_percent:.3f}%`
- Status: `{"PASS" if step_reference_passed else "REVIEW"}`

## Overall assessment

`{"PASS" if overall_passed else "REVIEW REQUIRED"}`

A review status does not necessarily indicate a software error. It indicates
that the selected numerical parameters have not yet satisfied the predefined
convergence or analytical-tolerance criteria.
"""

    REPORT_OUTPUT_PATH.write_text(
        report_text,
        encoding="utf-8",
    )

    return {
        "walker_reference_passed": (
            walker_reference_passed
        ),
        "step_reference_passed": (
            step_reference_passed
        ),
        "overall_passed": overall_passed,
        "walker_final_change_percent": (
            walker_final_change_percent
        ),
        "step_final_change_percent": (
            step_final_change_percent
        ),
    }


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    """
    Run the complete sensitivity and convergence analysis.
    """
    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 72)
    print("CemCT random-walk sensitivity analysis")
    print("=" * 72)

    print(
        "Creating fully open validation domain:"
        f" shape={DOMAIN_SHAPE},"
        f" phase label={PHASE_LABEL}"
    )

    volume = np.full(
        DOMAIN_SHAPE,
        fill_value=PHASE_LABEL,
        dtype=np.uint8,
    )

    results = run_sensitivity_analysis(
        volume=volume,
    )

    results.to_csv(
        CSV_OUTPUT_PATH,
        index=False,
    )

    create_sensitivity_figure(
        results=results,
    )

    assessment = create_summary_report(
        results=results,
    )

    print("\n" + "=" * 72)
    print("Sensitivity-analysis results")
    print("=" * 72)

    display_columns = [
        "study",
        "number_of_walkers",
        "number_of_steps",
        "mean_tortuosity",
        "mean_r_squared",
        "runtime_seconds",
        "difference_from_reference_percent",
        "change_from_previous_case_percent",
    ]

    print(
        results[display_columns].to_string(
            index=False
        )
    )

    print("\nOutput files:")
    print(f"CSV:    {CSV_OUTPUT_PATH}")
    print(f"Figure: {FIGURE_OUTPUT_PATH}")
    print(f"Report: {REPORT_OUTPUT_PATH}")

    print("\nReference assessment:")
    print(
        "Walker-count reference:",
        (
            "PASS"
            if assessment["walker_reference_passed"]
            else "REVIEW"
        ),
    )

    print(
        "Step-count reference:",
        (
            "PASS"
            if assessment["step_reference_passed"]
            else "REVIEW"
        ),
    )

    print(
        "Overall:",
        (
            "PASS"
            if assessment["overall_passed"]
            else "REVIEW REQUIRED"
        ),
    )


if __name__ == "__main__":
    main()