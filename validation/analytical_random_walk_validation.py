"""
Analytical Random-Walk Validation for CemCT

Purpose
-------
This script validates the CemCT random-walk transport implementation using
two geometries with analytically expected diffusion tortuosity close to 1:

1. A fully open three-dimensional domain
2. A straight channel aligned with the X direction

The validation examines:

- Diffusion tortuosity
- Linearity of mean-square displacement
- Directional MSD fit R-squared
- Reproducibility using a deterministic random seed
- Compatibility with different MSD DataFrame column-name conventions

Expected result
---------------
For unobstructed diffusion:

    tau = 1

Small deviations are expected because the calculation uses a finite number
of walkers, a finite number of steps and a finite simulation domain.

Outputs
-------
1. CSV table containing validation results
2. PNG figure showing directional MSD results
3. Terminal summary reporting PASS or FAIL
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import matplotlib

# Use a non-interactive backend because this script is executed from CMD.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemct.random_walk import calculate_connected_random_walk


# =============================================================================
# Validation configuration
# =============================================================================

OUTPUT_DIRECTORY = Path("validation") / "results"

CSV_OUTPUT_PATH = (
    OUTPUT_DIRECTORY
    / "analytical_random_walk_validation.csv"
)

FIGURE_OUTPUT_PATH = (
    OUTPUT_DIRECTORY
    / "analytical_random_walk_validation.png"
)

PHASE_LABEL = 3

NUMBER_OF_WALKERS = 3000
NUMBER_OF_STEPS = 600
STRIDE = 10
NUMBER_OF_PROCESSES = 1

FIT_START_FRACTION = 0.05
FIT_END_FRACTION = 0.35

EXPECTED_TORTUOSITY = 1.0

# A 15% tolerance is appropriate for this finite stochastic validation.
MAXIMUM_RELATIVE_TORTUOSITY_ERROR = 0.15

# MSD should be strongly linear within the selected fitting interval.
MINIMUM_R_SQUARED = 0.97

# Identical deterministic runs should agree to numerical precision.
MAXIMUM_DETERMINISTIC_DIFFERENCE = 1.0e-12


# =============================================================================
# Synthetic validation geometries
# =============================================================================

def create_fully_open_domain(
    shape: tuple[int, int, int] = (96, 96, 96),
    phase_label: int = PHASE_LABEL,
) -> np.ndarray:
    """
    Create a fully open labelled three-dimensional domain.

    Array order is:

        axis 0 = Z
        axis 1 = Y
        axis 2 = X
    """
    return np.full(
        shape,
        fill_value=phase_label,
        dtype=np.uint8,
    )


def create_straight_x_channel(
    shape: tuple[int, int, int] = (31, 31, 160),
    channel_width: int = 9,
    phase_label: int = PHASE_LABEL,
) -> np.ndarray:
    """
    Create a straight square channel aligned with the X direction.

    Voxels inside the channel contain ``phase_label``. All other voxels
    contain label 0.
    """
    if channel_width <= 0:
        raise ValueError("Channel width must be greater than zero.")

    z_size, y_size, _ = shape

    if channel_width > min(z_size, y_size):
        raise ValueError(
            "Channel width cannot exceed the transverse volume dimensions."
        )

    volume = np.zeros(
        shape,
        dtype=np.uint8,
    )

    z_start = (z_size - channel_width) // 2
    y_start = (y_size - channel_width) // 2

    z_end = z_start + channel_width
    y_end = y_start + channel_width

    volume[
        z_start:z_end,
        y_start:y_end,
        :,
    ] = phase_label

    return volume


# =============================================================================
# Result compatibility helpers
# =============================================================================

def normalise_name(value: Any) -> str:
    """
    Convert a name to a lower-case alphanumeric comparison string.
    """
    return "".join(
        character.lower()
        for character in str(value)
        if character.isalnum()
    )


def resolve_dataframe_column(
    table: pd.DataFrame,
    candidate_names: tuple[str, ...],
) -> str:
    """
    Resolve a DataFrame column without requiring exact capitalisation.
    """
    normalised_columns = {
        normalise_name(column): column
        for column in table.columns
    }

    for candidate_name in candidate_names:
        normalised_candidate = normalise_name(candidate_name)

        if normalised_candidate in normalised_columns:
            return normalised_columns[normalised_candidate]

    raise KeyError(
        "Could not resolve a required DataFrame column.\n"
        f"Candidates: {list(candidate_names)}\n"
        f"Available columns: {table.columns.tolist()}"
    )


def resolve_time_column(
    msd_table: pd.DataFrame,
) -> str:
    """
    Resolve the random-walk time or step column.
    """
    return resolve_dataframe_column(
        table=msd_table,
        candidate_names=(
            "time_step",
            "time",
            "step",
            "steps",
            "saved_step",
            "saved_time_step",
        ),
    )


def resolve_axial_msd_column(
    msd_table: pd.DataFrame,
    direction: str,
) -> str:
    """
    Resolve a directional axial-MSD column.

    This function supports different historical or capitalisation conventions,
    including:

        scaled_axial_msd_x
        scaled_axial_msd_X
        axial_msd_x
        X_msd
        msd_x
    """
    direction_lower = str(direction).strip().lower()

    candidate_names = (
        f"scaled_axial_msd_{direction_lower}",
        f"scaled_axial_msd_{direction_lower.upper()}",
        f"axial_msd_{direction_lower}",
        f"axial_msd_{direction_lower.upper()}",
        f"scaled_msd_{direction_lower}",
        f"directional_msd_{direction_lower}",
        f"msd_{direction_lower}",
        f"{direction_lower}_msd",
        f"{direction_lower.upper()}_msd",
    )

    normalised_columns = {
        normalise_name(column): column
        for column in msd_table.columns
    }

    for candidate_name in candidate_names:
        normalised_candidate = normalise_name(candidate_name)

        if normalised_candidate in normalised_columns:
            return normalised_columns[normalised_candidate]

    # Controlled fallback for package versions using another descriptive name.
    matching_columns: list[str] = []

    for column in msd_table.columns:
        normalised_column = normalise_name(column)

        if "msd" not in normalised_column:
            continue

        if (
            normalised_column.endswith(direction_lower)
            or normalised_column.startswith(direction_lower)
        ):
            matching_columns.append(column)

    if len(matching_columns) == 1:
        return matching_columns[0]

    raise KeyError(
        f"Could not identify the {direction.upper()}-direction axial-MSD "
        "column.\n"
        f"Available columns: {msd_table.columns.tolist()}"
    )


def resolve_mapping_key(
    mapping: Mapping[str, Any],
    requested_key: str,
) -> str:
    """
    Resolve a mapping key without requiring exact capitalisation.
    """
    requested_normalised = normalise_name(requested_key)

    for existing_key in mapping:
        if normalise_name(existing_key) == requested_normalised:
            return existing_key

    raise KeyError(
        f"Could not find key '{requested_key}'. "
        f"Available keys: {list(mapping.keys())}"
    )


def get_directional_fit(
    result: Mapping[str, Any],
    direction: str,
) -> dict[str, float]:
    """
    Extract directional slope, tortuosity and R-squared.

    The function first checks ``fit_results`` and then falls back to fields
    stored in ``summary``.
    """
    direction_upper = str(direction).strip().upper()

    fit_results = result.get("fit_results", {})

    if isinstance(fit_results, Mapping):
        try:
            direction_key = resolve_mapping_key(
                fit_results,
                direction_upper,
            )

            directional_fit = fit_results[direction_key]

            if isinstance(directional_fit, Mapping):
                return {
                    "slope": float(
                        directional_fit.get("slope", np.nan)
                    ),
                    "tortuosity": float(
                        directional_fit.get("tortuosity", np.nan)
                    ),
                    "r_squared": float(
                        directional_fit.get("r_squared", np.nan)
                    ),
                    "fit_start_index": float(
                        directional_fit.get(
                            "fit_start_index",
                            np.nan,
                        )
                    ),
                    "fit_end_index": float(
                        directional_fit.get(
                            "fit_end_index",
                            np.nan,
                        )
                    ),
                }

        except KeyError:
            pass

    summary = result.get("summary", {})

    if not isinstance(summary, Mapping):
        raise TypeError(
            "Random-walk result 'summary' must be a mapping."
        )

    return {
        "slope": float(
            summary.get(
                f"{direction_upper}_msd_slope",
                np.nan,
            )
        ),
        "tortuosity": float(
            summary.get(
                f"{direction_upper}_tortuosity",
                np.nan,
            )
        ),
        "r_squared": float(
            summary.get(
                f"{direction_upper}_r_squared",
                np.nan,
            )
        ),
        "fit_start_index": np.nan,
        "fit_end_index": np.nan,
    }


# =============================================================================
# Random-walk execution
# =============================================================================

def run_random_walk_case(
    volume: np.ndarray,
    connectivity_mode: str,
) -> dict[str, Any]:
    """
    Run a deterministic CemCT random-walk calculation.
    """
    return calculate_connected_random_walk(
        volume=volume,
        phase_label=PHASE_LABEL,
        connectivity_mode=connectivity_mode,
        number_of_steps=NUMBER_OF_STEPS,
        number_of_walkers=NUMBER_OF_WALKERS,
        stride=STRIDE,
        same_start=False,
        number_of_processes=NUMBER_OF_PROCESSES,
        deterministic_seed=True,
        fit_start_fraction=FIT_START_FRACTION,
        fit_end_fraction=FIT_END_FRACTION,
        return_coordinates=False,
    )


def calculate_relative_error(
    measured_value: float,
    expected_value: float,
) -> float:
    """
    Calculate absolute relative error.
    """
    if expected_value == 0:
        return np.nan

    return abs(
        measured_value - expected_value
    ) / abs(expected_value)


def evaluate_direction(
    case_name: str,
    result: Mapping[str, Any],
    direction: str,
) -> dict[str, Any]:
    """
    Evaluate one directional random-walk result.
    """
    fit = get_directional_fit(
        result=result,
        direction=direction,
    )

    tortuosity = fit["tortuosity"]
    r_squared = fit["r_squared"]
    slope = fit["slope"]

    relative_error = calculate_relative_error(
        measured_value=tortuosity,
        expected_value=EXPECTED_TORTUOSITY,
    )

    tortuosity_passed = bool(
        np.isfinite(tortuosity)
        and relative_error
        <= MAXIMUM_RELATIVE_TORTUOSITY_ERROR
    )

    linearity_passed = bool(
        np.isfinite(r_squared)
        and r_squared >= MINIMUM_R_SQUARED
    )

    overall_passed = bool(
        tortuosity_passed
        and linearity_passed
    )

    return {
        "case": case_name,
        "direction": direction.upper(),
        "expected_tortuosity": EXPECTED_TORTUOSITY,
        "calculated_tortuosity": tortuosity,
        "relative_tortuosity_error": relative_error,
        "maximum_allowed_relative_error": (
            MAXIMUM_RELATIVE_TORTUOSITY_ERROR
        ),
        "msd_slope": slope,
        "r_squared": r_squared,
        "minimum_required_r_squared": MINIMUM_R_SQUARED,
        "tortuosity_passed": tortuosity_passed,
        "linearity_passed": linearity_passed,
        "overall_passed": overall_passed,
    }


# =============================================================================
# Deterministic-seed validation
# =============================================================================

def compare_deterministic_results(
    first_result: Mapping[str, Any],
    second_result: Mapping[str, Any],
    directions: tuple[str, ...],
) -> tuple[bool, float]:
    """
    Compare two calculations performed using the deterministic seed.
    """
    first_msd = first_result["msd"]
    second_msd = second_result["msd"]

    if not isinstance(first_msd, pd.DataFrame):
        raise TypeError("The first MSD result must be a DataFrame.")

    if not isinstance(second_msd, pd.DataFrame):
        raise TypeError("The second MSD result must be a DataFrame.")

    maximum_difference = 0.0

    for direction in directions:
        first_column = resolve_axial_msd_column(
            msd_table=first_msd,
            direction=direction,
        )

        second_column = resolve_axial_msd_column(
            msd_table=second_msd,
            direction=direction,
        )

        first_values = first_msd[
            first_column
        ].to_numpy(dtype=float)

        second_values = second_msd[
            second_column
        ].to_numpy(dtype=float)

        if first_values.shape != second_values.shape:
            return False, np.inf

        if first_values.size == 0:
            continue

        difference = float(
            np.nanmax(
                np.abs(
                    first_values - second_values
                )
            )
        )

        maximum_difference = max(
            maximum_difference,
            difference,
        )

    passed = bool(
        np.isfinite(maximum_difference)
        and maximum_difference
        <= MAXIMUM_DETERMINISTIC_DIFFERENCE
    )

    return passed, maximum_difference


# =============================================================================
# Plotting
# =============================================================================

def plot_msd_case(
    axis: plt.Axes,
    result: Mapping[str, Any],
    directions: tuple[str, ...],
    title: str,
) -> None:
    """
    Plot directional MSD curves using automatically resolved column names.
    """
    msd_table = result["msd"]

    if not isinstance(msd_table, pd.DataFrame):
        raise TypeError(
            "Random-walk result 'msd' must be a pandas DataFrame."
        )

    time_column = resolve_time_column(msd_table)

    time_values = msd_table[
        time_column
    ].to_numpy(dtype=float)

    colours = {
        "X": "#D95F59",
        "Y": "#4C78A8",
        "Z": "#59A14F",
    }

    for direction in directions:
        msd_column = resolve_axial_msd_column(
            msd_table=msd_table,
            direction=direction,
        )

        msd_values = msd_table[
            msd_column
        ].to_numpy(dtype=float)

        fit = get_directional_fit(
            result=result,
            direction=direction,
        )

        label = (
            f"{direction.upper()}: "
            f"$\\tau$={fit['tortuosity']:.3f}, "
            f"$R^2$={fit['r_squared']:.4f}"
        )

        axis.plot(
            time_values,
            msd_values,
            marker="o",
            markersize=3.5,
            linewidth=1.4,
            color=colours.get(
                direction.upper(),
                None,
            ),
            label=label,
        )

        # Calculate and display an independent linear fit over the same
        # configured fitting interval.
        number_of_points = len(time_values)

        fit_start_index = max(
            0,
            int(
                np.floor(
                    FIT_START_FRACTION
                    * number_of_points
                )
            ),
        )

        fit_end_index = min(
            number_of_points,
            max(
                fit_start_index + 2,
                int(
                    np.ceil(
                        FIT_END_FRACTION
                        * number_of_points
                    )
                ),
            ),
        )

        fit_time = time_values[
            fit_start_index:fit_end_index
        ]

        fit_msd = msd_values[
            fit_start_index:fit_end_index
        ]

        finite_values = (
            np.isfinite(fit_time)
            & np.isfinite(fit_msd)
        )

        if np.count_nonzero(finite_values) >= 2:
            linear_coefficients = np.polyfit(
                fit_time[finite_values],
                fit_msd[finite_values],
                deg=1,
            )

            fitted_values = np.polyval(
                linear_coefficients,
                fit_time,
            )

            axis.plot(
                fit_time,
                fitted_values,
                linestyle="--",
                linewidth=1.2,
                color=colours.get(
                    direction.upper(),
                    None,
                ),
                alpha=0.8,
            )

    axis.set_title(
        title,
        fontsize=11,
        fontweight="bold",
    )

    axis.set_xlabel("Time step")

    axis.set_ylabel(
        "Scaled directional mean-square displacement"
    )

    axis.grid(
        True,
        which="major",
        linestyle="--",
        linewidth=0.6,
        alpha=0.4,
    )

    axis.legend(
        frameon=False,
        fontsize=8,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def create_validation_figure(
    open_result: Mapping[str, Any],
    channel_result: Mapping[str, Any],
    validation_table: pd.DataFrame,
    deterministic_passed: bool,
    deterministic_maximum_difference: float,
    output_path: Path,
) -> None:
    """
    Create and save the complete random-walk validation figure.
    """
    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(12, 8),
    )

    plot_msd_case(
        axis=axes[0, 0],
        result=open_result,
        directions=("X", "Y", "Z"),
        title="Fully Open Domain",
    )

    plot_msd_case(
        axis=axes[0, 1],
        result=channel_result,
        directions=("X",),
        title="Straight X Channel",
    )

    plotted_rows = validation_table.copy()

    plotted_rows["label"] = (
        plotted_rows["case"]
        + "\n"
        + plotted_rows["direction"]
    )

    bar_colours = [
        "#59A14F" if passed else "#E15759"
        for passed
        in plotted_rows["overall_passed"]
    ]

    axes[1, 0].bar(
        plotted_rows["label"],
        plotted_rows["calculated_tortuosity"],
        color=bar_colours,
        edgecolor="#333333",
        linewidth=0.8,
    )

    axes[1, 0].axhline(
        EXPECTED_TORTUOSITY,
        color="#222222",
        linestyle="--",
        linewidth=1.2,
        label="Analytical value",
    )

    axes[1, 0].set_ylabel("Diffusion tortuosity, $\\tau$")

    axes[1, 0].set_title(
        "Analytical Tortuosity Comparison",
        fontsize=11,
        fontweight="bold",
    )

    axes[1, 0].grid(
        True,
        axis="y",
        linestyle="--",
        linewidth=0.6,
        alpha=0.4,
    )

    axes[1, 0].legend(
        frameon=False,
        fontsize=8,
    )

    axes[1, 0].spines["top"].set_visible(False)
    axes[1, 0].spines["right"].set_visible(False)

    axes[1, 1].axis("off")

    overall_rows_passed = bool(
        validation_table["overall_passed"].all()
    )

    overall_status = (
        overall_rows_passed
        and deterministic_passed
    )

    status_text = (
        "Random-Walk Validation Summary\n\n"
        f"Expected tortuosity: {EXPECTED_TORTUOSITY:.3f}\n"
        "Maximum relative tortuosity error: "
        f"{MAXIMUM_RELATIVE_TORTUOSITY_ERROR:.1%}\n"
        f"Minimum MSD $R^2$: {MINIMUM_R_SQUARED:.3f}\n\n"
        "Deterministic seed reproducibility: "
        f"{'PASS' if deterministic_passed else 'FAIL'}\n"
        "Maximum repeat difference: "
        f"{deterministic_maximum_difference:.3e}\n\n"
        f"Overall status: {'PASS' if overall_status else 'FAIL'}"
    )

    axes[1, 1].text(
        0.05,
        0.95,
        status_text,
        transform=axes[1, 1].transAxes,
        horizontalalignment="left",
        verticalalignment="top",
        fontsize=10,
        linespacing=1.5,
        bbox={
            "boxstyle": "round,pad=0.7",
            "facecolor": (
                "#E8F5E9"
                if overall_status
                else "#FFEBEE"
            ),
            "edgecolor": (
                "#59A14F"
                if overall_status
                else "#E15759"
            ),
        },
    )

    figure.suptitle(
        "CemCT Analytical Random-Walk Validation",
        fontsize=15,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=(0, 0, 1, 0.96),
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)


# =============================================================================
# Terminal reporting
# =============================================================================

def print_validation_summary(
    validation_table: pd.DataFrame,
    deterministic_passed: bool,
    deterministic_maximum_difference: float,
) -> None:
    """
    Print a concise validation summary.
    """
    display_columns = [
        "case",
        "direction",
        "expected_tortuosity",
        "calculated_tortuosity",
        "relative_tortuosity_error",
        "msd_slope",
        "r_squared",
        "overall_passed",
    ]

    print("\nAnalytical random-walk validation results")
    print("=" * 80)

    print(
        validation_table[
            display_columns
        ].to_string(index=False)
    )

    print("\nDeterministic-seed validation")
    print("-" * 80)

    print(
        "Maximum difference between repeated runs:",
        f"{deterministic_maximum_difference:.6e}",
    )

    print(
        "Deterministic repeat status:",
        "PASS" if deterministic_passed else "FAIL",
    )


# =============================================================================
# Main validation workflow
# =============================================================================

def main() -> None:
    """
    Run the complete analytical random-walk validation.
    """
    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("Creating fully open domain...")

    fully_open_volume = create_fully_open_domain()

    print("Running first fully open-domain calculation...")

    open_result_first = run_random_walk_case(
        volume=fully_open_volume,
        connectivity_mode="any",
    )

    print("Running repeated deterministic open-domain calculation...")

    open_result_second = run_random_walk_case(
        volume=fully_open_volume,
        connectivity_mode="any",
    )

    print("Creating straight X channel...")

    straight_channel_volume = create_straight_x_channel()

    print("Running straight X-channel calculation...")

    channel_result = run_random_walk_case(
        volume=straight_channel_volume,
        connectivity_mode="X",
    )

    print("\nFully open-domain MSD columns:")

    print(
        open_result_first["msd"].columns.tolist()
    )

    print("\nStraight-channel MSD columns:")

    print(
        channel_result["msd"].columns.tolist()
    )

    validation_rows: list[dict[str, Any]] = []

    for direction in ("X", "Y", "Z"):
        validation_rows.append(
            evaluate_direction(
                case_name="Fully open domain",
                result=open_result_first,
                direction=direction,
            )
        )

    validation_rows.append(
        evaluate_direction(
            case_name="Straight X channel",
            result=channel_result,
            direction="X",
        )
    )

    deterministic_passed, maximum_difference = (
        compare_deterministic_results(
            first_result=open_result_first,
            second_result=open_result_second,
            directions=("X", "Y", "Z"),
        )
    )

    validation_table = pd.DataFrame(
        validation_rows
    )

    validation_table[
        "deterministic_repeat_passed"
    ] = deterministic_passed

    validation_table[
        "deterministic_maximum_difference"
    ] = maximum_difference

    validation_table.to_csv(
        CSV_OUTPUT_PATH,
        index=False,
    )

    create_validation_figure(
        open_result=open_result_first,
        channel_result=channel_result,
        validation_table=validation_table,
        deterministic_passed=deterministic_passed,
        deterministic_maximum_difference=maximum_difference,
        output_path=FIGURE_OUTPUT_PATH,
    )

    print_validation_summary(
        validation_table=validation_table,
        deterministic_passed=deterministic_passed,
        deterministic_maximum_difference=maximum_difference,
    )

    print("\nOutput files")
    print("-" * 80)
    print("CSV:", CSV_OUTPUT_PATH.resolve())
    print("PNG:", FIGURE_OUTPUT_PATH.resolve())

    directional_results_passed = bool(
        validation_table["overall_passed"].all()
    )

    overall_validation_passed = bool(
        directional_results_passed
        and deterministic_passed
    )

    print(
        "\nOverall validation status:",
        "PASS" if overall_validation_passed else "FAIL",
    )

    if not overall_validation_passed:
        raise RuntimeError(
            "Analytical random-walk validation failed. "
            "Inspect the CSV table and PNG figure before changing any "
            "acceptance criteria."
        )


if __name__ == "__main__":
    main()