"""
Finite-Difference Diffusion Solver Sensitivity Analysis for CemCT

Purpose
-------
This validation script evaluates the numerical sensitivity of CemCT's
finite-difference diffusion analysis to:

    1. Solver tolerance
    2. Maximum solver iterations

A deterministic synthetic three-dimensional porous domain is generated
internally. The pore phase is represented by label 3 and contains a guaranteed
face-connected pathway spanning the X direction.

Only X-direction diffusion is evaluated so that changes in the numerical
results can be attributed to solver settings rather than directional
differences in the pore structure.

The highest-accuracy configuration is used as the reference:

    solver_tolerance = 1e-10
    maximum_iterations = 5000

Outputs
-------
    1. CSV table containing all sensitivity results
    2. PNG figure summarising convergence and runtime
    3. Markdown validation summary
    4. Synthetic labelled TIFF used for the validation

Coordinate convention
---------------------
    Array axis 0 = Z direction
    Array axis 1 = Y direction
    Array axis 2 = X direction
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

from cemct.diffusion import calculate_directional_diffusion


# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIRECTORY = PROJECT_ROOT / "validation" / "results"

CSV_OUTPUT_PATH = (
    RESULTS_DIRECTORY
    / "finite_difference_sensitivity_results.csv"
)

FIGURE_OUTPUT_PATH = (
    RESULTS_DIRECTORY
    / "finite_difference_sensitivity_analysis.png"
)

SUMMARY_OUTPUT_PATH = (
    RESULTS_DIRECTORY
    / "finite_difference_sensitivity_summary.md"
)

TIFF_OUTPUT_PATH = (
    RESULTS_DIRECTORY
    / "finite_difference_sensitivity_domain.tif"
)

PHASE_LABEL = 3
DIRECTION = "X"

# The voxel size does not change the dimensionless diffusion results,
# but it is included for complete metadata.
VOXEL_SIZE_UM = 0.7

# Synthetic volume dimensions in NumPy (Z, Y, X) order.
VOLUME_SHAPE = (40, 40, 56)

# Deterministic seed ensures reproducibility.
RANDOM_SEED = 20260726

# Tolerance sensitivity at fixed maximum iterations.
TOLERANCE_VALUES = (
    1.0e-4,
    1.0e-6,
    1.0e-8,
    1.0e-10,
)

FIXED_MAXIMUM_ITERATIONS = 5000

# Iteration sensitivity at fixed solver tolerance.
ITERATION_VALUES = (
    250,
    500,
    1000,
    5000,
)

FIXED_SOLVER_TOLERANCE = 1.0e-10

# Reference configuration.
REFERENCE_TOLERANCE = 1.0e-10
REFERENCE_MAXIMUM_ITERATIONS = 5000

# Difference below this level is classified as converged.
CONVERGENCE_THRESHOLD_PERCENT = 1.0


# =============================================================================
# Synthetic validation domain
# =============================================================================

def create_synthetic_diffusion_domain(
    shape: tuple[int, int, int] = VOLUME_SHAPE,
    random_seed: int = RANDOM_SEED,
    phase_label: int = PHASE_LABEL,
) -> np.ndarray:
    """
    Create a deterministic, non-trivial, X-through-connected porous domain.

    The domain combines:

        1. A smooth random pore structure
        2. A guaranteed sinusoidal X-spanning pore backbone
        3. Connected inlet and outlet regions

    Parameters
    ----------
    shape
        Volume shape in (Z, Y, X) order.
    random_seed
        Seed used to generate the reproducible random field.
    phase_label
        Integer label assigned to the pore phase.

    Returns
    -------
    numpy.ndarray
        Integer-labelled 3D validation volume.
    """
    if len(shape) != 3:
        raise ValueError("The validation volume must be three-dimensional.")

    if min(shape) < 20:
        raise ValueError(
            "Every validation-volume dimension must be at least 20 voxels."
        )

    rng = np.random.default_rng(random_seed)

    random_field = rng.normal(
        loc=0.0,
        scale=1.0,
        size=shape,
    )

    smooth_field = gaussian_filter(
        random_field,
        sigma=2.2,
        mode="reflect",
    )

    threshold = np.quantile(
        smooth_field,
        0.58,
    )

    pore_mask = smooth_field >= threshold

    z_size, y_size, x_size = shape

    z_coordinates, y_coordinates = np.ogrid[
        :z_size,
        :y_size,
    ]

    channel_radius = max(
        3,
        min(z_size, y_size) // 9,
    )

    # Add a smoothly varying channel through the X direction.
    # Consecutive cross-sections overlap, ensuring face connectivity.
    for x_index in range(x_size):
        phase = 2.0 * np.pi * x_index / max(x_size - 1, 1)

        centre_z = int(
            round(
                z_size / 2
                + 0.16 * z_size * np.sin(phase)
            )
        )

        centre_y = int(
            round(
                y_size / 2
                + 0.14 * y_size * np.sin(
                    2.0 * phase + np.pi / 5.0
                )
            )
        )

        channel_cross_section = (
            (z_coordinates - centre_z) ** 2
            + (y_coordinates - centre_y) ** 2
            <= channel_radius**2
        )

        pore_mask[:, :, x_index] |= channel_cross_section

    # Widen inlet and outlet regions to produce stable boundary conditions.
    boundary_width = 2
    inlet_outlet_radius = channel_radius + 1

    for x_index in range(boundary_width):
        phase = 2.0 * np.pi * x_index / max(x_size - 1, 1)

        centre_z = int(
            round(
                z_size / 2
                + 0.16 * z_size * np.sin(phase)
            )
        )

        centre_y = int(
            round(
                y_size / 2
                + 0.14 * y_size * np.sin(
                    2.0 * phase + np.pi / 5.0
                )
            )
        )

        inlet_cross_section = (
            (z_coordinates - centre_z) ** 2
            + (y_coordinates - centre_y) ** 2
            <= inlet_outlet_radius**2
        )

        pore_mask[:, :, x_index] |= inlet_cross_section

    for x_index in range(
        x_size - boundary_width,
        x_size,
    ):
        phase = 2.0 * np.pi * x_index / max(x_size - 1, 1)

        centre_z = int(
            round(
                z_size / 2
                + 0.16 * z_size * np.sin(phase)
            )
        )

        centre_y = int(
            round(
                y_size / 2
                + 0.14 * y_size * np.sin(
                    2.0 * phase + np.pi / 5.0
                )
            )
        )

        outlet_cross_section = (
            (z_coordinates - centre_z) ** 2
            + (y_coordinates - centre_y) ** 2
            <= inlet_outlet_radius**2
        )

        pore_mask[:, :, x_index] |= outlet_cross_section

    labelled_volume = np.zeros(
        shape,
        dtype=np.uint8,
    )

    labelled_volume[pore_mask] = np.uint8(phase_label)

    return labelled_volume


# =============================================================================
# Result handling
# =============================================================================

def normalise_column_name(column_name: Any) -> str:
    """Return a lower-case underscore-separated column name."""
    return (
        str(column_name)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
    )


def find_column(
    dataframe: pd.DataFrame,
    candidates: tuple[str, ...],
) -> str:
    """
    Find a DataFrame column using several possible names.

    This makes the validation script tolerant of small changes in the
    presentation names used by the CemCT result table.
    """
    normalised_columns = {
        normalise_column_name(column): column
        for column in dataframe.columns
    }

    for candidate in candidates:
        normalised_candidate = normalise_column_name(candidate)

        if normalised_candidate in normalised_columns:
            return normalised_columns[normalised_candidate]

    raise KeyError(
        "Could not find any of the expected columns "
        f"{candidates}. Available columns: "
        f"{dataframe.columns.tolist()}"
    )


def extract_summary_values(
    result: dict[str, Any],
) -> dict[str, Any]:
    """Extract the X-direction numerical values from a CemCT result."""
    if "summary" not in result:
        raise KeyError(
            "The diffusion result does not contain a 'summary' table."
        )

    summary = result["summary"]

    if not isinstance(summary, pd.DataFrame):
        raise TypeError(
            "result['summary'] must be a pandas DataFrame."
        )

    if summary.empty:
        raise ValueError(
            "The diffusion summary table is empty."
        )

    direction_column = find_column(
        summary,
        (
            "direction",
            "Direction",
        ),
    )

    direction_rows = summary[
        summary[direction_column]
        .astype(str)
        .str.upper()
        .eq(DIRECTION)
    ]

    if direction_rows.empty:
        raise ValueError(
            f"No {DIRECTION}-direction row was found "
            "in the diffusion summary."
        )

    row = direction_rows.iloc[0]

    status_column = find_column(
        summary,
        (
            "status",
            "Status",
        ),
    )

    formation_factor_column = find_column(
        summary,
        (
            "formation_factor",
            "Formation factor",
        ),
    )

    relative_diffusivity_column = find_column(
        summary,
        (
            "relative_effective_diffusivity",
            "Relative effective diffusivity",
            "Deff / D0",
        ),
    )

    tortuosity_column = find_column(
        summary,
        (
            "diffusion_tortuosity",
            "Diffusion tortuosity",
        ),
    )

    runtime_column = find_column(
        summary,
        (
            "run_time_seconds",
            "runtime_seconds",
            "Run time seconds",
        ),
    )

    effective_porosity_column = find_column(
        summary,
        (
            "effective_porosity",
            "Effective porosity",
        ),
    )

    connectivity_column = find_column(
        summary,
        (
            "directional_connectivity",
            "Directional connectivity",
        ),
    )

    return {
        "status": str(row[status_column]),
        "formation_factor": float(
            row[formation_factor_column]
        ),
        "relative_effective_diffusivity": float(
            row[relative_diffusivity_column]
        ),
        "diffusion_tortuosity": float(
            row[tortuosity_column]
        ),
        "reported_run_time_seconds": float(
            row[runtime_column]
        ),
        "effective_porosity": float(
            row[effective_porosity_column]
        ),
        "directional_connectivity": float(
            row[connectivity_column]
        ),
    }


def run_diffusion_configuration(
    volume: np.ndarray,
    solver_tolerance: float,
    maximum_iterations: int,
) -> dict[str, Any]:
    """Run one X-direction finite-difference configuration."""
    print(
        "\nRunning configuration:"
        f" tolerance={solver_tolerance:.1e},"
        f" maximum_iterations={maximum_iterations}"
    )

    wall_clock_start = time.perf_counter()

    result = calculate_directional_diffusion(
        volume=volume,
        phase_label=PHASE_LABEL,
        directions=(DIRECTION,),
        solver=None,
        solver_tolerance=solver_tolerance,
        maximum_iterations=maximum_iterations,
    )

    wall_clock_seconds = (
        time.perf_counter() - wall_clock_start
    )

    extracted_values = extract_summary_values(result)

    status = extracted_values["status"].strip().lower()

    successful_statuses = {
        "calculated",
        "solved",
        "success",
        "successful",
    }

    if status not in successful_statuses:
        raise RuntimeError(
            "The diffusion calculation did not report a successful status. "
            f"Returned status: {extracted_values['status']}"
        )

    formation_factor = extracted_values["formation_factor"]
    relative_diffusivity = extracted_values[
        "relative_effective_diffusivity"
    ]

    reciprocal_consistency_error = abs(
        formation_factor * relative_diffusivity - 1.0
    )

    return {
        "solver_tolerance": float(solver_tolerance),
        "maximum_iterations": int(maximum_iterations),
        "status": extracted_values["status"],
        "effective_porosity": extracted_values[
            "effective_porosity"
        ],
        "directional_connectivity": extracted_values[
            "directional_connectivity"
        ],
        "formation_factor": formation_factor,
        "relative_effective_diffusivity": relative_diffusivity,
        "diffusion_tortuosity": extracted_values[
            "diffusion_tortuosity"
        ],
        "reported_run_time_seconds": extracted_values[
            "reported_run_time_seconds"
        ],
        "wall_clock_seconds": float(wall_clock_seconds),
        "reciprocal_consistency_error": float(
            reciprocal_consistency_error
        ),
    }


def calculate_relative_difference_percent(
    value: float,
    reference_value: float,
) -> float:
    """Calculate absolute relative difference from a reference value."""
    if not np.isfinite(value):
        return float("nan")

    if not np.isfinite(reference_value):
        return float("nan")

    if reference_value == 0.0:
        return float("nan")

    return float(
        abs(value - reference_value)
        / abs(reference_value)
        * 100.0
    )


# =============================================================================
# Validation execution
# =============================================================================

def run_sensitivity_analysis(
    volume: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Run tolerance and maximum-iteration sensitivity studies.

    Duplicate numerical configurations are calculated only once.
    """
    configuration_cache: dict[
        tuple[float, int],
        dict[str, Any],
    ] = {}

    result_rows: list[dict[str, Any]] = []

    def get_configuration_result(
        tolerance: float,
        maximum_iterations: int,
    ) -> dict[str, Any]:
        cache_key = (
            float(tolerance),
            int(maximum_iterations),
        )

        if cache_key not in configuration_cache:
            configuration_cache[cache_key] = (
                run_diffusion_configuration(
                    volume=volume,
                    solver_tolerance=tolerance,
                    maximum_iterations=maximum_iterations,
                )
            )

        return configuration_cache[cache_key].copy()

    # -------------------------------------------------------------------------
    # Tolerance sensitivity
    # -------------------------------------------------------------------------

    for tolerance in TOLERANCE_VALUES:
        configuration_result = get_configuration_result(
            tolerance=tolerance,
            maximum_iterations=FIXED_MAXIMUM_ITERATIONS,
        )

        result_rows.append(
            {
                "study": "solver_tolerance",
                "varied_parameter": "solver_tolerance",
                "varied_value": float(tolerance),
                **configuration_result,
            }
        )

    # -------------------------------------------------------------------------
    # Maximum-iteration sensitivity
    # -------------------------------------------------------------------------

    for maximum_iterations in ITERATION_VALUES:
        configuration_result = get_configuration_result(
            tolerance=FIXED_SOLVER_TOLERANCE,
            maximum_iterations=maximum_iterations,
        )

        result_rows.append(
            {
                "study": "maximum_iterations",
                "varied_parameter": "maximum_iterations",
                "varied_value": int(maximum_iterations),
                **configuration_result,
            }
        )

    results = pd.DataFrame(result_rows)

    reference_configuration = get_configuration_result(
        tolerance=REFERENCE_TOLERANCE,
        maximum_iterations=REFERENCE_MAXIMUM_ITERATIONS,
    )

    reference_tortuosity = reference_configuration[
        "diffusion_tortuosity"
    ]

    reference_formation_factor = reference_configuration[
        "formation_factor"
    ]

    reference_relative_diffusivity = reference_configuration[
        "relative_effective_diffusivity"
    ]

    results[
        "tortuosity_relative_difference_percent"
    ] = results["diffusion_tortuosity"].apply(
        lambda value: calculate_relative_difference_percent(
            value,
            reference_tortuosity,
        )
    )

    results[
        "formation_factor_relative_difference_percent"
    ] = results["formation_factor"].apply(
        lambda value: calculate_relative_difference_percent(
            value,
            reference_formation_factor,
        )
    )

    results[
        "relative_diffusivity_difference_percent"
    ] = results["relative_effective_diffusivity"].apply(
        lambda value: calculate_relative_difference_percent(
            value,
            reference_relative_diffusivity,
        )
    )

    results["within_one_percent_of_reference"] = (
        results[
            "tortuosity_relative_difference_percent"
        ]
        <= CONVERGENCE_THRESHOLD_PERCENT
    )

    return results, reference_configuration


# =============================================================================
# Figure
# =============================================================================

def create_sensitivity_figure(
    results: pd.DataFrame,
    reference_configuration: dict[str, Any],
    output_path: Path,
) -> None:
    """Create the finite-difference sensitivity summary figure."""
    tolerance_results = (
        results[
            results["study"].eq("solver_tolerance")
        ]
        .sort_values("solver_tolerance")
        .copy()
    )

    iteration_results = (
        results[
            results["study"].eq("maximum_iterations")
        ]
        .sort_values("maximum_iterations")
        .copy()
    )

    reference_tortuosity = reference_configuration[
        "diffusion_tortuosity"
    ]

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(12.5, 8.5),
    )

    # -------------------------------------------------------------------------
    # Panel A: tortuosity versus tolerance
    # -------------------------------------------------------------------------

    axes[0, 0].plot(
        tolerance_results["solver_tolerance"],
        tolerance_results["diffusion_tortuosity"],
        color="#D9534F",
        marker="o",
        linewidth=1.8,
        markersize=6,
    )

    axes[0, 0].axhline(
        reference_tortuosity,
        color="#333333",
        linestyle="--",
        linewidth=1.2,
        label="Reference",
    )

    axes[0, 0].set_xscale("log")
    axes[0, 0].invert_xaxis()

    axes[0, 0].set_xlabel("Solver tolerance")
    axes[0, 0].set_ylabel("Diffusion tortuosity, τ")
    axes[0, 0].set_title(
        "(a) Tortuosity sensitivity to solver tolerance"
    )
    axes[0, 0].legend(frameon=False)

    # -------------------------------------------------------------------------
    # Panel B: relative difference versus tolerance
    # -------------------------------------------------------------------------

    axes[0, 1].plot(
        tolerance_results["solver_tolerance"],
        tolerance_results[
            "tortuosity_relative_difference_percent"
        ],
        color="#4472C4",
        marker="o",
        linewidth=1.8,
        markersize=6,
        label="Tortuosity",
    )

    axes[0, 1].plot(
        tolerance_results["solver_tolerance"],
        tolerance_results[
            "formation_factor_relative_difference_percent"
        ],
        color="#ED7D31",
        marker="s",
        linewidth=1.6,
        markersize=5,
        label="Formation factor",
    )

    axes[0, 1].axhline(
        CONVERGENCE_THRESHOLD_PERCENT,
        color="#333333",
        linestyle="--",
        linewidth=1.1,
        label="1% criterion",
    )

    axes[0, 1].set_xscale("log")
    axes[0, 1].invert_xaxis()

    axes[0, 1].set_xlabel("Solver tolerance")
    axes[0, 1].set_ylabel(
        "Relative difference from reference (%)"
    )
    axes[0, 1].set_title(
        "(b) Numerical difference from reference"
    )
    axes[0, 1].legend(frameon=False)

    # -------------------------------------------------------------------------
    # Panel C: tortuosity versus maximum iterations
    # -------------------------------------------------------------------------

    axes[1, 0].plot(
        iteration_results["maximum_iterations"],
        iteration_results["diffusion_tortuosity"],
        color="#70AD47",
        marker="o",
        linewidth=1.8,
        markersize=6,
    )

    axes[1, 0].axhline(
        reference_tortuosity,
        color="#333333",
        linestyle="--",
        linewidth=1.2,
        label="Reference",
    )

    axes[1, 0].set_xscale("log")

    axes[1, 0].set_xlabel("Maximum iterations")
    axes[1, 0].set_ylabel("Diffusion tortuosity, τ")
    axes[1, 0].set_title(
        "(c) Tortuosity sensitivity to iteration limit"
    )
    axes[1, 0].legend(frameon=False)

    # -------------------------------------------------------------------------
    # Panel D: runtime
    # -------------------------------------------------------------------------

    tolerance_labels = [
        f"tol={value:.0e}"
        for value in tolerance_results["solver_tolerance"]
    ]

    iteration_labels = [
        f"iter={int(value)}"
        for value in iteration_results["maximum_iterations"]
    ]

    runtime_labels = tolerance_labels + iteration_labels

    runtime_values = np.concatenate(
        [
            tolerance_results[
                "wall_clock_seconds"
            ].to_numpy(dtype=float),
            iteration_results[
                "wall_clock_seconds"
            ].to_numpy(dtype=float),
        ]
    )

    runtime_colours = (
        ["#D9534F"] * len(tolerance_labels)
        + ["#70AD47"] * len(iteration_labels)
    )

    bar_positions = np.arange(len(runtime_labels))

    axes[1, 1].bar(
        bar_positions,
        runtime_values,
        color=runtime_colours,
        edgecolor="#333333",
        linewidth=0.6,
    )

    axes[1, 1].set_xticks(bar_positions)
    axes[1, 1].set_xticklabels(
        runtime_labels,
        rotation=45,
        ha="right",
    )

    axes[1, 1].set_xlabel("Solver configuration")
    axes[1, 1].set_ylabel("Wall-clock time (s)")
    axes[1, 1].set_title(
        "(d) Computational runtime"
    )

    for axis in axes.ravel():
        axis.grid(
            True,
            which="both",
            linestyle="--",
            linewidth=0.55,
            alpha=0.35,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        "CemCT Finite-Difference Diffusion Solver Sensitivity",
        fontsize=15,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=(0.0, 0.0, 1.0, 0.96)
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)


# =============================================================================
# Markdown report
# =============================================================================

def write_markdown_summary(
    results: pd.DataFrame,
    reference_configuration: dict[str, Any],
    labelled_volume: np.ndarray,
    output_path: Path,
) -> None:
    """Write a human-readable Markdown sensitivity summary."""
    tolerance_results = (
        results[
            results["study"].eq("solver_tolerance")
        ]
        .sort_values(
            "solver_tolerance",
            ascending=False,
        )
        .copy()
    )

    iteration_results = (
        results[
            results["study"].eq("maximum_iterations")
        ]
        .sort_values("maximum_iterations")
        .copy()
    )

    pore_fraction = float(
        np.count_nonzero(
            labelled_volume == PHASE_LABEL
        )
        / labelled_volume.size
    )

    all_statuses_successful = (
        results["status"]
        .astype(str)
        .str.lower()
        .isin(
            {
                "calculated",
                "solved",
                "success",
                "successful",
            }
        )
        .all()
    )

    maximum_tortuosity_difference = float(
        results[
            "tortuosity_relative_difference_percent"
        ].max()
    )

    converged_configuration_count = int(
        results["within_one_percent_of_reference"].sum()
    )

    total_configuration_count = int(len(results))

    reciprocal_error_maximum = float(
        results[
            "reciprocal_consistency_error"
        ].max()
    )

    validation_status = (
        "PASS"
        if all_statuses_successful
        and np.isfinite(maximum_tortuosity_difference)
        and reciprocal_error_maximum <= 1.0e-5
        else "REVIEW"
    )

    lines = [
        "# Finite-Difference Diffusion Solver Sensitivity",
        "",
        f"**Validation status:** {validation_status}",
        "",
        "## Purpose",
        "",
        (
            "This validation evaluates the numerical sensitivity of "
            "CemCT's X-direction finite-difference diffusion analysis "
            "to solver tolerance and maximum iteration count."
        ),
        "",
        "## Validation domain",
        "",
        f"- Volume shape (Z, Y, X): `{labelled_volume.shape}`",
        f"- Selected phase label: `{PHASE_LABEL}`",
        f"- Selected phase fraction: `{pore_fraction:.6f}`",
        f"- Direction: `{DIRECTION}`",
        f"- Random seed: `{RANDOM_SEED}`",
        (
            "- Connectivity definition: strict three-dimensional "
            "6-connectivity"
        ),
        "",
        "## Reference configuration",
        "",
        (
            f"- Solver tolerance: `{REFERENCE_TOLERANCE:.1e}`"
        ),
        (
            "- Maximum iterations: "
            f"`{REFERENCE_MAXIMUM_ITERATIONS}`"
        ),
        (
            "- Reference diffusion tortuosity: "
            f"`{reference_configuration['diffusion_tortuosity']:.10g}`"
        ),
        (
            "- Reference formation factor: "
            f"`{reference_configuration['formation_factor']:.10g}`"
        ),
        (
            "- Reference relative effective diffusivity: "
            f"`{reference_configuration['relative_effective_diffusivity']:.10g}`"
        ),
        "",
        "## Tolerance sensitivity",
        "",
        (
            "| Solver tolerance | Diffusion tortuosity | "
            "Formation factor | Deff / D0 | "
            "Tortuosity difference (%) | Runtime (s) |"
        ),
        (
            "|---:|---:|---:|---:|---:|---:|"
        ),
    ]

    for _, row in tolerance_results.iterrows():
        lines.append(
            "| "
            f"{row['solver_tolerance']:.1e} | "
            f"{row['diffusion_tortuosity']:.8g} | "
            f"{row['formation_factor']:.8g} | "
            f"{row['relative_effective_diffusivity']:.8g} | "
            f"{row['tortuosity_relative_difference_percent']:.6g} | "
            f"{row['wall_clock_seconds']:.4f} |"
        )

    lines.extend(
        [
            "",
            "## Maximum-iteration sensitivity",
            "",
            (
                "| Maximum iterations | Diffusion tortuosity | "
                "Formation factor | Deff / D0 | "
                "Tortuosity difference (%) | Runtime (s) |"
            ),
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for _, row in iteration_results.iterrows():
        lines.append(
            "| "
            f"{int(row['maximum_iterations'])} | "
            f"{row['diffusion_tortuosity']:.8g} | "
            f"{row['formation_factor']:.8g} | "
            f"{row['relative_effective_diffusivity']:.8g} | "
            f"{row['tortuosity_relative_difference_percent']:.6g} | "
            f"{row['wall_clock_seconds']:.4f} |"
        )

    lines.extend(
        [
            "",
            "## Numerical consistency",
            "",
            (
                "- Successful configurations: "
                f"`{int(all_statuses_successful)}`"
            ),
            (
                "- Configurations within 1% of the reference: "
                f"`{converged_configuration_count}/"
                f"{total_configuration_count}`"
            ),
            (
                "- Maximum tortuosity difference from the reference: "
                f"`{maximum_tortuosity_difference:.6g}%`"
            ),
            (
                "- Maximum reciprocal consistency error "
                "`abs(F × Deff/D0 - 1)`: "
                f"`{reciprocal_error_maximum:.6e}`"
            ),
            "",
            "## Interpretation",
            "",
            (
                "Configurations with a diffusion-tortuosity difference "
                f"below {CONVERGENCE_THRESHOLD_PERCENT:.1f}% are "
                "classified as numerically converged relative to the "
                "highest-accuracy reference configuration."
            ),
            "",
            (
                "Runtime measurements are wall-clock measurements and "
                "may vary with operating system, processor load and "
                "available memory."
            ),
            "",
            "## Output files",
            "",
            f"- Numerical results: `{CSV_OUTPUT_PATH.name}`",
            f"- Summary figure: `{FIGURE_OUTPUT_PATH.name}`",
            f"- Validation volume: `{TIFF_OUTPUT_PATH.name}`",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    """Run the complete finite-difference sensitivity validation."""
    RESULTS_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "Creating deterministic finite-difference "
        "validation domain..."
    )

    labelled_volume = create_synthetic_diffusion_domain()

    tifffile.imwrite(
        TIFF_OUTPUT_PATH,
        labelled_volume,
        photometric="minisblack",
    )

    pore_voxel_count = int(
        np.count_nonzero(
            labelled_volume == PHASE_LABEL
        )
    )

    pore_fraction = (
        pore_voxel_count / labelled_volume.size
    )

    print(f"Volume shape: {labelled_volume.shape}")
    print(f"Selected phase label: {PHASE_LABEL}")
    print(f"Selected pore voxels: {pore_voxel_count}")
    print(f"Selected pore fraction: {pore_fraction:.6f}")

    print(
        "\nRunning finite-difference solver "
        "sensitivity analysis..."
    )

    results, reference_configuration = (
        run_sensitivity_analysis(
            volume=labelled_volume,
        )
    )

    results.to_csv(
        CSV_OUTPUT_PATH,
        index=False,
    )

    create_sensitivity_figure(
        results=results,
        reference_configuration=reference_configuration,
        output_path=FIGURE_OUTPUT_PATH,
    )

    write_markdown_summary(
        results=results,
        reference_configuration=reference_configuration,
        labelled_volume=labelled_volume,
        output_path=SUMMARY_OUTPUT_PATH,
    )

    print("\nSensitivity analysis completed.")
    print(f"CSV: {CSV_OUTPUT_PATH}")
    print(f"Figure: {FIGURE_OUTPUT_PATH}")
    print(f"Summary: {SUMMARY_OUTPUT_PATH}")
    print(f"Validation TIFF: {TIFF_OUTPUT_PATH}")

    print("\nResults:")
    print(
        results[
            [
                "study",
                "solver_tolerance",
                "maximum_iterations",
                "status",
                "formation_factor",
                "relative_effective_diffusivity",
                "diffusion_tortuosity",
                "tortuosity_relative_difference_percent",
                "wall_clock_seconds",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()