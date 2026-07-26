"""
Permeability Parameter Sensitivity Analysis for CemCT

Purpose
-------
This script evaluates the sensitivity of pore-network-based intrinsic
permeability estimates to three SNOW2 and boundary-condition parameters:

    1. boundary_width
    2. sigma
    3. r_max

A one-factor-at-a-time approach is used. For each study, only one parameter
is varied while the remaining parameters are held at their baseline values.

The script evaluates:

    - Directional connectivity
    - Number of extracted network pores
    - Number of extracted network throats
    - Flow rate
    - Intrinsic permeability
    - Relative change from the baseline permeability
    - Wall-clock runtime

The analysis is performed in the X direction using a deterministic synthetic
labelled volume. Label 3 represents the pore phase.

Array convention
----------------
    Axis 0 = Z
    Axis 1 = Y
    Axis 2 = X

Outputs
-------
    1. Validation TIFF volume
    2. CSV table containing all sensitivity results
    3. PNG summary figure
    4. Markdown summary report
"""

from __future__ import annotations

import re
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

from cemct.permeability import analyse_permeability


# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RESULTS_DIRECTORY = PROJECT_ROOT / "validation" / "results"
RESULTS_DIRECTORY.mkdir(parents=True, exist_ok=True)

OUTPUT_STEM = "permeability_parameter_sensitivity"

VALIDATION_VOLUME_PATH = (
    RESULTS_DIRECTORY
    / f"{OUTPUT_STEM}_domain.tif"
)

RESULTS_CSV_PATH = (
    RESULTS_DIRECTORY
    / f"{OUTPUT_STEM}_results.csv"
)

FIGURE_PATH = (
    RESULTS_DIRECTORY
    / f"{OUTPUT_STEM}_analysis.png"
)

SUMMARY_PATH = (
    RESULTS_DIRECTORY
    / f"{OUTPUT_STEM}_summary.md"
)

PHASE_LABEL = 3
VOXEL_SIZE_UM = 0.7
DIRECTION = "X"

BASELINE_BOUNDARY_WIDTH = 3
BASELINE_SIGMA = 0.4
BASELINE_R_MAX = 4

SOLVER_TOLERANCE = 1.0e-10
MAXIMUM_ITERATIONS = 5000

DYNAMIC_VISCOSITY_PA_S = 1.0
INLET_PRESSURE_PA = 1.0
OUTLET_PRESSURE_PA = 0.0


# =============================================================================
# Synthetic validation volume
# =============================================================================

def create_deterministic_validation_volume(
    shape: tuple[int, int, int] = (40, 40, 56),
    phase_label: int = PHASE_LABEL,
    random_seed: int = 20260726,
) -> np.ndarray:
    """
    Create a deterministic labelled porous volume.

    The volume contains a spatially correlated pore structure together with
    several guaranteed X-direction through-channels.

    Parameters
    ----------
    shape
        Volume shape in (Z, Y, X) order.
    phase_label
        Integer label representing pore space.
    random_seed
        Random seed used to make the volume reproducible.

    Returns
    -------
    numpy.ndarray
        Labelled uint8 volume. Label 1 is solid and the selected phase label
        represents pore space.
    """
    if len(shape) != 3:
        raise ValueError("The validation volume must be three-dimensional.")

    if min(shape) < 16:
        raise ValueError(
            "Every validation-volume dimension must be at least 16 voxels."
        )

    rng = np.random.default_rng(random_seed)

    random_field = rng.normal(
        loc=0.0,
        scale=1.0,
        size=shape,
    )

    smooth_field = gaussian_filter(
        random_field,
        sigma=2.0,
        mode="reflect",
    )

    threshold = np.quantile(
        smooth_field,
        0.68,
    )

    pore_mask = smooth_field >= threshold

    z_coordinates, y_coordinates, _ = np.indices(shape)

    channel_definitions = (
        (
            int(round(shape[0] * 0.28)),
            int(round(shape[1] * 0.30)),
            3,
        ),
        (
            int(round(shape[0] * 0.50)),
            int(round(shape[1] * 0.52)),
            4,
        ),
        (
            int(round(shape[0] * 0.72)),
            int(round(shape[1] * 0.70)),
            3,
        ),
    )

    for z_centre, y_centre, radius in channel_definitions:
        channel_mask = (
            (z_coordinates - z_centre) ** 2
            + (y_coordinates - y_centre) ** 2
            <= radius**2
        )

        pore_mask |= channel_mask

    # Add several connector regions between the guaranteed channels and
    # surrounding pore clusters.
    connector_x_positions = (
        int(round(shape[2] * 0.25)),
        int(round(shape[2] * 0.50)),
        int(round(shape[2] * 0.75)),
    )

    for x_position in connector_x_positions:
        x_min = max(0, x_position - 2)
        x_max = min(shape[2], x_position + 3)

        pore_mask[
            int(shape[0] * 0.25):int(shape[0] * 0.75),
            int(shape[1] * 0.27):int(shape[1] * 0.73),
            x_min:x_max,
        ] |= smooth_field[
            int(shape[0] * 0.25):int(shape[0] * 0.75),
            int(shape[1] * 0.27):int(shape[1] * 0.73),
            x_min:x_max,
        ] >= np.quantile(smooth_field, 0.55)

    labelled_volume = np.ones(
        shape,
        dtype=np.uint8,
    )

    labelled_volume[pore_mask] = np.uint8(phase_label)

    return labelled_volume


# =============================================================================
# Result-column utilities
# =============================================================================

def normalise_column_name(value: Any) -> str:
    """Convert a DataFrame column name to a simplified comparison key."""
    return re.sub(
        pattern=r"[^a-z0-9]+",
        repl="",
        string=str(value).lower(),
    )


def get_row_value(
    row: pd.Series,
    candidate_names: tuple[str, ...],
    default: Any = np.nan,
) -> Any:
    """
    Retrieve a value using several possible column-name conventions.
    """
    normalised_columns = {
        normalise_column_name(column): column
        for column in row.index
    }

    for candidate_name in candidate_names:
        key = normalise_column_name(candidate_name)

        if key in normalised_columns:
            return row[normalised_columns[key]]

    return default


def status_is_successful(status: Any) -> bool:
    """Return True when the permeability calculation status is successful."""
    normalised_status = str(status).strip().lower()

    return normalised_status in {
        "solved",
        "calculated",
        "success",
        "successful",
        "completed",
    }


# =============================================================================
# Permeability calculation
# =============================================================================

def run_single_configuration(
    volume: np.ndarray,
    boundary_width: int,
    sigma: float,
    r_max: int,
) -> dict[str, Any]:
    """
    Run one X-direction permeability calculation.
    """
    start_time = time.perf_counter()

    try:
        analysis_result = analyse_permeability(
            volume=volume,
            phase_label=PHASE_LABEL,
            voxel_size_um=VOXEL_SIZE_UM,
            directions=(DIRECTION,),
            boundary_width=int(boundary_width),
            dynamic_viscosity_pa_s=DYNAMIC_VISCOSITY_PA_S,
            inlet_pressure_pa=INLET_PRESSURE_PA,
            outlet_pressure_pa=OUTLET_PRESSURE_PA,
            accuracy="standard",
            sigma=float(sigma),
            r_max=int(r_max),
            solver_tolerance=SOLVER_TOLERANCE,
            maximum_iterations=MAXIMUM_ITERATIONS,
        )

        elapsed_time = time.perf_counter() - start_time

        summary = analysis_result["summary"]

        if summary.empty:
            raise RuntimeError(
                "The permeability calculation returned an empty summary."
            )

        summary_row = summary.iloc[0]

        status = get_row_value(
            summary_row,
            ("Status", "status"),
            default="unknown",
        )

        permeability_m2 = get_row_value(
            summary_row,
            (
                "Intrinsic permeability (m^2)",
                "intrinsic_permeability_m2",
                "permeability_m2",
            ),
        )

        permeability_darcy = get_row_value(
            summary_row,
            (
                "Intrinsic permeability (Darcy)",
                "intrinsic_permeability_darcy",
                "permeability_darcy",
            ),
        )

        permeability_md = get_row_value(
            summary_row,
            (
                "Intrinsic permeability (mD)",
                "intrinsic_permeability_md",
                "permeability_md",
            ),
        )

        network_pores = get_row_value(
            summary_row,
            (
                "Network pores",
                "network_pores",
            ),
        )

        network_throats = get_row_value(
            summary_row,
            (
                "Network throats",
                "network_throats",
            ),
        )

        connectivity_percent = get_row_value(
            summary_row,
            (
                "Directional connectivity (%)",
                "directional_connectivity_percent",
                "directional_connectivity",
            ),
        )

        flow_rate_m3_s = get_row_value(
            summary_row,
            (
                "Flow rate (m^3/s)",
                "flow_rate_m3_s",
                "flow_rate",
            ),
        )

        selected_phase_voxels = get_row_value(
            summary_row,
            (
                "Selected phase voxels",
                "selected_phase_voxels",
            ),
        )

        through_connected_voxels = get_row_value(
            summary_row,
            (
                "Through-connected voxels",
                "through_connected_voxels",
            ),
        )

        if np.isfinite(float(permeability_m2)):
            permeability_m2 = abs(float(permeability_m2))

        if np.isfinite(float(permeability_darcy)):
            permeability_darcy = abs(float(permeability_darcy))

        if np.isfinite(float(permeability_md)):
            permeability_md = abs(float(permeability_md))

        if np.isfinite(float(flow_rate_m3_s)):
            flow_rate_m3_s = abs(float(flow_rate_m3_s))

        return {
            "status": str(status),
            "successful": status_is_successful(status),
            "boundary_width": int(boundary_width),
            "sigma": float(sigma),
            "r_max": int(r_max),
            "selected_phase_voxels": selected_phase_voxels,
            "through_connected_voxels": through_connected_voxels,
            "directional_connectivity_percent": connectivity_percent,
            "network_pores": network_pores,
            "network_throats": network_throats,
            "flow_rate_m3_s": flow_rate_m3_s,
            "intrinsic_permeability_m2": permeability_m2,
            "intrinsic_permeability_darcy": permeability_darcy,
            "intrinsic_permeability_md": permeability_md,
            "wall_clock_seconds": float(elapsed_time),
            "error_message": "",
        }

    except Exception as error:
        elapsed_time = time.perf_counter() - start_time

        return {
            "status": "failed",
            "successful": False,
            "boundary_width": int(boundary_width),
            "sigma": float(sigma),
            "r_max": int(r_max),
            "selected_phase_voxels": np.nan,
            "through_connected_voxels": np.nan,
            "directional_connectivity_percent": np.nan,
            "network_pores": np.nan,
            "network_throats": np.nan,
            "flow_rate_m3_s": np.nan,
            "intrinsic_permeability_m2": np.nan,
            "intrinsic_permeability_darcy": np.nan,
            "intrinsic_permeability_md": np.nan,
            "wall_clock_seconds": float(elapsed_time),
            "error_message": (
                f"{type(error).__name__}: {error}"
            ),
        }


# =============================================================================
# Sensitivity-study configuration
# =============================================================================

def create_study_configurations() -> list[dict[str, Any]]:
    """
    Create one-factor-at-a-time parameter configurations.

    The baseline configuration is:

        boundary_width = 3
        sigma = 0.4
        r_max = 4
    """
    configurations: list[dict[str, Any]] = []

    for value in (2, 3, 4):
        configurations.append(
            {
                "study": "boundary_width",
                "parameter_value": float(value),
                "boundary_width": int(value),
                "sigma": BASELINE_SIGMA,
                "r_max": BASELINE_R_MAX,
            }
        )

    for value in (0.2, 0.4, 0.6):
        configurations.append(
            {
                "study": "sigma",
                "parameter_value": float(value),
                "boundary_width": BASELINE_BOUNDARY_WIDTH,
                "sigma": float(value),
                "r_max": BASELINE_R_MAX,
            }
        )

    for value in (2, 4, 6):
        configurations.append(
            {
                "study": "r_max",
                "parameter_value": float(value),
                "boundary_width": BASELINE_BOUNDARY_WIDTH,
                "sigma": BASELINE_SIGMA,
                "r_max": int(value),
            }
        )

    return configurations


def run_sensitivity_analysis(
    volume: np.ndarray,
) -> pd.DataFrame:
    """
    Run all unique parameter configurations.

    Duplicate baseline calculations are cached so that the baseline is only
    calculated once.
    """
    configurations = create_study_configurations()

    calculation_cache: dict[
        tuple[int, float, int],
        dict[str, Any],
    ] = {}

    result_rows: list[dict[str, Any]] = []

    for index, configuration in enumerate(
        configurations,
        start=1,
    ):
        cache_key = (
            int(configuration["boundary_width"]),
            float(configuration["sigma"]),
            int(configuration["r_max"]),
        )

        print()
        print(
            f"Study {index}/{len(configurations)}: "
            f"{configuration['study']} = "
            f"{configuration['parameter_value']}"
        )
        print(
            "Configuration: "
            f"boundary_width={cache_key[0]}, "
            f"sigma={cache_key[1]}, "
            f"r_max={cache_key[2]}"
        )

        if cache_key not in calculation_cache:
            calculation_cache[cache_key] = (
                run_single_configuration(
                    volume=volume,
                    boundary_width=cache_key[0],
                    sigma=cache_key[1],
                    r_max=cache_key[2],
                )
            )
        else:
            print(
                "Using cached baseline result for this configuration."
            )

        result_row = {
            "study": configuration["study"],
            "parameter_value": configuration["parameter_value"],
            **calculation_cache[cache_key],
        }

        result_rows.append(result_row)

        print(
            f"Status: {result_row['status']}"
        )

        if result_row["successful"]:
            print(
                "Intrinsic permeability: "
                f"{result_row['intrinsic_permeability_m2']:.6e} m^2"
            )
            print(
                "Network size: "
                f"{result_row['network_pores']} pores, "
                f"{result_row['network_throats']} throats"
            )
        else:
            print(
                f"Error: {result_row['error_message']}"
            )

    results = pd.DataFrame(result_rows)

    baseline_mask = (
        (results["boundary_width"] == BASELINE_BOUNDARY_WIDTH)
        & np.isclose(results["sigma"], BASELINE_SIGMA)
        & (results["r_max"] == BASELINE_R_MAX)
        & results["successful"]
    )

    if not baseline_mask.any():
        raise RuntimeError(
            "The baseline permeability configuration did not complete "
            "successfully."
        )

    baseline_permeability = float(
        results.loc[
            baseline_mask,
            "intrinsic_permeability_m2",
        ].iloc[0]
    )

    results["baseline_permeability_m2"] = baseline_permeability

    results["relative_change_percent"] = (
        (
            results["intrinsic_permeability_m2"]
            - baseline_permeability
        )
        / baseline_permeability
        * 100.0
    )

    results["absolute_relative_change_percent"] = (
        results["relative_change_percent"].abs()
    )

    return results


# =============================================================================
# Figure
# =============================================================================

def create_sensitivity_figure(
    results: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Create a publication-style permeability sensitivity figure.
    """
    successful_results = results.loc[
        results["successful"]
    ].copy()

    if successful_results.empty:
        raise RuntimeError(
            "No successful permeability results are available for plotting."
        )

    study_settings = (
        ("boundary_width", "Boundary width (voxels)"),
        ("sigma", r"SNOW2 smoothing parameter, $\sigma$"),
        ("r_max", r"SNOW2 maximum radius, $r_{\max}$"),
    )

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(12.5, 9.0),
        dpi=150,
    )

    colours = {
        "boundary_width": "#4C78A8",
        "sigma": "#F58518",
        "r_max": "#54A24B",
    }

    for axis, (study_name, x_label) in zip(
        axes.flat[:3],
        study_settings,
    ):
        subset = successful_results.loc[
            successful_results["study"] == study_name
        ].sort_values("parameter_value")

        axis.plot(
            subset["parameter_value"],
            subset["intrinsic_permeability_m2"],
            marker="o",
            markersize=7,
            linewidth=2.0,
            color=colours[study_name],
        )

        axis.set_xlabel(x_label)
        axis.set_ylabel(
            r"Intrinsic permeability, $k$ (m$^2$)"
        )

        axis.ticklabel_format(
            axis="y",
            style="scientific",
            scilimits=(0, 0),
        )

        axis.grid(
            True,
            linestyle="--",
            linewidth=0.6,
            alpha=0.45,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

        for _, row in subset.iterrows():
            axis.annotate(
                f"{row['relative_change_percent']:+.2f}%",
                (
                    row["parameter_value"],
                    row["intrinsic_permeability_m2"],
                ),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=8,
            )

    network_axis = axes[1, 1]

    plot_index = 0

    for study_name, x_label in study_settings:
        subset = successful_results.loc[
            successful_results["study"] == study_name
        ].sort_values("parameter_value")

        x_positions = (
            np.arange(len(subset))
            + plot_index * (len(subset) + 1)
        )

        network_axis.plot(
            x_positions,
            subset["network_pores"],
            marker="o",
            linewidth=1.8,
            color=colours[study_name],
            label=f"{study_name}: pores",
        )

        network_axis.plot(
            x_positions,
            subset["network_throats"],
            marker="s",
            linestyle="--",
            linewidth=1.5,
            color=colours[study_name],
            alpha=0.75,
            label=f"{study_name}: throats",
        )

        plot_index += 1

    network_axis.set_xlabel("Sensitivity configurations")
    network_axis.set_ylabel("Network element count")
    network_axis.set_title("Extracted Network Stability")

    network_axis.grid(
        True,
        linestyle="--",
        linewidth=0.6,
        alpha=0.45,
    )

    network_axis.spines["top"].set_visible(False)
    network_axis.spines["right"].set_visible(False)

    network_axis.legend(
        fontsize=8,
        frameon=False,
        ncol=2,
    )

    figure.suptitle(
        "CemCT Permeability Parameter Sensitivity Analysis",
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

def write_summary_report(
    results: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a compact Markdown summary."""
    successful_results = results.loc[
        results["successful"]
    ].copy()

    failed_results = results.loc[
        ~results["successful"]
    ].copy()

    baseline_permeability = float(
        results["baseline_permeability_m2"].iloc[0]
    )

    lines = [
        "# CemCT Permeability Parameter Sensitivity Analysis",
        "",
        "## Baseline configuration",
        "",
        f"- Direction: {DIRECTION}",
        f"- Phase label: {PHASE_LABEL}",
        f"- Voxel size: {VOXEL_SIZE_UM} um",
        f"- Boundary width: {BASELINE_BOUNDARY_WIDTH} voxels",
        f"- Sigma: {BASELINE_SIGMA}",
        f"- r_max: {BASELINE_R_MAX}",
        (
            "- Baseline intrinsic permeability: "
            f"{baseline_permeability:.8e} m^2"
        ),
        "",
        "## Completion status",
        "",
        (
            f"- Successful study rows: "
            f"{len(successful_results)}/{len(results)}"
        ),
        f"- Failed study rows: {len(failed_results)}",
        "",
        "## Maximum sensitivity by parameter",
        "",
    ]

    for study_name in (
        "boundary_width",
        "sigma",
        "r_max",
    ):
        subset = successful_results.loc[
            successful_results["study"] == study_name
        ]

        if subset.empty:
            lines.append(
                f"- {study_name}: no successful result"
            )
            continue

        maximum_change = float(
            subset["absolute_relative_change_percent"].max()
        )

        lines.append(
            f"- {study_name}: "
            f"{maximum_change:.4f}% maximum absolute change "
            "from the baseline permeability"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "This study uses a one-factor-at-a-time design. "
                "Differences in permeability can arise from changes in "
                "boundary-network identification, image smoothing, and "
                "maximum feature size during SNOW2 network extraction."
            ),
            "",
            (
                "The selected baseline parameters should be retained for "
                "routine CemCT analysis unless an independent calibration "
                "or resolution-sensitivity study supports alternative "
                "values."
            ),
            "",
            "## Failed configurations",
            "",
        ]
    )

    if failed_results.empty:
        lines.append("- None")
    else:
        for _, row in failed_results.iterrows():
            lines.append(
                f"- {row['study']}={row['parameter_value']}: "
                f"{row['error_message']}"
            )

    output_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    """Run the complete permeability sensitivity study."""
    print(
        "Creating deterministic permeability-sensitivity domain..."
    )

    volume = create_deterministic_validation_volume()

    pore_voxel_count = int(
        np.count_nonzero(volume == PHASE_LABEL)
    )

    pore_fraction = (
        pore_voxel_count
        / volume.size
    )

    print(f"Volume shape: {volume.shape}")
    print(f"Selected phase label: {PHASE_LABEL}")
    print(f"Pore voxels: {pore_voxel_count}")
    print(f"Pore fraction: {pore_fraction:.6f}")

    tifffile.imwrite(
        VALIDATION_VOLUME_PATH,
        volume,
        photometric="minisblack",
    )

    print()
    print(
        "Running permeability parameter sensitivity analysis..."
    )
    print(
        "SNOW2 and OpenPNM calculations may take several minutes."
    )

    results = run_sensitivity_analysis(
        volume=volume,
    )

    results.to_csv(
        RESULTS_CSV_PATH,
        index=False,
    )

    create_sensitivity_figure(
        results=results,
        output_path=FIGURE_PATH,
    )

    write_summary_report(
        results=results,
        output_path=SUMMARY_PATH,
    )

    print()
    print("Permeability sensitivity analysis completed.")
    print(f"CSV: {RESULTS_CSV_PATH}")
    print(f"Figure: {FIGURE_PATH}")
    print(f"Summary: {SUMMARY_PATH}")
    print(f"Validation TIFF: {VALIDATION_VOLUME_PATH}")

    display_columns = [
        "study",
        "parameter_value",
        "boundary_width",
        "sigma",
        "r_max",
        "status",
        "network_pores",
        "network_throats",
        "intrinsic_permeability_m2",
        "relative_change_percent",
        "wall_clock_seconds",
    ]

    print()
    print("Results:")
    print(
        results[display_columns].to_string(
            index=False,
        )
    )

    failed_results = results.loc[
        ~results["successful"]
    ]

    if not failed_results.empty:
        print()
        print("WARNING: Some configurations failed:")
        print(
            failed_results[
                [
                    "study",
                    "parameter_value",
                    "error_message",
                ]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()