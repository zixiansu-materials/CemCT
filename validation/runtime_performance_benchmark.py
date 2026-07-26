"""
Runtime and memory performance benchmark for CemCT.

This script benchmarks the six principal CemCT analysis modules:

    1. Multi-phase fraction analysis
    2. Pore connectivity analysis
    3. Local-thickness pore-size analysis
    4. Finite-difference diffusion analysis
    5. Random-walk transport analysis
    6. Pore-network permeability analysis

The benchmark uses deterministic synthetic labelled volumes at several
volume sizes. The same generated volume is supplied to every module for a
given size.

The transport calculations are limited to the X direction to keep the
benchmark computationally practical.

Outputs
-------
    1. CSV file containing every benchmark run
    2. CSV file containing aggregated performance results
    3. JSON file containing system and benchmark metadata
    4. PNG figure showing runtime and memory scaling
    5. Markdown performance summary

Array convention
----------------
    Array axis 0 = Z direction
    Array axis 1 = Y direction
    Array axis 2 = X direction

Synthetic phase convention
--------------------------
    Label 1 = solid phase
    Label 3 = connected pore phase
"""

from __future__ import annotations

import gc
import json
import platform
import sys
import threading
import time
import traceback
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from scipy.ndimage import gaussian_filter

from cemct.connectivity import analyse_connectivity
from cemct.diffusion import calculate_directional_diffusion
from cemct.permeability import analyse_permeability
from cemct.phase_fraction import calculate_phase_fractions
from cemct.pore_size import calculate_pore_size_distribution
from cemct.random_walk import calculate_connected_random_walk


# =============================================================================
# User-adjustable benchmark configuration
# =============================================================================

# Start with these three moderate sizes.
# After confirming successful execution, 100 can optionally be added.
VOLUME_SIZES = (40, 60, 80)

# Use one repeat for the initial benchmark.
# Change to 3 for the final manuscript benchmark.
NUMBER_OF_REPEATS = 3

RANDOM_SEED = 20260727

PHASE_LABEL = 3
VOXEL_SIZE_UM = 0.7

# Random-walk settings selected to keep the benchmark practical.
RANDOM_WALK_STEPS = 2000
RANDOM_WALK_WALKERS = 1000
RANDOM_WALK_STRIDE = 20

# Diffusion-solver settings.
SOLVER_TOLERANCE = 1.0e-10
MAXIMUM_ITERATIONS = 5000

# Permeability-network extraction settings.
BOUNDARY_WIDTH = 3
SNOW_SIGMA = 0.4
SNOW_R_MAX = 4

# Memory polling interval.
MEMORY_POLLING_INTERVAL_SECONDS = 0.02

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = PROJECT_ROOT / "validation" / "results"

OUTPUT_DIRECTORY.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# Package-version utilities
# =============================================================================

def get_package_version(package_name: str) -> str:
    """Return an installed package version or 'not installed'."""
    try:
        return version(package_name)
    except PackageNotFoundError:
        return "not installed"


# =============================================================================
# Deterministic benchmark-volume generation
# =============================================================================

def create_benchmark_volume(
    size: int,
    random_seed: int,
) -> np.ndarray:
    """
    Create a deterministic labelled porous volume.

    The stochastic porous field is supplemented with guaranteed X-, Y-, and
    Z-directional face-connected pore channels. This ensures that directional
    transport calculations receive an appropriate through-connected phase.

    Parameters
    ----------
    size
        Length of every side of the cubic volume.
    random_seed
        Seed used for deterministic volume generation.

    Returns
    -------
    numpy.ndarray
        Labelled uint8 volume in (Z, Y, X) order.
    """
    if size < 20:
        raise ValueError(
            "Benchmark volume size must be at least 20 voxels."
        )

    rng = np.random.default_rng(
        random_seed + size
    )

    random_field = rng.random(
        (size, size, size)
    )

    smoothing_sigma = max(
        1.0,
        size / 32.0,
    )

    smoothed_field = gaussian_filter(
        random_field,
        sigma=smoothing_sigma,
        mode="reflect",
    )

    # Approximately 38% stochastic pore fraction before channels are added.
    threshold = np.quantile(
        smoothed_field,
        0.62,
    )

    pore_mask = smoothed_field >= threshold

    centre = size // 2
    channel_half_width = max(
        2,
        size // 16,
    )

    lower = centre - channel_half_width
    upper = centre + channel_half_width + 1

    # Guaranteed X-directional pore channel.
    pore_mask[
        lower:upper,
        lower:upper,
        :
    ] = True

    # Guaranteed Y-directional pore channel.
    pore_mask[
        lower:upper,
        :,
        lower:upper,
    ] = True

    # Guaranteed Z-directional pore channel.
    pore_mask[
        :,
        lower:upper,
        lower:upper,
    ] = True

    labelled_volume = np.ones(
        pore_mask.shape,
        dtype=np.uint8,
    )

    labelled_volume[pore_mask] = PHASE_LABEL

    return labelled_volume


# =============================================================================
# Peak-memory monitoring
# =============================================================================

class PeakMemoryMonitor:
    """Monitor the peak resident memory of the current Python process."""

    def __init__(
        self,
        polling_interval_seconds: float = 0.02,
    ) -> None:
        self.process = psutil.Process()
        self.polling_interval_seconds = polling_interval_seconds
        self.stop_event = threading.Event()
        self.thread: threading.Thread | None = None

        self.baseline_rss_bytes = int(
            self.process.memory_info().rss
        )

        self.peak_rss_bytes = self.baseline_rss_bytes

    def _monitor(self) -> None:
        while not self.stop_event.is_set():
            current_rss = int(
                self.process.memory_info().rss
            )

            self.peak_rss_bytes = max(
                self.peak_rss_bytes,
                current_rss,
            )

            self.stop_event.wait(
                self.polling_interval_seconds
            )

    def start(self) -> None:
        self.thread = threading.Thread(
            target=self._monitor,
            daemon=True,
        )

        self.thread.start()

    def stop(self) -> None:
        self.stop_event.set()

        if self.thread is not None:
            self.thread.join(
                timeout=2.0
            )

        current_rss = int(
            self.process.memory_info().rss
        )

        self.peak_rss_bytes = max(
            self.peak_rss_bytes,
            current_rss,
        )

    @property
    def baseline_rss_mb(self) -> float:
        return self.baseline_rss_bytes / (1024.0 ** 2)

    @property
    def peak_rss_mb(self) -> float:
        return self.peak_rss_bytes / (1024.0 ** 2)

    @property
    def additional_peak_rss_mb(self) -> float:
        difference = (
            self.peak_rss_bytes
            - self.baseline_rss_bytes
        )

        return max(
            difference,
            0,
        ) / (1024.0 ** 2)


# =============================================================================
# Individual benchmark execution
# =============================================================================

def benchmark_function(
    module_name: str,
    operation: Callable[[], Any],
    volume_size: int,
    repeat_number: int,
    pore_fraction: float,
) -> dict[str, Any]:
    """Benchmark one CemCT operation."""
    gc.collect()

    memory_monitor = PeakMemoryMonitor(
        polling_interval_seconds=(
            MEMORY_POLLING_INTERVAL_SECONDS
        )
    )

    start_time = time.perf_counter()
    memory_monitor.start()

    status = "calculated"
    error_message = ""

    try:
        result = operation()

        # Keep the result alive until after peak memory has been sampled.
        _ = result

    except Exception as error:
        status = "failed"

        error_message = (
            f"{type(error).__name__}: {error}"
        )

        traceback.print_exc()

    finally:
        elapsed_seconds = (
            time.perf_counter()
            - start_time
        )

        memory_monitor.stop()

    row = {
        "module": module_name,
        "volume_size": int(volume_size),
        "volume_shape": (
            f"{volume_size} x "
            f"{volume_size} x "
            f"{volume_size}"
        ),
        "voxel_count": int(volume_size ** 3),
        "pore_fraction": float(pore_fraction),
        "repeat": int(repeat_number),
        "status": status,
        "wall_time_seconds": float(
            elapsed_seconds
        ),
        "baseline_rss_mb": float(
            memory_monitor.baseline_rss_mb
        ),
        "peak_rss_mb": float(
            memory_monitor.peak_rss_mb
        ),
        "additional_peak_rss_mb": float(
            memory_monitor.additional_peak_rss_mb
        ),
        "error": error_message,
    }

    del operation

    if "result" in locals():
        del result

    gc.collect()

    return row


# =============================================================================
# Module-operation definitions
# =============================================================================

def create_module_operations(
    volume: np.ndarray,
) -> list[tuple[str, Callable[[], Any]]]:
    """Create the six module calls for one benchmark volume."""

    operations = [
        (
            "phase_fraction",
            lambda: calculate_phase_fractions(
                volume=volume,
                phases={
                    1: "Solid phase",
                    3: "Pore phase",
                },
                voxel_size_um=VOXEL_SIZE_UM,
            ),
        ),
        (
            "connectivity",
            lambda: analyse_connectivity(
                volume=volume,
                phase_label=PHASE_LABEL,
            ),
        ),
        (
            "pore_size",
            lambda: calculate_pore_size_distribution(
                volume=volume,
                phase_label=PHASE_LABEL,
                voxel_size_um=VOXEL_SIZE_UM,
                number_of_bins=25,
                method="dt",
                sizes=25,
            ),
        ),
        (
            "finite_difference_diffusion",
            lambda: calculate_directional_diffusion(
                volume=volume,
                phase_label=PHASE_LABEL,
                directions=("X",),
                solver=None,
                solver_tolerance=SOLVER_TOLERANCE,
                maximum_iterations=MAXIMUM_ITERATIONS,
            ),
        ),
        (
            "random_walk",
            lambda: calculate_connected_random_walk(
                volume=volume,
                phase_label=PHASE_LABEL,
                connectivity_mode="X",
                number_of_steps=RANDOM_WALK_STEPS,
                number_of_walkers=RANDOM_WALK_WALKERS,
                stride=RANDOM_WALK_STRIDE,
                same_start=False,
                number_of_processes=1,
                deterministic_seed=True,
                fit_start_fraction=0.2,
                fit_end_fraction=0.8,
                return_coordinates=False,
            ),
        ),
        (
            "permeability",
            lambda: analyse_permeability(
                volume=volume,
                phase_label=PHASE_LABEL,
                voxel_size_um=VOXEL_SIZE_UM,
                directions=("X",),
                boundary_width=BOUNDARY_WIDTH,
                dynamic_viscosity_pa_s=1.0,
                inlet_pressure_pa=1.0,
                outlet_pressure_pa=0.0,
                accuracy="standard",
                sigma=SNOW_SIGMA,
                r_max=SNOW_R_MAX,
                solver_tolerance=SOLVER_TOLERANCE,
                maximum_iterations=MAXIMUM_ITERATIONS,
            ),
        ),
    ]

    return operations


# =============================================================================
# Result aggregation
# =============================================================================

def aggregate_results(
    results: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate median runtime and memory for successful runs."""
    successful = results.loc[
        results["status"] == "calculated"
    ].copy()

    if successful.empty:
        return pd.DataFrame()

    aggregated = (
        successful.groupby(
            [
                "module",
                "volume_size",
                "voxel_count",
            ],
            as_index=False,
        )
        .agg(
            median_wall_time_seconds=(
                "wall_time_seconds",
                "median",
            ),
            minimum_wall_time_seconds=(
                "wall_time_seconds",
                "min",
            ),
            maximum_wall_time_seconds=(
                "wall_time_seconds",
                "max",
            ),
            median_peak_rss_mb=(
                "peak_rss_mb",
                "median",
            ),
            median_additional_peak_rss_mb=(
                "additional_peak_rss_mb",
                "median",
            ),
            successful_runs=(
                "status",
                "count",
            ),
        )
    )

    return aggregated


# =============================================================================
# Figure generation
# =============================================================================

def create_performance_figure(
    aggregated: pd.DataFrame,
    output_path: Path,
) -> None:
    """Create runtime and memory-scaling plots."""
    if aggregated.empty:
        print(
            "No successful benchmark results were available "
            "for plotting."
        )
        return

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(14, 5.5),
    )

    for module_name, module_results in aggregated.groupby(
        "module"
    ):
        module_results = module_results.sort_values(
            "voxel_count"
        )

        axes[0].plot(
            module_results["voxel_count"],
            module_results[
                "median_wall_time_seconds"
            ],
            marker="o",
            linewidth=1.8,
            label=module_name,
        )

        memory_values = np.maximum(
            module_results[
                "median_additional_peak_rss_mb"
            ].to_numpy(dtype=float),
            0.01,
        )

        axes[1].plot(
            module_results["voxel_count"],
            memory_values,
            marker="o",
            linewidth=1.8,
            label=module_name,
        )

    axes[0].set_title(
        "Runtime Scaling"
    )

    axes[0].set_xlabel(
        "Number of voxels"
    )

    axes[0].set_ylabel(
        "Median wall time (s)"
    )

    axes[0].set_xscale(
        "log"
    )

    axes[0].set_yscale(
        "log"
    )

    axes[1].set_title(
        "Additional Peak Memory"
    )

    axes[1].set_xlabel(
        "Number of voxels"
    )

    axes[1].set_ylabel(
        "Additional peak RSS (MB)"
    )

    axes[1].set_xscale(
        "log"
    )

    axes[1].set_yscale(
        "log"
    )

    for axis in axes:
        axis.grid(
            True,
            which="both",
            linestyle="--",
            linewidth=0.6,
            alpha=0.45,
        )

        axis.spines["top"].set_visible(
            False
        )

        axis.spines["right"].set_visible(
            False
        )

    handles, labels = axes[0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )

    figure.suptitle(
        "CemCT Runtime and Memory Performance Benchmark",
        fontsize=14,
        fontweight="bold",
    )

    figure.tight_layout(
        rect=(0.0, 0.08, 1.0, 0.95)
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )


# =============================================================================
# Metadata and Markdown output
# =============================================================================

def create_system_metadata() -> dict[str, Any]:
    """Collect system and benchmark configuration metadata."""
    virtual_memory = psutil.virtual_memory()

    return {
        "benchmark_name": (
            "CemCT runtime and performance benchmark"
        ),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "operating_system": platform.platform(),
        "processor": platform.processor(),
        "physical_cpu_count": psutil.cpu_count(
            logical=False
        ),
        "logical_cpu_count": psutil.cpu_count(
            logical=True
        ),
        "total_system_memory_gb": (
            virtual_memory.total
            / (1024.0 ** 3)
        ),
        "cemct_version": get_package_version(
            "cemct"
        ),
        "numpy_version": get_package_version(
            "numpy"
        ),
        "scipy_version": get_package_version(
            "scipy"
        ),
        "pandas_version": get_package_version(
            "pandas"
        ),
        "porespy_version": get_package_version(
            "porespy"
        ),
        "openpnm_version": get_package_version(
            "openpnm"
        ),
        "pyamg_version": get_package_version(
            "pyamg"
        ),
        "pytrax_version": get_package_version(
            "pytrax"
        ),
        "psutil_version": get_package_version(
            "psutil"
        ),
        "volume_sizes": list(
            VOLUME_SIZES
        ),
        "number_of_repeats": NUMBER_OF_REPEATS,
        "voxel_size_um": VOXEL_SIZE_UM,
        "phase_label": PHASE_LABEL,
        "random_seed": RANDOM_SEED,
        "random_walk_steps": RANDOM_WALK_STEPS,
        "random_walk_walkers": RANDOM_WALK_WALKERS,
        "random_walk_stride": RANDOM_WALK_STRIDE,
        "solver_tolerance": SOLVER_TOLERANCE,
        "maximum_iterations": MAXIMUM_ITERATIONS,
        "permeability_boundary_width": BOUNDARY_WIDTH,
        "permeability_sigma": SNOW_SIGMA,
        "permeability_r_max": SNOW_R_MAX,
        "transport_directions": ["X"],
    }


def write_markdown_summary(
    results: pd.DataFrame,
    aggregated: pd.DataFrame,
    metadata: dict[str, Any],
    output_path: Path,
) -> None:
    """Write a publication-oriented Markdown summary."""
    successful_count = int(
        (results["status"] == "calculated").sum()
    )

    failed_count = int(
        (results["status"] == "failed").sum()
    )

    lines = [
        "# CemCT Runtime and Performance Benchmark",
        "",
        "## Purpose",
        "",
        (
            "This benchmark evaluates the wall-clock runtime and "
            "resident-memory requirements of the six principal CemCT "
            "analysis modules."
        ),
        "",
        "## Benchmark configuration",
        "",
        f"- Volume sizes: `{list(VOLUME_SIZES)}`",
        f"- Repeats per case: `{NUMBER_OF_REPEATS}`",
        f"- Voxel size: `{VOXEL_SIZE_UM} um`",
        "- Transport direction: `X`",
        f"- Random-walk steps: `{RANDOM_WALK_STEPS}`",
        f"- Random-walk walkers: `{RANDOM_WALK_WALKERS}`",
        f"- Solver tolerance: `{SOLVER_TOLERANCE:.1e}`",
        f"- Maximum solver iterations: `{MAXIMUM_ITERATIONS}`",
        "",
        "## Computing environment",
        "",
        f"- Operating system: `{metadata['operating_system']}`",
        f"- Python: `{metadata['python_version'].split()[0]}`",
        (
            f"- Logical CPU count: "
            f"`{metadata['logical_cpu_count']}`"
        ),
        (
            f"- Total system memory: "
            f"`{metadata['total_system_memory_gb']:.2f} GB`"
        ),
        f"- CemCT version: `{metadata['cemct_version']}`",
        "",
        "## Completion status",
        "",
        f"- Successful benchmark runs: `{successful_count}`",
        f"- Failed benchmark runs: `{failed_count}`",
        "",
        "## Aggregated results",
        "",
    ]

    if aggregated.empty:
        lines.append(
            "No successful benchmark results were available."
        )
    else:
        display_columns = [
            "module",
            "volume_size",
            "voxel_count",
            "median_wall_time_seconds",
            "median_additional_peak_rss_mb",
            "successful_runs",
        ]

        display_table = aggregated[
            display_columns
        ].copy()

        lines.append(
            display_table.to_markdown(
                index=False,
                floatfmt=".4g",
            )
        )

    failed_rows = results.loc[
        results["status"] == "failed",
        [
            "module",
            "volume_size",
            "repeat",
            "error",
        ],
    ]

    lines.extend(
        [
            "",
            "## Failed calculations",
            "",
        ]
    )

    if failed_rows.empty:
        lines.append(
            "No module failures were recorded."
        )
    else:
        lines.append(
            failed_rows.to_markdown(
                index=False
            )
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "Runtime and memory values are implementation- and "
                "hardware-specific. They should therefore be interpreted "
                "as performance measurements for the computing environment "
                "reported above, rather than universal execution times."
            ),
            "",
            (
                "Finite-difference diffusion and pore-network permeability "
                "were calculated only in the X direction to keep the "
                "benchmark computationally practical."
            ),
            "",
            (
                "The reported additional peak RSS is the maximum observed "
                "increase in process resident memory relative to the "
                "baseline immediately before each module call."
            ),
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


# =============================================================================
# Main benchmark
# =============================================================================

def main() -> None:
    """Run the complete CemCT performance benchmark."""
    print(
        "Starting CemCT runtime and performance benchmark..."
    )

    print(
        f"Volume sizes: {VOLUME_SIZES}"
    )

    print(
        f"Repeats: {NUMBER_OF_REPEATS}"
    )

    all_rows: list[dict[str, Any]] = []

    for volume_size in VOLUME_SIZES:
        print()
        print(
            "=" * 72
        )

        print(
            f"Creating benchmark volume: "
            f"{volume_size} x "
            f"{volume_size} x "
            f"{volume_size}"
        )

        volume = create_benchmark_volume(
            size=volume_size,
            random_seed=RANDOM_SEED,
        )

        pore_fraction = float(
            np.mean(
                volume == PHASE_LABEL
            )
        )

        print(
            f"Pore fraction: {pore_fraction:.6f}"
        )

        for repeat_number in range(
            1,
            NUMBER_OF_REPEATS + 1,
        ):
            print(
                f"Repeat: {repeat_number}/"
                f"{NUMBER_OF_REPEATS}"
            )

            module_operations = create_module_operations(
                volume=volume
            )

            for module_name, operation in module_operations:
                print(
                    f"Running {module_name}..."
                )

                row = benchmark_function(
                    module_name=module_name,
                    operation=operation,
                    volume_size=volume_size,
                    repeat_number=repeat_number,
                    pore_fraction=pore_fraction,
                )

                all_rows.append(
                    row
                )

                print(
                    f"  status: {row['status']}"
                )

                print(
                    f"  wall time: "
                    f"{row['wall_time_seconds']:.4f} s"
                )

                print(
                    f"  additional peak RSS: "
                    f"{row['additional_peak_rss_mb']:.2f} MB"
                )

                if row["error"]:
                    print(
                        f"  error: {row['error']}"
                    )

        del volume

        gc.collect()

    results = pd.DataFrame(
        all_rows
    )

    aggregated = aggregate_results(
        results
    )

    results_path = (
        OUTPUT_DIRECTORY
        / "runtime_performance_benchmark_results.csv"
    )

    aggregated_path = (
        OUTPUT_DIRECTORY
        / "runtime_performance_benchmark_aggregated.csv"
    )

    metadata_path = (
        OUTPUT_DIRECTORY
        / "runtime_performance_benchmark_metadata.json"
    )

    figure_path = (
        OUTPUT_DIRECTORY
        / "runtime_performance_benchmark.png"
    )

    summary_path = (
        OUTPUT_DIRECTORY
        / "runtime_performance_benchmark_summary.md"
    )

    results.to_csv(
        results_path,
        index=False,
    )

    aggregated.to_csv(
        aggregated_path,
        index=False,
    )

    metadata = create_system_metadata()

    metadata_path.write_text(
        json.dumps(
            metadata,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    create_performance_figure(
        aggregated=aggregated,
        output_path=figure_path,
    )

    write_markdown_summary(
        results=results,
        aggregated=aggregated,
        metadata=metadata,
        output_path=summary_path,
    )

    print()
    print(
        "=" * 72
    )

    print(
        "Runtime and performance benchmark completed."
    )

    print(
        f"Detailed CSV: {results_path}"
    )

    print(
        f"Aggregated CSV: {aggregated_path}"
    )

    print(
        f"Metadata JSON: {metadata_path}"
    )

    print(
        f"Figure: {figure_path}"
    )

    print(
        f"Summary: {summary_path}"
    )

    print()
    print(
        "Results:"
    )

    if aggregated.empty:
        print(
            "No successful benchmark results."
        )
    else:
        print(
            aggregated.to_string(
                index=False
            )
        )


if __name__ == "__main__":
    main()