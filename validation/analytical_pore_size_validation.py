"""
Analytical pore-size validation for CemCT.

This script validates the local-thickness-based pore-diameter calculation
using straight square channels with known cross-sectional widths.

For a straight square channel:

    expected local diameter = channel width × voxel size

Only the central region of each channel is evaluated to minimise the influence
of the inlet and outlet boundaries.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemct.pore_size import calculate_local_diameter_map


PHASE_LABEL = 3
VOXEL_SIZE_UM = 1.0

VOLUME_SHAPE = (48, 48, 64)
CHANNEL_WIDTHS_VOXELS = (6, 10, 14)

RELATIVE_ERROR_LIMIT_PERCENT = 10.0

OUTPUT_DIRECTORY = Path("validation") / "results"
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)


def create_straight_x_channel(
    shape: tuple[int, int, int],
    channel_width: int,
    phase_label: int = PHASE_LABEL,
) -> np.ndarray:
    """
    Create a labelled straight square channel spanning the X direction.

    Array convention
    ----------------
    axis 0 = Z
    axis 1 = Y
    axis 2 = X
    """
    if channel_width <= 0:
        raise ValueError("Channel width must be greater than zero.")

    if channel_width >= min(shape[0], shape[1]):
        raise ValueError(
            "Channel width must be smaller than the Z and Y dimensions."
        )

    volume = np.zeros(shape, dtype=np.uint8)

    z_start = (shape[0] - channel_width) // 2
    y_start = (shape[1] - channel_width) // 2

    z_stop = z_start + channel_width
    y_stop = y_start + channel_width

    volume[
        z_start:z_stop,
        y_start:y_stop,
        :,
    ] = phase_label

    return volume


def extract_diameter_map(
    result: tuple[np.ndarray, np.ndarray],
) -> np.ndarray:
    """
    Extract the local-diameter array returned by CemCT.

    This helper also checks that the result contains two 3D arrays.
    """
    if not isinstance(result, tuple) or len(result) != 2:
        raise TypeError(
            "calculate_local_diameter_map must return a tuple of two arrays."
        )

    first_array = np.asarray(result[0])
    second_array = np.asarray(result[1])

    if first_array.ndim != 3 or second_array.ndim != 3:
        raise ValueError("Both returned arrays must be three-dimensional.")

    # The selected mask is normally Boolean, while the diameter map is float.
    if np.issubdtype(second_array.dtype, np.floating):
        return second_array

    if np.issubdtype(first_array.dtype, np.floating):
        return first_array

    # Fallback for unexpected data types.
    return second_array.astype(np.float64)


def calculate_channel_result(
    channel_width_voxels: int,
) -> dict[str, float | int | str]:
    """Run local-thickness analysis for one straight channel."""
    volume = create_straight_x_channel(
        shape=VOLUME_SHAPE,
        channel_width=channel_width_voxels,
    )

    calculation_result = calculate_local_diameter_map(
        volume=volume,
        phase_label=PHASE_LABEL,
        voxel_size_um=VOXEL_SIZE_UM,
        method="dt",
        sizes=50,
    )

    local_diameter_um = extract_diameter_map(calculation_result)

    x_size = volume.shape[2]
    x_start = x_size // 4
    x_stop = 3 * x_size // 4

    central_selected_mask = (
        volume[:, :, x_start:x_stop] == PHASE_LABEL
    )

    central_diameters = local_diameter_um[
        :, :, x_start:x_stop
    ][central_selected_mask]

    central_diameters = central_diameters[
        np.isfinite(central_diameters)
        & (central_diameters > 0)
    ]

    if central_diameters.size == 0:
        raise RuntimeError(
            f"No positive local diameters were calculated for "
            f"channel width {channel_width_voxels}."
        )

    expected_diameter_um = (
        channel_width_voxels * VOXEL_SIZE_UM
    )

    measured_median_um = float(
        np.median(central_diameters)
    )

    measured_mean_um = float(
        np.mean(central_diameters)
    )

    absolute_error_um = abs(
        measured_median_um - expected_diameter_um
    )

    relative_error_percent = (
        absolute_error_um
        / expected_diameter_um
        * 100.0
    )

    status = (
        "passed"
        if relative_error_percent
        <= RELATIVE_ERROR_LIMIT_PERCENT
        else "review"
    )

    return {
        "channel_width_voxels": channel_width_voxels,
        "voxel_size_um": VOXEL_SIZE_UM,
        "expected_diameter_um": expected_diameter_um,
        "measured_median_diameter_um": measured_median_um,
        "measured_mean_diameter_um": measured_mean_um,
        "absolute_error_um": absolute_error_um,
        "relative_error_percent": relative_error_percent,
        "evaluated_voxel_count": int(
            central_diameters.size
        ),
        "acceptance_limit_percent": (
            RELATIVE_ERROR_LIMIT_PERCENT
        ),
        "status": status,
    }


def create_validation_figure(
    results: pd.DataFrame,
    output_path: Path,
) -> None:
    """Create a publication-style validation figure."""
    expected = results["expected_diameter_um"].to_numpy(
        dtype=float
    )

    measured = results[
        "measured_median_diameter_um"
    ].to_numpy(dtype=float)

    maximum_value = max(
        float(expected.max()),
        float(measured.max()),
    )

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(10, 4.2),
        dpi=150,
    )

    axes[0].plot(
        expected,
        expected,
        linestyle="--",
        linewidth=1.5,
        color="black",
        label="Ideal agreement",
    )

    axes[0].scatter(
        expected,
        measured,
        s=70,
        color="#2F6B9A",
        edgecolor="black",
        linewidth=0.6,
        zorder=3,
        label="CemCT result",
    )

    axes[0].set_xlim(0, maximum_value * 1.1)
    axes[0].set_ylim(0, maximum_value * 1.1)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("Expected diameter (µm)")
    axes[0].set_ylabel("Measured median diameter (µm)")
    axes[0].set_title("Expected versus measured diameter")
    axes[0].legend(frameon=False)

    axes[1].bar(
        results["channel_width_voxels"].astype(str),
        results["relative_error_percent"],
        color="#E07A5F",
        edgecolor="black",
        linewidth=0.6,
    )

    axes[1].axhline(
        RELATIVE_ERROR_LIMIT_PERCENT,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="Acceptance limit",
    )

    axes[1].set_xlabel("Channel width (voxels)")
    axes[1].set_ylabel("Relative error (%)")
    axes[1].set_title("Local-diameter error")
    axes[1].legend(frameon=False)

    for axis in axes:
        axis.grid(
            True,
            linestyle="--",
            linewidth=0.5,
            alpha=0.35,
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        "CemCT Analytical Pore-Size Validation",
        fontsize=13,
        fontweight="bold",
    )

    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def write_summary(
    results: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a Markdown summary of the validation."""
    all_passed = bool(
        results["status"].eq("passed").all()
    )

    maximum_error = float(
        results["relative_error_percent"].max()
    )

    lines = [
        "# Analytical pore-size validation",
        "",
        "## Objective",
        "",
        "Validate the CemCT local-thickness-based pore-diameter "
        "calculation using straight square channels of known width.",
        "",
        "## Analytical reference",
        "",
        "Expected diameter = channel width × voxel size.",
        "",
        "## Acceptance criterion",
        "",
        f"Relative median-diameter error <= "
        f"{RELATIVE_ERROR_LIMIT_PERCENT:.1f}%.",
        "",
        "## Results",
        "",
        results.to_markdown(index=False),
        "",
        "## Overall assessment",
        "",
        f"- All configurations passed: `{all_passed}`",
        f"- Maximum relative error: `{maximum_error:.6f}%`",
        "",
    ]

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    """Run the complete analytical pore-size validation."""
    print("Running analytical pore-size validation...")

    result_rows = []

    for channel_width in CHANNEL_WIDTHS_VOXELS:
        print(
            f"Testing channel width: "
            f"{channel_width} voxels"
        )

        result_rows.append(
            calculate_channel_result(channel_width)
        )

    results = pd.DataFrame(result_rows)

    csv_path = (
        OUTPUT_DIRECTORY
        / "analytical_pore_size_validation_results.csv"
    )

    figure_path = (
        OUTPUT_DIRECTORY
        / "analytical_pore_size_validation.png"
    )

    summary_path = (
        OUTPUT_DIRECTORY
        / "analytical_pore_size_validation_summary.md"
    )

    results.to_csv(csv_path, index=False)

    create_validation_figure(
        results=results,
        output_path=figure_path,
    )

    write_summary(
        results=results,
        output_path=summary_path,
    )

    print()
    print(results.to_string(index=False))
    print()
    print("Validation completed.")
    print("CSV:", csv_path.resolve())
    print("Figure:", figure_path.resolve())
    print("Summary:", summary_path.resolve())

    if not results["status"].eq("passed").all():
        raise RuntimeError(
            "One or more pore-size validation cases exceeded "
            "the acceptance limit."
        )


if __name__ == "__main__":
    main()