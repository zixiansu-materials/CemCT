"""
Local-thickness-based size analysis for labelled 3D CT data.

This module analyses one user-selected label from an integer-labelled
three-dimensional TIFF volume.

Method
------
PoreSpy local_thickness assigns each selected-phase voxel the radius of the
largest sphere that overlaps that voxel while remaining inside the selected
phase.

The local physical diameter is calculated as:

    local diameter = 2 x local radius x isotropic voxel size

The resulting distribution is volume weighted because every selected-phase
voxel contributes one local-diameter value.

Coordinate convention
---------------------
The NumPy array is interpreted in (Z, Y, X) order.

The current local-thickness implementation requires isotropic voxels.
"""

from typing import Any

import numpy as np
import pandas as pd
import porespy as ps


def calculate_local_diameter_map(
    volume: np.ndarray,
    phase_label: int,
    voxel_size_um: float,
    method: str = "dt",
    sizes: int = 25,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the local-diameter map for one selected label.

    Parameters
    ----------
    volume
        Integer-labelled three-dimensional array in (Z, Y, X) order.
    phase_label
        Integer label representing the phase to analyse.
    voxel_size_um
        Isotropic physical voxel size in micrometres.
    method
        PoreSpy local-thickness method. The default is ``"dt"``.
    sizes
        Number of sphere radii used by the PoreSpy calculation.

    Returns
    -------
    selected_mask
        Boolean mask of the selected label.
    local_diameter_um
        Float32 local-diameter map in micrometres. Voxels outside the selected
        phase are assigned zero.
    """

    volume = np.asarray(volume)

    if volume.ndim != 3:
        raise ValueError(
            f"The input volume must be three-dimensional; received "
            f"shape {volume.shape}."
        )

    if not np.issubdtype(volume.dtype, np.integer):
        raise TypeError(
            f"The input volume must contain integer labels; received "
            f"dtype {volume.dtype}."
        )

    phase_label = int(phase_label)
    voxel_size_um = float(voxel_size_um)
    sizes = int(sizes)

    if voxel_size_um <= 0:
        raise ValueError("Voxel size must be greater than zero.")

    if sizes < 2:
        raise ValueError("The number of local-thickness sizes must be at least 2.")

    selected_mask = volume == phase_label

    if not np.any(selected_mask):
        raise ValueError(
            f"Selected label {phase_label} is absent from the input volume."
        )

    local_radius_voxels = ps.filters.local_thickness(
        im=selected_mask,
        method=method,
        sizes=sizes,
    )

    local_diameter_um = (
        2.0
        * np.asarray(local_radius_voxels, dtype=np.float32)
        * voxel_size_um
    )

    local_diameter_um[~selected_mask] = 0.0

    return selected_mask, local_diameter_um


def _create_logarithmic_bin_edges(
    diameter_values_um: np.ndarray,
    number_of_bins: int,
) -> np.ndarray:
    """Create logarithmically spaced pore-diameter bin edges."""

    minimum_diameter = float(np.min(diameter_values_um))
    maximum_diameter = float(np.max(diameter_values_um))

    if minimum_diameter <= 0:
        raise ValueError("All analysed diameter values must be positive.")

    if np.isclose(minimum_diameter, maximum_diameter):
        lower_limit = minimum_diameter / 1.5
        upper_limit = maximum_diameter * 1.5
    else:
        lower_limit = minimum_diameter
        upper_limit = maximum_diameter

    return np.geomspace(
        lower_limit,
        upper_limit,
        number_of_bins + 1,
    )


def _calculate_size_classes(
    diameter_values_um: np.ndarray,
    total_volume_voxels: int,
    voxel_volume_um3: float,
) -> pd.DataFrame:
    """Calculate porosity contributions for three diameter classes."""

    class_definitions = [
        {
            "class_name": "Small capillary pores",
            "diameter_range_um": "<1",
            "mask": diameter_values_um < 1.0,
        },
        {
            "class_name": "Large capillary pores",
            "diameter_range_um": "1-10",
            "mask": (
                (diameter_values_um >= 1.0)
                & (diameter_values_um <= 10.0)
            ),
        },
        {
            "class_name": "Macropores",
            "diameter_range_um": ">10",
            "mask": diameter_values_um > 10.0,
        },
    ]

    selected_voxel_count = int(diameter_values_um.size)
    rows = []

    for definition in class_definitions:
        voxel_count = int(np.count_nonzero(definition["mask"]))

        rows.append(
            {
                "class_name": definition["class_name"],
                "diameter_range_um": definition["diameter_range_um"],
                "voxel_count": voxel_count,
                "fraction_of_selected_phase": (
                    voxel_count / selected_voxel_count
                ),
                "fraction_of_selected_phase_percent": (
                    100.0 * voxel_count / selected_voxel_count
                ),
                "porosity_contribution_percent": (
                    100.0 * voxel_count / total_volume_voxels
                ),
                "physical_volume_um3": (
                    voxel_count * voxel_volume_um3
                ),
                "physical_volume_mm3": (
                    voxel_count * voxel_volume_um3 / 1.0e9
                ),
            }
        )

    return pd.DataFrame(rows)


def calculate_pore_size_distribution(
    volume: np.ndarray,
    phase_label: int,
    voxel_size_um: float,
    number_of_bins: int = 30,
    method: str = "dt",
    sizes: int = 25,
) -> dict[str, Any]:
    """
    Calculate local diameter, pore-size distribution, and cumulative volume.

    Parameters
    ----------
    volume
        Integer-labelled three-dimensional array in (Z, Y, X) order.
    phase_label
        Integer label representing the phase to analyse.
    voxel_size_um
        Isotropic voxel size in micrometres.
    number_of_bins
        Number of logarithmic diameter bins.
    method
        PoreSpy local-thickness calculation method.
    sizes
        Number of sphere radii used by PoreSpy.

    Returns
    -------
    dict
        Dictionary containing the summary, local-diameter map, class table,
        distribution table, and cumulative-volume table.
    """

    number_of_bins = int(number_of_bins)

    if number_of_bins < 2:
        raise ValueError("The number of histogram bins must be at least 2.")

    selected_mask, local_diameter_um = calculate_local_diameter_map(
        volume=volume,
        phase_label=phase_label,
        voxel_size_um=voxel_size_um,
        method=method,
        sizes=sizes,
    )

    diameter_values_um = local_diameter_um[selected_mask]
    diameter_values_um = diameter_values_um[diameter_values_um > 0]

    if diameter_values_um.size == 0:
        raise RuntimeError(
            "The local-thickness calculation returned no positive diameters."
        )

    total_volume_voxels = int(volume.size)
    selected_voxel_count = int(np.count_nonzero(selected_mask))
    voxel_volume_um3 = float(voxel_size_um) ** 3

    bin_edges_um = _create_logarithmic_bin_edges(
        diameter_values_um=diameter_values_um,
        number_of_bins=number_of_bins,
    )

    histogram_counts, _ = np.histogram(
        diameter_values_um,
        bins=bin_edges_um,
    )

    bin_centres_um = np.sqrt(
        bin_edges_um[:-1] * bin_edges_um[1:]
    )

    logarithmic_bin_widths = np.diff(
        np.log10(bin_edges_um)
    )

    fraction_of_selected_phase = (
        histogram_counts / selected_voxel_count
    )

    porosity_contribution_percent = (
        100.0 * histogram_counts / total_volume_voxels
    )

    differential_porosity = (
        porosity_contribution_percent / logarithmic_bin_widths
    )

    cumulative_voxel_counts = np.cumsum(
        histogram_counts[::-1]
    )[::-1]

    cumulative_phase_fraction = (
        cumulative_voxel_counts / selected_voxel_count
    )

    cumulative_pore_volume_percent = (
        100.0 * cumulative_voxel_counts / total_volume_voxels
    )

    distribution_table = pd.DataFrame(
        {
            "bin_lower_diameter_um": bin_edges_um[:-1],
            "bin_upper_diameter_um": bin_edges_um[1:],
            "bin_centre_diameter_um": bin_centres_um,
            "log10_bin_centre": np.log10(bin_centres_um),
            "voxel_count": histogram_counts,
            "fraction_of_selected_phase": fraction_of_selected_phase,
            "porosity_contribution_percent": (
                porosity_contribution_percent
            ),
            "differential_porosity_percent_per_log10D": (
                differential_porosity
            ),
        }
    )

    cumulative_table = pd.DataFrame(
        {
            "diameter_threshold_um": bin_centres_um,
            "cumulative_voxel_count_D_greater_equal": (
                cumulative_voxel_counts
            ),
            "cumulative_fraction_of_selected_phase": (
                cumulative_phase_fraction
            ),
            "cumulative_pore_volume_percent_of_sample": (
                cumulative_pore_volume_percent
            ),
        }
    )

    class_table = _calculate_size_classes(
        diameter_values_um=diameter_values_um,
        total_volume_voxels=total_volume_voxels,
        voxel_volume_um3=voxel_volume_um3,
    )

    summary = {
        "phase_label": int(phase_label),
        "volume_shape_ZYX": tuple(int(value) for value in volume.shape),
        "voxel_size_um": float(voxel_size_um),
        "total_volume_voxel_count": total_volume_voxels,
        "selected_voxel_count": selected_voxel_count,
        "selected_phase_fraction": (
            selected_voxel_count / total_volume_voxels
        ),
        "selected_phase_fraction_percent": (
            100.0 * selected_voxel_count / total_volume_voxels
        ),
        "selected_physical_volume_um3": (
            selected_voxel_count * voxel_volume_um3
        ),
        "minimum_local_diameter_um": float(
            np.min(diameter_values_um)
        ),
        "median_local_diameter_um": float(
            np.median(diameter_values_um)
        ),
        "mean_local_diameter_um": float(
            np.mean(diameter_values_um)
        ),
        "maximum_local_diameter_um": float(
            np.max(diameter_values_um)
        ),
        "number_of_bins": number_of_bins,
        "local_thickness_method": str(method),
        "local_thickness_sizes": int(sizes),
    }

    return {
        "summary": summary,
        "selected_mask": selected_mask,
        "local_diameter_um": local_diameter_um,
        "size_classes": class_table,
        "distribution": distribution_table,
        "cumulative": cumulative_table,
    }