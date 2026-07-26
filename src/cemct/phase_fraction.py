"""User-defined multi-phase fraction analysis."""

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


VoxelSize = float | tuple[float, float, float]


def normalise_voxel_size(
    voxel_size_um: VoxelSize,
) -> tuple[float, float, float]:
    """
    Convert voxel-size input to (Z, Y, X) order.
    """
    if np.isscalar(voxel_size_um):
        value = float(voxel_size_um)

        if value <= 0:
            raise ValueError("Voxel size must be greater than zero.")

        return value, value, value

    values = tuple(float(value) for value in voxel_size_um)

    if len(values) != 3:
        raise ValueError(
            "Voxel size must be a scalar or a three-value (Z, Y, X) tuple."
        )

    if any(value <= 0 for value in values):
        raise ValueError("All voxel dimensions must be greater than zero.")

    return values


def calculate_phase_fractions(
    volume: np.ndarray,
    phases: Mapping[int, str],
    voxel_size_um: VoxelSize,
) -> pd.DataFrame:
    """
    Calculate user-selected phase fractions and physical volumes.

    Parameters
    ----------
    volume
        Integer-labelled 3D array in (Z, Y, X) order.
    phases
        Mapping between label values and phase names.

        Example:

        {
            1: "Unreacted particle",
            2: "Reaction product",
            3: "Pore",
            4: "Aggregate",
        }

    voxel_size_um
        Isotropic voxel size as a scalar, or anisotropic voxel dimensions
        in (Z, Y, X) order.

    Returns
    -------
    pandas.DataFrame
        Phase-fraction analysis results.

    Notes
    -----
    Only user-selected labels are included in the fraction denominator.
    """
    volume = np.asarray(volume)

    if volume.ndim != 3:
        raise ValueError("The input volume must be three-dimensional.")

    if not phases:
        raise ValueError("At least one phase must be selected.")

    labels = [int(label) for label in phases]

    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate phase labels are not permitted.")

    observed_labels = set(np.unique(volume).tolist())
    missing_labels = set(labels) - observed_labels

    if missing_labels:
        raise ValueError(
            f"Selected labels are absent from the volume: "
            f"{sorted(missing_labels)}"
        )

    dz, dy, dx = normalise_voxel_size(voxel_size_um)
    voxel_volume_um3 = dz * dy * dx

    analysed_mask = np.isin(volume, labels)
    analysed_voxel_count = int(np.count_nonzero(analysed_mask))

    if analysed_voxel_count == 0:
        raise ValueError("The selected phases contain no voxels.")

    rows = []

    for label, phase_name in phases.items():
        voxel_count = int(np.count_nonzero(volume == label))
        volume_fraction = voxel_count / analysed_voxel_count
        physical_volume_um3 = voxel_count * voxel_volume_um3

        rows.append(
            {
                "label": int(label),
                "phase_name": str(phase_name),
                "voxel_count": voxel_count,
                "volume_fraction": volume_fraction,
                "volume_fraction_percent": volume_fraction * 100.0,
                "physical_volume_um3": physical_volume_um3,
                "physical_volume_mm3": physical_volume_um3 / 1.0e9,
            }
        )

    return pd.DataFrame(rows)