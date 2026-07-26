"""
Random-walk transport analysis for labelled three-dimensional CT data.

This module uses PyTrax to estimate mean and directional tortuosity from the
mean-square displacement of random walkers travelling through a selected
phase.

Input convention
----------------
Any non-zero mask value is interpreted as accessible transport space.
Therefore, Boolean, 0/1, and 0/255 binary masks are accepted.

The NumPy array is interpreted in (Z, Y, X) order:

    Array axis 0 = physical Z direction
    Array axis 1 = physical Y direction
    Array axis 2 = physical X direction

Connectivity
------------
The high-level ``calculate_connected_random_walk`` function first extracts
through-connected regions using strict three-dimensional 6-connectivity.

PyTrax applies reflected periodic image boundaries during random walking.
The resulting tortuosity is obtained from the inverse slope of scaled
mean-square displacement against simulation time.
"""

from typing import Any

import numpy as np
import pandas as pd
import pytrax as pt

from cemct.connectivity import analyse_connectivity


ARRAY_AXIS_TO_DIRECTION = {
    0: "Z",
    1: "Y",
    2: "X",
}


def _validate_random_walk_parameters(
    number_of_steps: int,
    number_of_walkers: int,
    stride: int,
    fit_start_fraction: float,
    fit_end_fraction: float,
) -> tuple[int, int, int, float, float]:
    """Validate and normalise random-walk simulation parameters."""

    number_of_steps = int(number_of_steps)
    number_of_walkers = int(number_of_walkers)
    stride = int(stride)
    fit_start_fraction = float(fit_start_fraction)
    fit_end_fraction = float(fit_end_fraction)

    if number_of_steps < 2:
        raise ValueError(
            "Number of random-walk steps must be at least two."
        )

    if number_of_walkers < 1:
        raise ValueError(
            "Number of random walkers must be at least one."
        )

    if stride < 1:
        raise ValueError("Stride must be at least one.")

    if stride > number_of_steps:
        raise ValueError(
            "Stride cannot exceed the number of steps."
        )

    if number_of_steps % stride != 0:
        raise ValueError(
            "Number of steps must be exactly divisible by stride."
        )

    if not 0.0 <= fit_start_fraction < 1.0:
        raise ValueError(
            "fit_start_fraction must be greater than or equal to zero "
            "and less than one."
        )

    if not 0.0 < fit_end_fraction <= 1.0:
        raise ValueError(
            "fit_end_fraction must be greater than zero and no greater "
            "than one."
        )

    if fit_end_fraction <= fit_start_fraction:
        raise ValueError(
            "fit_end_fraction must be greater than fit_start_fraction."
        )

    return (
        number_of_steps,
        number_of_walkers,
        stride,
        fit_start_fraction,
        fit_end_fraction,
    )


def _fit_tortuosity(
    time_steps: np.ndarray,
    square_displacement: np.ndarray,
    fit_start_fraction: float,
    fit_end_fraction: float,
) -> dict[str, float]:
    """
    Fit square displacement through the origin and calculate tortuosity.

    Tortuosity is calculated as:

        tau = 1 / slope
    """

    time_steps = np.asarray(time_steps, dtype=float)
    square_displacement = np.asarray(
        square_displacement,
        dtype=float,
    )

    if time_steps.ndim != 1:
        raise ValueError("Time-step data must be one-dimensional.")

    if square_displacement.shape != time_steps.shape:
        raise ValueError(
            "Square-displacement data must match the time-step data."
        )

    number_of_points = int(time_steps.size)

    start_index = int(
        np.floor(number_of_points * fit_start_fraction)
    )

    end_index = int(
        np.ceil(number_of_points * fit_end_fraction)
    )

    start_index = max(0, start_index)
    end_index = min(number_of_points, end_index)

    fit_time = time_steps[start_index:end_index]
    fit_displacement = square_displacement[
        start_index:end_index
    ]

    if fit_time.size < 2:
        raise ValueError(
            "At least two saved time points are required for fitting."
        )

    denominator = float(np.dot(fit_time, fit_time))

    if denominator <= 0:
        raise RuntimeError(
            "The selected fitting interval contains no positive time."
        )

    slope = float(
        np.dot(fit_time, fit_displacement) / denominator
    )

    if not np.isfinite(slope) or slope <= 0:
        raise RuntimeError(
            "Random-walk MSD fitting returned a non-positive or "
            "non-finite slope."
        )

    fitted_displacement = slope * fit_time
    residual_sum_squares = float(
        np.sum(
            (fit_displacement - fitted_displacement) ** 2
        )
    )

    total_sum_squares = float(
        np.sum(
            (fit_displacement - np.mean(fit_displacement)) ** 2
        )
    )

    if total_sum_squares > 0:
        coefficient_of_determination = (
            1.0 - residual_sum_squares / total_sum_squares
        )
    else:
        coefficient_of_determination = np.nan

    tortuosity = 1.0 / slope

    return {
        "slope": slope,
        "tortuosity": tortuosity,
        "r_squared": coefficient_of_determination,
        "fit_start_index": start_index,
        "fit_end_index": end_index,
    }


def run_random_walk(
    mask: np.ndarray,
    number_of_steps: int = 10000,
    number_of_walkers: int = 5000,
    stride: int = 10,
    same_start: bool = False,
    number_of_processes: int | None = 1,
    deterministic_seed: bool = False,
    fit_start_fraction: float = 0.0,
    fit_end_fraction: float = 1.0,
    return_coordinates: bool = False,
) -> dict[str, Any]:
    """
    Run PyTrax random-walk tortuosity analysis on a 3D transport mask.

    Parameters
    ----------
    mask
        Three-dimensional Boolean, 0/1, or 0/255 transport mask.
        Non-zero values represent accessible transport space.
    number_of_steps
        Number of random-walk time steps.
    number_of_walkers
        Number of random walkers.
    stride
        Interval between saved walker coordinates.
    same_start
        If True, all walkers start at the same randomly selected location.
    number_of_processes
        Number of PyTrax processes. The default of one is safest on Windows
        and inside Jupyter.
    deterministic_seed
        Enable PyTrax's deterministic debugging seed.
    fit_start_fraction
        Fractional start of the MSD fitting interval.
    fit_end_fraction
        Fractional end of the MSD fitting interval.
    return_coordinates
        If True, include saved walker coordinates in the returned results.

    Returns
    -------
    dict
        Simulation parameters, tortuosity summary, MSD table, and optionally
        saved walker coordinates.
    """

    (
        number_of_steps,
        number_of_walkers,
        stride,
        fit_start_fraction,
        fit_end_fraction,
    ) = _validate_random_walk_parameters(
        number_of_steps=number_of_steps,
        number_of_walkers=number_of_walkers,
        stride=stride,
        fit_start_fraction=fit_start_fraction,
        fit_end_fraction=fit_end_fraction,
    )

    mask = np.asarray(mask)

    if mask.ndim != 3:
        raise ValueError(
            f"The random-walk mask must be three-dimensional; received "
            f"shape {mask.shape}."
        )

    transport_mask = mask != 0
    accessible_voxel_count = int(
        np.count_nonzero(transport_mask)
    )

    if accessible_voxel_count == 0:
        raise ValueError(
            "The random-walk mask contains no accessible voxels."
        )

    pytrax_image = transport_mask.astype(np.uint8)

    random_walk = pt.RandomWalk(
        pytrax_image,
        seed=bool(deterministic_seed),
    )

    random_walk.run(
        nt=number_of_steps,
        nw=number_of_walkers,
        same_start=bool(same_start),
        stride=stride,
        num_proc=number_of_processes,
    )

    random_walk.calc_msd()

    mean_square_displacement = np.asarray(
        random_walk.msd,
        dtype=float,
    )

    axial_mean_square_displacement = np.asarray(
        random_walk.axial_msd,
        dtype=float,
    )

    if axial_mean_square_displacement.ndim != 2:
        raise RuntimeError(
            "PyTrax returned invalid axial MSD data."
        )

    if axial_mean_square_displacement.shape[1] != 3:
        raise RuntimeError(
            "PyTrax did not return three-dimensional axial MSD data."
        )

    saved_point_count = int(mean_square_displacement.size)

    time_steps = (
        np.arange(saved_point_count, dtype=float) * stride
    )

    if axial_mean_square_displacement.shape[0] != saved_point_count:
        raise RuntimeError(
            "PyTrax mean and axial MSD lengths do not match."
        )

    mean_fit = _fit_tortuosity(
        time_steps=time_steps,
        square_displacement=mean_square_displacement,
        fit_start_fraction=fit_start_fraction,
        fit_end_fraction=fit_end_fraction,
    )

    directional_fits: dict[str, dict[str, float]] = {}
    scaled_axial_msd: dict[str, np.ndarray] = {}

    number_of_dimensions = 3

    for array_axis, physical_direction in (
        ARRAY_AXIS_TO_DIRECTION.items()
    ):
        scaled_displacement = (
            axial_mean_square_displacement[:, array_axis]
            * number_of_dimensions
        )

        scaled_axial_msd[physical_direction] = (
            scaled_displacement
        )

        directional_fits[physical_direction] = (
            _fit_tortuosity(
                time_steps=time_steps,
                square_displacement=scaled_displacement,
                fit_start_fraction=fit_start_fraction,
                fit_end_fraction=fit_end_fraction,
            )
        )

    msd_table = pd.DataFrame(
        {
            "time_step": time_steps,
            "mean_square_displacement": (
                mean_square_displacement
            ),
            "scaled_axial_msd_X": scaled_axial_msd["X"],
            "scaled_axial_msd_Y": scaled_axial_msd["Y"],
            "scaled_axial_msd_Z": scaled_axial_msd["Z"],
        }
    )

    summary = {
        "accessible_voxel_count": accessible_voxel_count,
        "accessible_fraction": (
            accessible_voxel_count / transport_mask.size
        ),
        "number_of_steps": number_of_steps,
        "number_of_walkers": number_of_walkers,
        "stride": stride,
        "saved_time_point_count": saved_point_count,
        "same_start": bool(same_start),
        "number_of_processes": number_of_processes,
        "deterministic_seed": bool(deterministic_seed),
        "fit_start_fraction": fit_start_fraction,
        "fit_end_fraction": fit_end_fraction,
        "mean_tortuosity": mean_fit["tortuosity"],
        "mean_msd_slope": mean_fit["slope"],
        "mean_r_squared": mean_fit["r_squared"],
        "X_tortuosity": directional_fits["X"][
            "tortuosity"
        ],
        "Y_tortuosity": directional_fits["Y"][
            "tortuosity"
        ],
        "Z_tortuosity": directional_fits["Z"][
            "tortuosity"
        ],
        "X_msd_slope": directional_fits["X"]["slope"],
        "Y_msd_slope": directional_fits["Y"]["slope"],
        "Z_msd_slope": directional_fits["Z"]["slope"],
        "X_r_squared": directional_fits["X"]["r_squared"],
        "Y_r_squared": directional_fits["Y"]["r_squared"],
        "Z_r_squared": directional_fits["Z"]["r_squared"],
    }

    output = {
        "summary": summary,
        "transport_mask": transport_mask,
        "msd": msd_table,
        "fit_results": {
            "mean": mean_fit,
            "X": directional_fits["X"],
            "Y": directional_fits["Y"],
            "Z": directional_fits["Z"],
        },
    }

    if return_coordinates:
        output["walker_coordinates"] = np.asarray(
            random_walk.real_coords
        )

    return output


def calculate_connected_random_walk(
    volume: np.ndarray,
    phase_label: int,
    connectivity_mode: str = "any",
    number_of_steps: int = 10000,
    number_of_walkers: int = 5000,
    stride: int = 10,
    same_start: bool = False,
    number_of_processes: int | None = 1,
    deterministic_seed: bool = False,
    fit_start_fraction: float = 0.0,
    fit_end_fraction: float = 1.0,
    return_coordinates: bool = False,
) -> dict[str, Any]:
    """
    Run random-walk analysis on through-connected regions of a labelled volume.

    Parameters
    ----------
    volume
        Integer-labelled 3D array in (Z, Y, X) order.
    phase_label
        Integer label representing the transport phase.
    connectivity_mode
        Connected mask to analyse: ``"any"``, ``"X"``, ``"Y"``, or ``"Z"``.
        ``"any"`` is the union of X-, Y-, and Z-through-connected regions.
    number_of_steps
        Number of random-walk steps.
    number_of_walkers
        Number of random walkers.
    stride
        Interval between stored coordinates.
    same_start
        Whether all walkers use the same starting point.
    number_of_processes
        Number of processes. One is recommended on Windows and Jupyter.
    deterministic_seed
        Enable PyTrax's deterministic debugging seed.
    fit_start_fraction
        Fractional beginning of the MSD fit.
    fit_end_fraction
        Fractional end of the MSD fit.
    return_coordinates
        Whether to return the full walker-coordinate array.

    Returns
    -------
    dict
        Connectivity information and random-walk analysis results.
    """

    connectivity_mode = str(
        connectivity_mode
    ).strip().upper()

    valid_modes = {
        "ANY": "through_connected_any_direction",
        "X": "through_connected_X",
        "Y": "through_connected_Y",
        "Z": "through_connected_Z",
    }

    if connectivity_mode not in valid_modes:
        raise ValueError(
            "connectivity_mode must be 'any', 'X', 'Y', or 'Z'."
        )

    connectivity_result = analyse_connectivity(
        volume=volume,
        phase_label=phase_label,
    )

    connected_mask = connectivity_result["masks"][
        valid_modes[connectivity_mode]
    ]

    connected_voxel_count = int(
        np.count_nonzero(connected_mask)
    )

    if connected_voxel_count == 0:
        raise ValueError(
            f"Selected label {phase_label} has no "
            f"{connectivity_mode}-mode through-connected voxels."
        )

    random_walk_result = run_random_walk(
        mask=connected_mask,
        number_of_steps=number_of_steps,
        number_of_walkers=number_of_walkers,
        stride=stride,
        same_start=same_start,
        number_of_processes=number_of_processes,
        deterministic_seed=deterministic_seed,
        fit_start_fraction=fit_start_fraction,
        fit_end_fraction=fit_end_fraction,
        return_coordinates=return_coordinates,
    )

    random_walk_result["summary"].update(
        {
            "phase_label": int(phase_label),
            "connectivity_mode": connectivity_mode,
            "connected_voxel_count": connected_voxel_count,
            "X_connectivity": connectivity_result["summary"][
                "X_connectivity"
            ],
            "Y_connectivity": connectivity_result["summary"][
                "Y_connectivity"
            ],
            "Z_connectivity": connectivity_result["summary"][
                "Z_connectivity"
            ],
        }
    )

    random_walk_result["connectivity_summary"] = (
        connectivity_result["summary"]
    )

    return random_walk_result