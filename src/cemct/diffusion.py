"""
Finite-difference diffusion tortuosity analysis for labelled 3D CT data.

This module calculates direction-dependent diffusion transport properties
from directionally through-connected regions of a user-selected phase.

Connectivity
------------
Directional through-connected masks are extracted using strict
three-dimensional 6-connectivity:

    Face connection   = included
    Edge connection   = excluded
    Corner connection = excluded

Direction convention
--------------------
The NumPy array is interpreted in (Z, Y, X) order:

    X direction = array axis 2
    Y direction = array axis 1
    Z direction = array axis 0

Transport definitions
---------------------
Formation factor:

    F = D0 / Deff

Relative effective diffusivity:

    Deff / D0 = 1 / F

Diffusion tortuosity:

    tau = effective porosity x F

Only a phase that connects the two opposing faces in the requested direction
can be analysed in that direction.
"""

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
import porespy as ps

from cemct.connectivity import analyse_connectivity


DIRECTION_TO_AXIS = {
    "X": 2,
    "Y": 1,
    "Z": 0,
}


def create_pyamg_solver(
    tolerance: float = 1.0e-10,
    maximum_iterations: int = 5000,
) -> Any:
    """
    Create an OpenPNM PyAMG Ruge-Stuben iterative solver.

    Parameters
    ----------
    tolerance
        Numerical convergence tolerance.
    maximum_iterations
        Maximum number of solver iterations.

    Returns
    -------
    object
        OpenPNM PyamgRugeStubenSolver instance.
    """

    tolerance = float(tolerance)
    maximum_iterations = int(maximum_iterations)

    if tolerance <= 0:
        raise ValueError("Solver tolerance must be greater than zero.")

    if maximum_iterations < 1:
        raise ValueError(
            "Maximum solver iterations must be at least one."
        )

    try:
        import openpnm as op
    except ImportError as error:
        raise ImportError(
            "OpenPNM is required for finite-difference diffusion analysis."
        ) from error

    try:
        solver = op.solvers.PyamgRugeStubenSolver(
            tol=tolerance,
            maxiter=maximum_iterations,
        )
    except AttributeError as error:
        raise RuntimeError(
            "The installed OpenPNM version does not provide "
            "PyamgRugeStubenSolver."
        ) from error

    return solver


def _validate_direction(direction: str) -> tuple[str, int]:
    """Validate a physical direction and return its NumPy array axis."""

    direction = str(direction).strip().upper()

    if direction not in DIRECTION_TO_AXIS:
        raise ValueError(
            "Direction must be one of 'X', 'Y', or 'Z'."
        )

    return direction, DIRECTION_TO_AXIS[direction]


def _touches_opposing_faces(
    mask: np.ndarray,
    axis: int,
) -> bool:
    """Return True when a mask touches both opposing faces."""

    minimum_face = np.take(mask, indices=0, axis=axis)
    maximum_face = np.take(mask, indices=-1, axis=axis)

    return bool(np.any(minimum_face) and np.any(maximum_face))


def calculate_diffusion_tortuosity(
    connected_mask: np.ndarray,
    direction: str,
    solver: Any | None = None,
    solver_tolerance: float = 1.0e-10,
    maximum_iterations: int = 5000,
) -> dict[str, Any]:
    """
    Calculate finite-difference diffusion properties in one direction.

    Parameters
    ----------
    connected_mask
        Three-dimensional Boolean or binary mask containing the
        directionally through-connected phase.
    direction
        Physical transport direction: ``"X"``, ``"Y"``, or ``"Z"``.
    solver
        Optional OpenPNM solver. If omitted, a PyAMG Ruge-Stuben solver is
        created automatically.
    solver_tolerance
        Convergence tolerance used when creating the default solver.
    maximum_iterations
        Maximum number of iterations used by the default solver.

    Returns
    -------
    dict
        Direction, array axis, tortuosity, formation factor, relative
        effective diffusivity, effective porosity, original porosity,
        concentration field, and solver run time.
    """

    connected_mask = np.asarray(connected_mask, dtype=bool)
    direction, axis = _validate_direction(direction)

    if connected_mask.ndim != 3:
        raise ValueError(
            f"The connected mask must be three-dimensional; received "
            f"shape {connected_mask.shape}."
        )

    if connected_mask.shape[axis] < 2:
        raise ValueError(
            f"The mask must contain at least two voxels along "
            f"the {direction} direction."
        )

    connected_voxel_count = int(
        np.count_nonzero(connected_mask)
    )

    if connected_voxel_count == 0:
        raise ValueError(
            f"The {direction}-direction connected mask is empty."
        )

    if not _touches_opposing_faces(
        mask=connected_mask,
        axis=axis,
    ):
        raise ValueError(
            f"The supplied mask does not touch both opposing "
            f"{direction} boundary faces."
        )

    if solver is None:
        solver = create_pyamg_solver(
            tolerance=solver_tolerance,
            maximum_iterations=maximum_iterations,
        )

    simulation_result = ps.simulations.tortuosity_fd(
        im=connected_mask,
        axis=axis,
        solver=solver,
    )

    tortuosity = float(simulation_result.tortuosity)
    formation_factor = float(
        simulation_result.formation_factor
    )
    effective_porosity = float(
        simulation_result.effective_porosity
    )
    original_porosity = float(
        simulation_result.original_porosity
    )

    if not np.isfinite(formation_factor):
        raise RuntimeError(
            "The finite-difference solver returned a non-finite "
            "formation factor."
        )

    if formation_factor <= 0:
        raise RuntimeError(
            "The finite-difference solver returned a non-positive "
            "formation factor."
        )

    relative_effective_diffusivity = (
        1.0 / formation_factor
    )

    concentration_field = np.asarray(
        simulation_result.im_conc,
        dtype=np.float32,
    )

    if concentration_field.shape != connected_mask.shape:
        raise RuntimeError(
            "The returned concentration field has an unexpected shape."
        )

    run_time_seconds = float(
        getattr(simulation_result, "time", np.nan)
    )

    return {
        "status": "calculated",
        "direction": direction,
        "axis": axis,
        "connected_voxel_count": connected_voxel_count,
        "original_porosity": original_porosity,
        "effective_porosity": effective_porosity,
        "formation_factor": formation_factor,
        "relative_effective_diffusivity": (
            relative_effective_diffusivity
        ),
        "diffusion_tortuosity": tortuosity,
        "run_time_seconds": run_time_seconds,
        "concentration_field": concentration_field,
    }


def calculate_directional_diffusion(
    volume: np.ndarray,
    phase_label: int,
    directions: Sequence[str] = ("X", "Y", "Z"),
    solver: Any | None = None,
    solver_tolerance: float = 1.0e-10,
    maximum_iterations: int = 5000,
) -> dict[str, Any]:
    """
    Calculate diffusion properties for selected physical directions.

    The function first performs strict 6-connected component analysis.
    A finite-difference simulation is run only when the selected phase
    through-connects the corresponding pair of boundary faces.

    Parameters
    ----------
    volume
        Integer-labelled three-dimensional array in (Z, Y, X) order.
    phase_label
        Integer label representing the transport phase, normally pores.
    directions
        Sequence containing any combination of ``"X"``, ``"Y"``, and ``"Z"``.
    solver
        Optional OpenPNM solver instance. If omitted, one PyAMG solver is
        created and reused.
    solver_tolerance
        Convergence tolerance for the default PyAMG solver.
    maximum_iterations
        Maximum number of iterations for the default PyAMG solver.

    Returns
    -------
    dict
        Connectivity results, directional simulation results, summary table,
        and concentration fields.
    """

    volume = np.asarray(volume)

    if volume.ndim != 3:
        raise ValueError(
            f"The input volume must be three-dimensional; received "
            f"shape {volume.shape}."
        )

    requested_directions = []

    for direction in directions:
        validated_direction, _ = _validate_direction(direction)

        if validated_direction not in requested_directions:
            requested_directions.append(validated_direction)

    if not requested_directions:
        raise ValueError(
            "At least one diffusion direction must be requested."
        )

    connectivity_result = analyse_connectivity(
        volume=volume,
        phase_label=phase_label,
    )

    connectivity_masks = connectivity_result["masks"]
    direction_results: dict[str, dict[str, Any]] = {}
    concentration_fields: dict[str, np.ndarray] = {}
    summary_rows = []

    calculated_solver = solver

    for direction in requested_directions:
        axis = DIRECTION_TO_AXIS[direction]

        connected_mask = connectivity_masks[
            f"through_connected_{direction}"
        ]

        connected_voxel_count = int(
            np.count_nonzero(connected_mask)
        )

        if connected_voxel_count == 0:
            direction_result = {
                "status": "not_percolating",
                "direction": direction,
                "axis": axis,
                "connected_voxel_count": 0,
                "original_porosity": np.nan,
                "effective_porosity": np.nan,
                "formation_factor": np.nan,
                "relative_effective_diffusivity": np.nan,
                "diffusion_tortuosity": np.nan,
                "run_time_seconds": np.nan,
                "concentration_field": None,
            }
        else:
            if calculated_solver is None:
                calculated_solver = create_pyamg_solver(
                    tolerance=solver_tolerance,
                    maximum_iterations=maximum_iterations,
                )

            direction_result = calculate_diffusion_tortuosity(
                connected_mask=connected_mask,
                direction=direction,
                solver=calculated_solver,
                solver_tolerance=solver_tolerance,
                maximum_iterations=maximum_iterations,
            )

            concentration_fields[direction] = direction_result[
                "concentration_field"
            ]

        direction_results[direction] = direction_result

        summary_rows.append(
            {
                "direction": direction,
                "axis": axis,
                "status": direction_result["status"],
                "connected_voxel_count": direction_result[
                    "connected_voxel_count"
                ],
                "directional_connectivity": (
                    connectivity_result["summary"][
                        f"{direction}_connectivity"
                    ]
                ),
                "original_porosity": direction_result[
                    "original_porosity"
                ],
                "effective_porosity": direction_result[
                    "effective_porosity"
                ],
                "formation_factor": direction_result[
                    "formation_factor"
                ],
                "relative_effective_diffusivity": direction_result[
                    "relative_effective_diffusivity"
                ],
                "diffusion_tortuosity": direction_result[
                    "diffusion_tortuosity"
                ],
                "run_time_seconds": direction_result[
                    "run_time_seconds"
                ],
            }
        )

    summary_table = pd.DataFrame(summary_rows)

    return {
        "phase_label": int(phase_label),
        "connectivity_summary": connectivity_result["summary"],
        "direction_results": direction_results,
        "summary": summary_table,
        "concentration_fields": concentration_fields,
    }