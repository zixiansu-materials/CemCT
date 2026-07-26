"""
Intrinsic permeability analysis for labelled 3D CT data.

This module extracts direction-specific pore networks from previously
identified through-connected regions and calculates intrinsic permeability
using steady-state Stokes-flow simulations.

Workflow
--------
1. Select one integer label from a labelled 3D volume.
2. Identify X-, Y-, and Z-through-connected regions using strict
   face-connected 6-connectivity.
3. Extract a pore-throat network using PoreSpy SNOW2.
4. Import the extracted network into OpenPNM.
5. Assign hydraulic size-factor and hydraulic-conductance models.
6. Solve steady-state Stokes flow.
7. Calculate intrinsic permeability using Darcy's law.

Coordinate convention
---------------------
Input NumPy volume:

    axis 0 = Z
    axis 1 = Y
    axis 2 = X

Before SNOW2 extraction, the selected mask is transposed internally from
(Z, Y, X) to physical (X, Y, Z) order. OpenPNM boundary labels therefore
retain their physical meanings:

    pore.xmin / pore.xmax = X boundaries
    pore.ymin / pore.ymax = Y boundaries
    pore.zmin / pore.zmax = Z boundaries

The core functions contain no input(), print(), plt.show(), TIFF writing,
or Excel writing. User interaction and result export belong in separate
interface and export layers.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import openpnm as op
import pandas as pd
import porespy as ps
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from cemct.connectivity import analyse_connectivity


DIRECTION_TO_MASK_KEY = {
    "X": "through_connected_X",
    "Y": "through_connected_Y",
    "Z": "through_connected_Z",
}

DIRECTION_TO_BOUNDARY_LABELS = {
    "X": ("xmin", "xmax"),
    "Y": ("ymin", "ymax"),
    "Z": ("zmin", "zmax"),
}

DIRECTION_TO_BOUNDARY_WIDTH = {
    "X": [3, 0, 0],
    "Y": [0, 3, 0],
    "Z": [0, 0, 3],
}

DIRECTION_TO_PHYSICAL_SHAPE_INDEX = {
    "X": 0,
    "Y": 1,
    "Z": 2,
}

DARCY_IN_SQUARE_METRES = 9.869233e-13


def validate_labelled_volume(
    volume: np.ndarray,
    phase_label: int,
) -> np.ndarray:
    """
    Validate an integer-labelled three-dimensional volume.
    """
    array = np.asarray(volume)

    if array.ndim != 3:
        raise ValueError(
            f"Expected a 3D volume, but received shape {array.shape}."
        )

    if not np.issubdtype(array.dtype, np.integer):
        raise TypeError(
            f"Expected an integer-labelled volume, but received {array.dtype}."
        )

    if phase_label not in np.unique(array):
        raise ValueError(
            f"Selected label {phase_label} is absent from the input volume."
        )

    return array


def normalise_directions(
    directions: Sequence[str],
) -> tuple[str, ...]:
    """
    Validate and standardise requested physical directions.
    """
    normalised = tuple(str(direction).upper() for direction in directions)

    if not normalised:
        raise ValueError("At least one analysis direction must be provided.")

    invalid = set(normalised) - {"X", "Y", "Z"}

    if invalid:
        raise ValueError(
            f"Invalid directions: {sorted(invalid)}. "
            "Allowed directions are X, Y, and Z."
        )

    if len(set(normalised)) != len(normalised):
        raise ValueError("Analysis directions must not contain duplicates.")

    return normalised


def to_physical_xyz(mask_zyx: np.ndarray) -> np.ndarray:
    """
    Convert a mask from NumPy (Z, Y, X) order to physical (X, Y, Z) order.
    """
    mask = np.asarray(mask_zyx, dtype=bool)

    if mask.ndim != 3:
        raise ValueError("The connected-phase mask must be three-dimensional.")

    return np.transpose(mask, (2, 1, 0))


def get_boundary_width(
    direction: str,
    boundary_width: int,
) -> list[int]:
    """
    Return the SNOW2 boundary width in physical (X, Y, Z) order.
    """
    direction = direction.upper()

    if direction not in DIRECTION_TO_BOUNDARY_WIDTH:
        raise ValueError("Direction must be X, Y, or Z.")

    if not isinstance(boundary_width, int) or boundary_width < 1:
        raise ValueError("boundary_width must be a positive integer.")

    widths = [0, 0, 0]
    widths[DIRECTION_TO_PHYSICAL_SHAPE_INDEX[direction]] = boundary_width

    return widths


def get_domain_geometry(
    original_shape_zyx: Sequence[int],
    direction: str,
    voxel_size_m: float,
    boundary_width: int,
) -> tuple[float, float]:
    """
    Calculate simulation length and cross-sectional area.

    The SNOW2 boundary regions extend the simulation domain by
    ``boundary_width`` voxels at each opposing boundary.
    """
    if len(original_shape_zyx) != 3:
        raise ValueError("original_shape_zyx must contain Z, Y, and X sizes.")

    if voxel_size_m <= 0:
        raise ValueError("voxel_size_m must be greater than zero.")

    direction = direction.upper()

    if direction not in {"X", "Y", "Z"}:
        raise ValueError("Direction must be X, Y, or Z.")

    shape_zyx = tuple(int(value) for value in original_shape_zyx)
    shape_xyz = (
        shape_zyx[2],
        shape_zyx[1],
        shape_zyx[0],
    )

    flow_axis = DIRECTION_TO_PHYSICAL_SHAPE_INDEX[direction]
    transverse_axes = tuple(
        axis for axis in range(3) if axis != flow_axis
    )

    length_voxels = shape_xyz[flow_axis] + 2 * boundary_width

    area_voxels = (
        shape_xyz[transverse_axes[0]]
        * shape_xyz[transverse_axes[1]]
    )

    length_m = length_voxels * voxel_size_m
    area_m2 = area_voxels * voxel_size_m**2

    return float(length_m), float(area_m2)


def calculate_darcy_permeability(
    flow_rate_m3_s: float,
    dynamic_viscosity_pa_s: float,
    sample_length_m: float,
    cross_sectional_area_m2: float,
    pressure_drop_pa: float,
) -> float:
    """
    Calculate intrinsic permeability using Darcy's law.

    Returns
    -------
    float
        Intrinsic permeability in square metres.
    """
    values = {
        "flow_rate_m3_s": flow_rate_m3_s,
        "dynamic_viscosity_pa_s": dynamic_viscosity_pa_s,
        "sample_length_m": sample_length_m,
        "cross_sectional_area_m2": cross_sectional_area_m2,
        "pressure_drop_pa": pressure_drop_pa,
    }

    for name, value in values.items():
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite.")

    if flow_rate_m3_s < 0:
        raise ValueError("flow_rate_m3_s must not be negative.")

    if dynamic_viscosity_pa_s <= 0:
        raise ValueError(
            "dynamic_viscosity_pa_s must be greater than zero."
        )

    if sample_length_m <= 0:
        raise ValueError("sample_length_m must be greater than zero.")

    if cross_sectional_area_m2 <= 0:
        raise ValueError(
            "cross_sectional_area_m2 must be greater than zero."
        )

    if pressure_drop_pa <= 0:
        raise ValueError("pressure_drop_pa must be greater than zero.")

    permeability_m2 = (
        flow_rate_m3_s
        * dynamic_viscosity_pa_s
        * sample_length_m
        / (cross_sectional_area_m2 * pressure_drop_pa)
    )

    return float(permeability_m2)


def extract_directional_network(
    connected_mask_zyx: np.ndarray,
    direction: str,
    voxel_size_um: float,
    boundary_width: int = 3,
    accuracy: str = "standard",
    sigma: float = 0.4,
    r_max: int = 4,
) -> tuple[Any, Any]:
    """
    Extract an OpenPNM network from one directional through-connected mask.

    Returns
    -------
    tuple
        ``(openpnm_network, snow2_result)``
    """
    connected_mask = np.asarray(connected_mask_zyx, dtype=bool)

    if connected_mask.ndim != 3:
        raise ValueError("connected_mask_zyx must be three-dimensional.")

    if not np.any(connected_mask):
        raise ValueError(
            f"No through-connected voxels were found in the {direction} "
            "direction."
        )

    if voxel_size_um <= 0:
        raise ValueError("voxel_size_um must be greater than zero.")

    direction = direction.upper()
    boundary_width_xyz = get_boundary_width(direction, boundary_width)

    # PoreSpy/OpenPNM spatial coordinates are treated as physical X, Y, Z.
    mask_xyz = to_physical_xyz(connected_mask).astype(np.uint8)

    voxel_size_m = float(voxel_size_um) * 1e-6

    snow_result = ps.networks.snow2(
        phases=mask_xyz,
        phase_alias={1: "selected_phase"},
        boundary_width=boundary_width_xyz,
        accuracy=accuracy,
        voxel_size=voxel_size_m,
        sigma=float(sigma),
        r_max=int(r_max),
        parallel_kw=None,
    )

    network = op.io.network_from_porespy(snow_result.network)

    return network, snow_result


def retain_inlet_to_outlet_components(
    network: Any,
    direction: str,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Remove network components that do not connect inlet to outlet.

    Multiple independent through-connected pathways are retained because
    they represent parallel transport pathways.
    """
    direction = direction.upper()
    inlet_label, outlet_label = DIRECTION_TO_BOUNDARY_LABELS[direction]

    inlet_pores = np.asarray(
        network.pores(inlet_label),
        dtype=int,
    )
    outlet_pores = np.asarray(
        network.pores(outlet_label),
        dtype=int,
    )

    if inlet_pores.size == 0:
        raise RuntimeError(
            f"The extracted network has no pore.{inlet_label} pores."
        )

    if outlet_pores.size == 0:
        raise RuntimeError(
            f"The extracted network has no pore.{outlet_label} pores."
        )

    number_of_pores = int(network.Np)
    throat_connections = np.asarray(network["throat.conns"], dtype=int)

    if throat_connections.size == 0:
        raise RuntimeError("The extracted network contains no throats.")

    row = np.concatenate(
        [throat_connections[:, 0], throat_connections[:, 1]]
    )
    column = np.concatenate(
        [throat_connections[:, 1], throat_connections[:, 0]]
    )
    data = np.ones(row.size, dtype=np.uint8)

    adjacency = coo_matrix(
        (data, (row, column)),
        shape=(number_of_pores, number_of_pores),
    ).tocsr()

    number_of_components, component_labels = connected_components(
        adjacency,
        directed=False,
        return_labels=True,
    )

    inlet_components = set(component_labels[inlet_pores].tolist())
    outlet_components = set(component_labels[outlet_pores].tolist())

    retained_components = inlet_components.intersection(outlet_components)

    if not retained_components:
        raise RuntimeError(
            f"No pore-network component spans the {direction} direction."
        )

    retain_mask = np.isin(
        component_labels,
        list(retained_components),
    )
    pores_to_remove = np.flatnonzero(~retain_mask)

    if pores_to_remove.size:
        op.topotools.trim(
            network=network,
            pores=pores_to_remove,
        )

    inlet_pores = np.asarray(
        network.pores(inlet_label),
        dtype=int,
    )
    outlet_pores = np.asarray(
        network.pores(outlet_label),
        dtype=int,
    )

    if inlet_pores.size == 0 or outlet_pores.size == 0:
        raise RuntimeError(
            "Boundary pores were lost during network cleaning."
        )

    return inlet_pores, outlet_pores


def assign_hydraulic_models(
    network: Any,
    dynamic_viscosity_pa_s: float,
) -> Any:
    """
    Add only models required for hydraulic flow.

    Unrelated capillary-entry-pressure and diffusive-conductance models are
    deliberately not added.
    """
    if dynamic_viscosity_pa_s <= 0:
        raise ValueError(
            "dynamic_viscosity_pa_s must be greater than zero."
        )

    network["pore.diameter"] = np.asarray(
        network["pore.equivalent_diameter"],
        dtype=float,
    )

    network["throat.diameter"] = np.asarray(
        network["throat.inscribed_diameter"],
        dtype=float,
    )

    network["throat.spacing"] = np.asarray(
        network["throat.total_length"],
        dtype=float,
    )

    network.add_model(
        propname="throat.hydraulic_size_factors",
        model=(
            op.models.geometry.hydraulic_size_factors
            .pyramids_and_cuboids
        ),
    )
    network.regenerate_models()

    fluid = op.phase.Phase(network=network)

    fluid["pore.viscosity"] = float(dynamic_viscosity_pa_s)

    fluid.add_model(
        propname="throat.hydraulic_conductance",
        model=(
            op.models.physics.hydraulic_conductance
            .generic_hydraulic
        ),
    )
    fluid.regenerate_models()

    return fluid


def solve_stokes_flow(
    network: Any,
    fluid: Any,
    inlet_pores: np.ndarray,
    outlet_pores: np.ndarray,
    inlet_pressure_pa: float = 1.0,
    outlet_pressure_pa: float = 0.0,
    solver_tolerance: float = 1e-10,
    maximum_iterations: int = 5000,
) -> tuple[Any, float, np.ndarray]:
    """
    Solve steady-state Stokes flow using a PyAMG Ruge-Stuben solver.
    """
    if inlet_pressure_pa <= outlet_pressure_pa:
        raise ValueError(
            "inlet_pressure_pa must be greater than outlet_pressure_pa."
        )

    flow = op.algorithms.StokesFlow(
        network=network,
        phase=fluid,
    )

    flow.set_value_BC(
        pores=inlet_pores,
        values=float(inlet_pressure_pa),
    )
    flow.set_value_BC(
        pores=outlet_pores,
        values=float(outlet_pressure_pa),
    )

    solver = op.solvers.PyamgRugeStubenSolver(
        tol=float(solver_tolerance),
        maxiter=int(maximum_iterations),
    )

    flow.run(solver=solver)

    flow_rate = flow.rate(
        pores=inlet_pores,
        mode="group",
    )

    flow_rate_m3_s = abs(
        float(np.asarray(flow_rate).reshape(-1)[0])
    )

    pressure_field_pa = np.asarray(
        flow["pore.pressure"],
        dtype=float,
    )

    return flow, flow_rate_m3_s, pressure_field_pa


def analyse_permeability(
    volume: np.ndarray,
    phase_label: int,
    voxel_size_um: float,
    directions: Sequence[str] = ("X", "Y", "Z"),
    boundary_width: int = 3,
    dynamic_viscosity_pa_s: float = 1.0,
    inlet_pressure_pa: float = 1.0,
    outlet_pressure_pa: float = 0.0,
    accuracy: str = "standard",
    sigma: float = 0.4,
    r_max: int = 4,
    solver_tolerance: float = 1e-10,
    maximum_iterations: int = 5000,
) -> dict[str, Any]:
    """
    Calculate direction-dependent intrinsic permeability.

    Notes
    -----
    A viscosity of 1 Pa.s and pressure drop of 1 Pa are convenient numerical
    reference values. Intrinsic permeability is independent of the chosen
    reference viscosity when the hydraulic conductance and Darcy conversion
    are applied consistently.
    """
    labelled_volume = validate_labelled_volume(volume, phase_label)
    requested_directions = normalise_directions(directions)

    if voxel_size_um <= 0:
        raise ValueError("voxel_size_um must be greater than zero.")

    if dynamic_viscosity_pa_s <= 0:
        raise ValueError(
            "dynamic_viscosity_pa_s must be greater than zero."
        )

    pressure_drop_pa = inlet_pressure_pa - outlet_pressure_pa

    if pressure_drop_pa <= 0:
        raise ValueError(
            "Inlet pressure must be greater than outlet pressure."
        )

    connectivity_result = analyse_connectivity(
        labelled_volume,
        phase_label,
    )

    selected_phase_voxels = int(
        np.count_nonzero(labelled_volume == phase_label)
    )

    voxel_size_m = float(voxel_size_um) * 1e-6

    result_rows: list[dict[str, Any]] = []
    networks: dict[str, Any] = {}
    snow_results: dict[str, Any] = {}
    flow_algorithms: dict[str, Any] = {}
    pressure_fields: dict[str, np.ndarray] = {}

    for direction in requested_directions:
        mask_key = DIRECTION_TO_MASK_KEY[direction]

        connected_mask = np.asarray(
            connectivity_result["masks"][mask_key],
            dtype=bool,
        )

        connected_voxels = int(np.count_nonzero(connected_mask))

        if connected_voxels == 0:
            result_rows.append(
                {
                    "Direction": direction,
                    "Status": "No through-connected phase",
                    "Selected phase voxels": selected_phase_voxels,
                    "Through-connected voxels": 0,
                    "Directional connectivity (%)": 0.0,
                    "Network pores": 0,
                    "Network throats": 0,
                    "Flow rate (m^3/s)": np.nan,
                    "Pressure drop (Pa)": pressure_drop_pa,
                    "Simulation length (m)": np.nan,
                    "Cross-sectional area (m^2)": np.nan,
                    "Intrinsic permeability (m^2)": np.nan,
                    "Intrinsic permeability (Darcy)": np.nan,
                    "Intrinsic permeability (mD)": np.nan,
                }
            )
            continue

        network, snow_result = extract_directional_network(
            connected_mask_zyx=connected_mask,
            direction=direction,
            voxel_size_um=voxel_size_um,
            boundary_width=boundary_width,
            accuracy=accuracy,
            sigma=sigma,
            r_max=r_max,
        )

        inlet_pores, outlet_pores = retain_inlet_to_outlet_components(
            network=network,
            direction=direction,
        )

        fluid = assign_hydraulic_models(
            network=network,
            dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
        )

        flow, flow_rate_m3_s, pressure_field_pa = solve_stokes_flow(
            network=network,
            fluid=fluid,
            inlet_pores=inlet_pores,
            outlet_pores=outlet_pores,
            inlet_pressure_pa=inlet_pressure_pa,
            outlet_pressure_pa=outlet_pressure_pa,
            solver_tolerance=solver_tolerance,
            maximum_iterations=maximum_iterations,
        )

        length_m, area_m2 = get_domain_geometry(
            original_shape_zyx=labelled_volume.shape,
            direction=direction,
            voxel_size_m=voxel_size_m,
            boundary_width=boundary_width,
        )

        permeability_m2 = calculate_darcy_permeability(
            flow_rate_m3_s=flow_rate_m3_s,
            dynamic_viscosity_pa_s=dynamic_viscosity_pa_s,
            sample_length_m=length_m,
            cross_sectional_area_m2=area_m2,
            pressure_drop_pa=pressure_drop_pa,
        )

        permeability_darcy = (
            permeability_m2 / DARCY_IN_SQUARE_METRES
        )
        permeability_md = permeability_darcy * 1000.0

        connectivity_percent = (
            100.0 * connected_voxels / selected_phase_voxels
        )

        result_rows.append(
            {
                "Direction": direction,
                "Status": "Solved",
                "Selected phase voxels": selected_phase_voxels,
                "Through-connected voxels": connected_voxels,
                "Directional connectivity (%)": connectivity_percent,
                "Network pores": int(network.Np),
                "Network throats": int(network.Nt),
                "Flow rate (m^3/s)": flow_rate_m3_s,
                "Pressure drop (Pa)": pressure_drop_pa,
                "Simulation length (m)": length_m,
                "Cross-sectional area (m^2)": area_m2,
                "Intrinsic permeability (m^2)": permeability_m2,
                "Intrinsic permeability (Darcy)": permeability_darcy,
                "Intrinsic permeability (mD)": permeability_md,
            }
        )

        networks[direction] = network
        snow_results[direction] = snow_result
        flow_algorithms[direction] = flow
        pressure_fields[direction] = pressure_field_pa

    summary = pd.DataFrame(result_rows)

    return {
        "summary": summary,
        "connectivity": connectivity_result,
        "networks": networks,
        "snow_results": snow_results,
        "flow_algorithms": flow_algorithms,
        "pressure_fields": pressure_fields,
        "metadata": {
            "phase_label": int(phase_label),
            "voxel_size_um": float(voxel_size_um),
            "directions": requested_directions,
            "boundary_width": int(boundary_width),
            "dynamic_viscosity_pa_s": float(dynamic_viscosity_pa_s),
            "inlet_pressure_pa": float(inlet_pressure_pa),
            "outlet_pressure_pa": float(outlet_pressure_pa),
            "pressure_drop_pa": float(pressure_drop_pa),
            "accuracy": accuracy,
            "sigma": float(sigma),
            "r_max": int(r_max),
            "solver": "PyAMG Ruge-Stuben",
            "solver_tolerance": float(solver_tolerance),
            "maximum_iterations": int(maximum_iterations),
        },
    }