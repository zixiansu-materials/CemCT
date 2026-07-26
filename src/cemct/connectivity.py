"""
Connectivity analysis for selected phases in labelled 3D CT data.

Connectivity definition
-----------------------
Only face-connected voxels are considered connected.

Face connection   = included
Edge connection   = excluded
Corner connection = excluded

This corresponds to 6-connectivity in three dimensions.

Direction convention
--------------------
The input NumPy array is interpreted in (Z, Y, X) order:

    X direction = array axis 2
    Y direction = array axis 1
    Z direction = array axis 0

A connected component is considered directionally through-connected when it
touches both opposing boundary faces in that direction.
"""

from typing import Any

import numpy as np
from scipy import ndimage


def _component_ids_touching_face(
    labelled_components: np.ndarray,
    axis: int,
    index: int,
) -> np.ndarray:
    """Return non-zero component IDs touching one specified boundary face."""

    face = np.take(labelled_components, indices=index, axis=axis)
    component_ids = np.unique(face)

    return component_ids[component_ids != 0]


def _mask_from_component_ids(
    labelled_components: np.ndarray,
    component_ids: np.ndarray,
) -> np.ndarray:
    """Create a Boolean mask containing the selected component IDs."""

    if component_ids.size == 0:
        return np.zeros(labelled_components.shape, dtype=bool)

    return np.isin(labelled_components, component_ids)


def _through_connected_component_ids(
    labelled_components: np.ndarray,
    axis: int,
) -> np.ndarray:
    """
    Return component IDs touching both opposing faces of one direction.
    """

    minimum_face_ids = _component_ids_touching_face(
        labelled_components=labelled_components,
        axis=axis,
        index=0,
    )

    maximum_face_ids = _component_ids_touching_face(
        labelled_components=labelled_components,
        axis=axis,
        index=-1,
    )

    return np.intersect1d(
        minimum_face_ids,
        maximum_face_ids,
        assume_unique=False,
    )


def analyse_connectivity(
    volume: np.ndarray,
    phase_label: int,
) -> dict[str, Any]:
    """
    Analyse connectivity of one selected label in a labelled 3D volume.

    Parameters
    ----------
    volume
        Integer-labelled three-dimensional NumPy array in (Z, Y, X) order.
    phase_label
        Integer label representing the phase to be analysed, normally pores.

    Returns
    -------
    dict
        Dictionary containing numerical connectivity results and Boolean masks.

    Notes
    -----
    Directional connectivity is defined as:

        directionally through-connected selected voxels
        ------------------------------------------------
                 total selected-phase voxels

    Boundary-accessible pores touch at least one of the six external faces.

    Isolated pores do not touch any external boundary face.
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
    selected_mask = volume == phase_label
    selected_voxel_count = int(np.count_nonzero(selected_mask))

    if selected_voxel_count == 0:
        raise ValueError(
            f"Selected label {phase_label} is absent from the input volume."
        )

    # Three-dimensional 6-connectivity:
    # face neighbours are included; edge and corner neighbours are excluded.
    connectivity_structure = ndimage.generate_binary_structure(
        rank=3,
        connectivity=1,
    )

    labelled_components, component_count = ndimage.label(
        selected_mask,
        structure=connectivity_structure,
    )

    # NumPy volume order is (Z, Y, X).
    direction_axes = {
        "X": 2,
        "Y": 1,
        "Z": 0,
    }

    directional_masks: dict[str, np.ndarray] = {}
    directional_component_counts: dict[str, int] = {}

    for direction, axis in direction_axes.items():
        component_ids = _through_connected_component_ids(
            labelled_components=labelled_components,
            axis=axis,
        )

        directional_masks[direction] = _mask_from_component_ids(
            labelled_components=labelled_components,
            component_ids=component_ids,
        )

        directional_component_counts[direction] = int(component_ids.size)

    # Components touching any one of the six external faces.
    boundary_component_ids = np.unique(
        np.concatenate(
            [
                _component_ids_touching_face(labelled_components, 0, 0),
                _component_ids_touching_face(labelled_components, 0, -1),
                _component_ids_touching_face(labelled_components, 1, 0),
                _component_ids_touching_face(labelled_components, 1, -1),
                _component_ids_touching_face(labelled_components, 2, 0),
                _component_ids_touching_face(labelled_components, 2, -1),
            ]
        )
    )

    boundary_accessible_mask = _mask_from_component_ids(
        labelled_components=labelled_components,
        component_ids=boundary_component_ids,
    )

    isolated_mask = selected_mask & ~boundary_accessible_mask

    any_direction_mask = (
        directional_masks["X"]
        | directional_masks["Y"]
        | directional_masks["Z"]
    )

    boundary_accessible_non_through_mask = (
        boundary_accessible_mask & ~any_direction_mask
    )

    x_voxel_count = int(np.count_nonzero(directional_masks["X"]))
    y_voxel_count = int(np.count_nonzero(directional_masks["Y"]))
    z_voxel_count = int(np.count_nonzero(directional_masks["Z"]))
    any_direction_voxel_count = int(np.count_nonzero(any_direction_mask))
    accessible_voxel_count = int(np.count_nonzero(boundary_accessible_mask))
    isolated_voxel_count = int(np.count_nonzero(isolated_mask))

    summary = {
        "phase_label": phase_label,
        "selected_voxel_count": selected_voxel_count,
        "connected_component_count": int(component_count),
        "boundary_accessible_component_count": int(
            boundary_component_ids.size
        ),
        "boundary_accessible_voxel_count": accessible_voxel_count,
        "isolated_voxel_count": isolated_voxel_count,
        "X_connected_component_count": directional_component_counts["X"],
        "Y_connected_component_count": directional_component_counts["Y"],
        "Z_connected_component_count": directional_component_counts["Z"],
        "X_connected_voxel_count": x_voxel_count,
        "Y_connected_voxel_count": y_voxel_count,
        "Z_connected_voxel_count": z_voxel_count,
        "any_direction_connected_voxel_count": any_direction_voxel_count,
        "boundary_accessible_fraction": (
            accessible_voxel_count / selected_voxel_count
        ),
        "isolated_fraction": isolated_voxel_count / selected_voxel_count,
        "X_connectivity": x_voxel_count / selected_voxel_count,
        "Y_connectivity": y_voxel_count / selected_voxel_count,
        "Z_connectivity": z_voxel_count / selected_voxel_count,
        "any_direction_connectivity": (
            any_direction_voxel_count / selected_voxel_count
        ),
    }

    masks = {
        "selected": selected_mask,
        "boundary_accessible": boundary_accessible_mask,
        "boundary_accessible_non_through": (
            boundary_accessible_non_through_mask
        ),
        "isolated": isolated_mask,
        "through_connected_X": directional_masks["X"],
        "through_connected_Y": directional_masks["Y"],
        "through_connected_Z": directional_masks["Z"],
        "through_connected_any_direction": any_direction_mask,
    }

    return {
        "summary": summary,
        "masks": masks,
        "labelled_components": labelled_components,
    }