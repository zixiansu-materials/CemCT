"""
Unified result-export functions for CemCT.

This module provides a common output layer for all CemCT analysis modules.

Supported outputs
-----------------
1. Binary three-dimensional TIFF masks stored as uint8 using 0/255
2. Scalar three-dimensional TIFF volumes stored as float32
3. Excel workbooks containing one or more result tables
4. JSON files containing analysis metadata
5. Compressed NPZ archives containing numerical arrays
6. Publication-quality PNG figures

Design principles
-----------------
The export layer performs file writing only. Scientific calculations remain
inside the individual analysis modules.

This module contains no:

    input()
    print()
    plt.show()

All created file paths are returned to the caller.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile

from cemct.io import save_binary_tiff


INVALID_EXCEL_SHEET_CHARACTERS = re.compile(r"[\[\]:*?/\\]")
INVALID_FILENAME_CHARACTERS = re.compile(r"[^A-Za-z0-9_.-]+")


def make_safe_name(value: str) -> str:
    """
    Convert a user-supplied name into a filesystem-safe name.
    """
    safe_name = INVALID_FILENAME_CHARACTERS.sub(
        "_",
        str(value).strip(),
    )

    safe_name = safe_name.strip("._")

    if not safe_name:
        raise ValueError("The supplied name does not contain valid characters.")

    return safe_name


def create_output_directory(
    input_file: str | Path,
    output_directory: str | Path | None = None,
) -> Path:
    """
    Create and return the analysis output directory.

    If ``output_directory`` is not supplied, a folder is created beside the
    input TIFF using the following convention:

        input_stem_cemct_results
    """
    input_path = Path(input_file)

    if output_directory is None:
        output_path = (
            input_path.parent
            / f"{input_path.stem}_cemct_results"
        )
    else:
        output_path = Path(output_directory)

    output_path.mkdir(parents=True, exist_ok=True)

    return output_path


def _validate_three_dimensional_array(
    array: np.ndarray,
    array_name: str,
) -> np.ndarray:
    """
    Validate a three-dimensional NumPy array.
    """
    validated_array = np.asarray(array)

    if validated_array.ndim != 3:
        raise ValueError(
            f"{array_name} must be three-dimensional, "
            f"but received shape {validated_array.shape}."
        )

    return validated_array


def export_binary_mask(
    mask: np.ndarray,
    output_path: str | Path,
) -> Path:
    """
    Save a three-dimensional Boolean or binary mask as uint8 TIFF.

    Output convention
    -----------------
    0   = background
    255 = selected phase
    """
    validated_mask = _validate_three_dimensional_array(
        mask,
        "Binary mask",
    )

    return save_binary_tiff(
        mask=np.asarray(validated_mask, dtype=bool),
        output_path=output_path,
    )


def export_scalar_volume(
    volume: np.ndarray,
    output_path: str | Path,
) -> Path:
    """
    Save a three-dimensional scalar field as a float32 TIFF.
    """
    validated_volume = _validate_three_dimensional_array(
        volume,
        "Scalar volume",
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tifffile.imwrite(
        output_path,
        np.asarray(validated_volume, dtype=np.float32),
        photometric="minisblack",
        metadata={"axes": "ZYX"},
    )

    return output_path


def make_safe_excel_sheet_name(
    requested_name: str,
    used_names: set[str],
) -> str:
    """
    Create a valid and unique Excel sheet name.
    """
    sheet_name = INVALID_EXCEL_SHEET_CHARACTERS.sub(
        "_",
        str(requested_name).strip(),
    )

    if not sheet_name:
        sheet_name = "Results"

    sheet_name = sheet_name[:31]

    candidate = sheet_name
    counter = 2

    while candidate.lower() in used_names:
        suffix = f"_{counter}"
        candidate = f"{sheet_name[:31 - len(suffix)]}{suffix}"
        counter += 1

    used_names.add(candidate.lower())

    return candidate


def _metadata_value_for_excel(value: Any) -> Any:
    """
    Convert a metadata value into an Excel-compatible representation.
    """
    converted_value = _to_json_compatible(value)

    if isinstance(converted_value, (dict, list)):
        return json.dumps(
            converted_value,
            ensure_ascii=False,
        )

    return converted_value


def metadata_to_dataframe(
    metadata: Mapping[str, Any],
) -> pd.DataFrame:
    """
    Convert analysis metadata into a two-column DataFrame.
    """
    rows = [
        {
            "Parameter": str(key),
            "Value": _metadata_value_for_excel(value),
        }
        for key, value in metadata.items()
    ]

    return pd.DataFrame(
        rows,
        columns=["Parameter", "Value"],
    )


def export_excel_workbook(
    tables: Mapping[str, pd.DataFrame | pd.Series],
    output_path: str | Path,
    metadata: Mapping[str, Any] | None = None,
) -> Path:
    """
    Save multiple result tables into one Excel workbook.

    Parameters
    ----------
    tables
        Mapping from requested sheet name to pandas DataFrame or Series.
    output_path
        Output .xlsx path.
    metadata
        Optional analysis metadata written to a Metadata worksheet.
    """
    if not tables and metadata is None:
        raise ValueError(
            "At least one result table or metadata mapping is required."
        )

    output_path = Path(output_path)

    if output_path.suffix.lower() != ".xlsx":
        raise ValueError("Excel output path must have an .xlsx extension.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    used_sheet_names: set[str] = set()

    with pd.ExcelWriter(
        output_path,
        engine="openpyxl",
    ) as writer:
        for requested_name, table in tables.items():
            sheet_name = make_safe_excel_sheet_name(
                requested_name,
                used_sheet_names,
            )

            if isinstance(table, pd.Series):
                dataframe = table.to_frame()
            elif isinstance(table, pd.DataFrame):
                dataframe = table.copy()
            else:
                raise TypeError(
                    f"Excel table '{requested_name}' must be a pandas "
                    "DataFrame or Series."
                )

            dataframe.to_excel(
                writer,
                sheet_name=sheet_name,
                index=False,
            )

        if metadata is not None:
            metadata_sheet_name = make_safe_excel_sheet_name(
                "Metadata",
                used_sheet_names,
            )

            metadata_to_dataframe(metadata).to_excel(
                writer,
                sheet_name=metadata_sheet_name,
                index=False,
            )

    return output_path


def _to_json_compatible(value: Any) -> Any:
    """
    Recursively convert values into JSON-compatible Python objects.
    """
    if isinstance(value, Path):
        return str(value)

    if isinstance(value, np.generic):
        return value.item()

    if isinstance(value, np.ndarray):
        return value.tolist()

    if isinstance(value, Mapping):
        return {
            str(key): _to_json_compatible(item)
            for key, item in value.items()
        }

    if isinstance(value, (list, tuple, set)):
        return [
            _to_json_compatible(item)
            for item in value
        ]

    if isinstance(value, (str, int, float, bool)) or value is None:
        return value

    return str(value)


def export_metadata_json(
    metadata: Mapping[str, Any],
    output_path: str | Path,
) -> Path:
    """
    Save analysis metadata as a human-readable JSON file.
    """
    output_path = Path(output_path)

    if output_path.suffix.lower() != ".json":
        raise ValueError("Metadata output path must have a .json extension.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    json_compatible_metadata = _to_json_compatible(metadata)

    with output_path.open(
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(
            json_compatible_metadata,
            file,
            indent=4,
            ensure_ascii=False,
            allow_nan=False,
        )

    return output_path


def export_array_archive(
    arrays: Mapping[str, np.ndarray],
    output_path: str | Path,
) -> Path:
    """
    Save numerical arrays in a compressed NumPy NPZ archive.

    This is suitable for one-dimensional network fields that cannot be saved
    directly as voxel-based three-dimensional TIFF volumes.
    """
    if not arrays:
        raise ValueError("At least one numerical array must be supplied.")

    output_path = Path(output_path)

    if output_path.suffix.lower() != ".npz":
        raise ValueError("Array archive path must have an .npz extension.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    safe_arrays = {
        make_safe_name(name): np.asarray(array)
        for name, array in arrays.items()
    }

    np.savez_compressed(
        output_path,
        **safe_arrays,
    )

    return output_path


def export_figure(
    figure: Any,
    output_path: str | Path,
    dpi: int = 300,
) -> Path:
    """
    Save a Matplotlib Figure without calling plt.show().
    """
    if figure is None:
        raise ValueError("A Matplotlib Figure must be supplied.")

    if not hasattr(figure, "savefig"):
        raise TypeError(
            "figure must provide a Matplotlib-compatible savefig() method."
        )

    if dpi <= 0:
        raise ValueError("dpi must be greater than zero.")

    output_path = Path(output_path)

    if output_path.suffix.lower() != ".png":
        raise ValueError("Figure output path must have a .png extension.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure.savefig(
        output_path,
        dpi=int(dpi),
        bbox_inches="tight",
        facecolor="white",
    )

    return output_path


def export_analysis_bundle(
    input_file: str | Path,
    analysis_name: str,
    tables: Mapping[str, pd.DataFrame | pd.Series] | None = None,
    masks: Mapping[str, np.ndarray] | None = None,
    scalar_volumes: Mapping[str, np.ndarray] | None = None,
    arrays: Mapping[str, np.ndarray] | None = None,
    metadata: Mapping[str, Any] | None = None,
    figure: Any | None = None,
    output_directory: str | Path | None = None,
    figure_dpi: int = 300,
) -> dict[str, Any]:
    """
    Export a complete CemCT analysis-result bundle.

    Returns
    -------
    dict
        Paths to all files created by this function.
    """
    input_path = Path(input_file)
    safe_analysis_name = make_safe_name(analysis_name)

    output_directory_path = create_output_directory(
        input_file=input_path,
        output_directory=output_directory,
    )

    base_name = (
        f"{make_safe_name(input_path.stem)}_"
        f"{safe_analysis_name}"
    )

    created_paths: dict[str, Any] = {
        "output_directory": output_directory_path,
        "binary_masks": {},
        "scalar_volumes": {},
    }

    if tables or metadata is not None:
        excel_path = (
            output_directory_path
            / f"{base_name}_results.xlsx"
        )

        created_paths["excel"] = export_excel_workbook(
            tables={} if tables is None else tables,
            output_path=excel_path,
            metadata=metadata,
        )

    if metadata is not None:
        json_path = (
            output_directory_path
            / f"{base_name}_metadata.json"
        )

        created_paths["metadata_json"] = export_metadata_json(
            metadata=metadata,
            output_path=json_path,
        )

    if masks:
        for mask_name, mask in masks.items():
            safe_mask_name = make_safe_name(mask_name)

            mask_path = (
                output_directory_path
                / f"{base_name}_{safe_mask_name}.tif"
            )

            created_paths["binary_masks"][mask_name] = (
                export_binary_mask(
                    mask=mask,
                    output_path=mask_path,
                )
            )

    if scalar_volumes:
        for volume_name, scalar_volume in scalar_volumes.items():
            safe_volume_name = make_safe_name(volume_name)

            volume_path = (
                output_directory_path
                / f"{base_name}_{safe_volume_name}.tif"
            )

            created_paths["scalar_volumes"][volume_name] = (
                export_scalar_volume(
                    volume=scalar_volume,
                    output_path=volume_path,
                )
            )

    if arrays:
        array_archive_path = (
            output_directory_path
            / f"{base_name}_arrays.npz"
        )

        created_paths["array_archive"] = export_array_archive(
            arrays=arrays,
            output_path=array_archive_path,
        )

    if figure is not None:
        figure_path = (
            output_directory_path
            / f"{base_name}_summary.png"
        )

        created_paths["figure"] = export_figure(
            figure=figure,
            output_path=figure_path,
            dpi=figure_dpi,
        )

    return created_paths