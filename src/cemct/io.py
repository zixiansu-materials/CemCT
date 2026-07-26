"""Input and output functions for CemCT."""

from pathlib import Path

import numpy as np
import tifffile


def load_labelled_tiff(file_path: str | Path) -> np.ndarray:
    """
    Load an integer-labelled three-dimensional TIFF volume.

    Parameters
    ----------
    file_path
        Path to a multi-page 3D TIFF file.

    Returns
    -------
    numpy.ndarray
        Labelled volume in (Z, Y, X) array order.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"TIFF file not found: {file_path}")

    if file_path.suffix.lower() not in {".tif", ".tiff"}:
        raise ValueError("Input file must have a .tif or .tiff extension.")

    volume = tifffile.imread(file_path)

    if volume.ndim != 3:
        raise ValueError(
            f"Expected a 3D TIFF volume, but received shape {volume.shape}."
        )

    if not np.issubdtype(volume.dtype, np.integer):
        raise TypeError(
            f"Expected an integer-labelled volume, but received {volume.dtype}."
        )

    return volume


def save_binary_tiff(
    mask: np.ndarray,
    output_path: str | Path,
) -> Path:
    """
    Save a Boolean mask as an unsigned 8-bit binary TIFF.

    Output convention
    -----------------
    0   = background
    255 = selected phase
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    binary_volume = np.asarray(mask, dtype=bool).astype(np.uint8) * 255

    tifffile.imwrite(
        output_path,
        binary_volume,
        photometric="minisblack",
    )

    return output_path