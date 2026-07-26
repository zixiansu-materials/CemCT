# CemCT

[![CemCT automated tests](https://github.com/zixiansu-materials/CemCT/actions/workflows/tests.yml/badge.svg)](https://github.com/zixiansu-materials/CemCT/actions/workflows/tests.yml)

CemCT is a Python package for the post-segmentation analysis of
three-dimensional X-ray micro-computed tomography (XRµCT) datasets of
cementitious materials.

The package supports binary and integer-labelled 3D TIFF volumes and provides
a reproducible workflow for quantitative phase, pore-structure and
transport-property analysis.

> **Development status:** Release Candidate  
> **Current version:** 1.0.0  
> **Release status:** Unreleased working version  
>
> The present version is being prepared for open-source release. It should not
> be redistributed without permission from the development team.

---

## Overview

Quantitative analysis of segmented XRµCT datasets commonly requires multiple
software packages, manual file conversion and repeated implementation of
similar image-processing procedures.

CemCT integrates complementary voxel-based, image-based, random-walk and
pore-network modelling methods into a common Python workflow.

The package is intended for the analysis of:

- Cement paste
- Mortar
- Concrete
- Alkali-activated materials
- Supplementary cementitious materials
- Other segmented porous and multi-phase materials

CemCT starts from reconstructed and segmented image volumes. Image
acquisition, tomographic reconstruction and segmentation are outside the
scope of the current package.

---

## Main capabilities

CemCT currently provides:

- User-defined multi-phase fraction analysis
- Isolated pore analysis
- Boundary-accessible pore analysis
- Directionally through-connected pore analysis along X, Y and Z
- Local-thickness-based pore-size distribution analysis
- Pore-size-class contribution analysis
- Cumulative pore-volume analysis
- Directional finite-difference diffusion analysis
- Formation-factor calculation
- Relative-effective-diffusivity calculation
- Diffusion-tortuosity calculation
- Random-walk-based transport analysis
- Pore-network extraction using SNOW2
- Directional intrinsic-permeability estimation
- Unified export of TIFF, Excel, JSON, NPZ and PNG results

---

## Analysis modules

| Module | Main function | Description |
|---|---|---|
| Input/output | `load_labelled_tiff()` | Reads integer-labelled 3D TIFF volumes |
| Phase fraction | `calculate_phase_fractions()` | Calculates user-defined phase fractions and physical volumes |
| Connectivity | `analyse_connectivity()` | Classifies isolated, accessible and directionally connected regions |
| Pore size | `calculate_pore_size_distribution()` | Calculates local-thickness-based pore-size distributions |
| Diffusion | `calculate_directional_diffusion()` | Calculates directional diffusion transport properties |
| Random walk | `calculate_connected_random_walk()` | Performs random-walk-based transport analysis |
| Permeability | `analyse_permeability()` | Extracts pore networks and estimates intrinsic permeability |
| Export | `export_analysis_bundle()` | Exports numerical, image and metadata results |

---

## Input data

The primary input is a three-dimensional TIFF volume containing integer phase
labels.

Example:

| Label | Phase |
|---:|---|
| 1 | Unreacted particle |
| 2 | Reaction product |
| 3 | Pore |
| 4 | Other phase |

The number of labels, label values and phase names are defined by the user.
CemCT is therefore not restricted to a fixed four-phase material system.

Binary volumes using `0/1`, `False/True` or `0/255` conventions can also be
converted internally to Boolean masks after explicit selection of the phase
to be analysed.

### Input requirements

The input dataset should satisfy the following requirements:

- The image must be a three-dimensional multi-page TIFF volume.
- Label values must be stored using an integer data type.
- Each voxel must contain one segmented phase label.
- The physical voxel dimensions must be provided in micrometres.
- The user must identify the label representing the phase to be analysed.
- The volume should already be reconstructed, corrected and segmented.

---

## Coordinate convention

CemCT interprets three-dimensional NumPy arrays in the following order:

```text
Array axis 0 = Z direction
Array axis 1 = Y direction
Array axis 2 = X direction
```

The directional correspondence used throughout the package is:

| Physical direction | NumPy axis | Opposing faces |
|---|---:|---|
| X | 2 | X-min and X-max |
| Y | 1 | Y-min and Y-max |
| Z | 0 | Z-min and Z-max |

All directional image analysis, connectivity, diffusion and permeability
results use this convention.

---

## Connectivity definition

Connectivity is calculated using strict 6-connectivity in three dimensions:

```text
Face connection   = included
Edge connection   = excluded
Corner connection = excluded
```

A connected component is considered directionally through-connected when it
touches both opposing boundary faces:

```text
X-connected = touches both X-min and X-max
Y-connected = touches both Y-min and Y-max
Z-connected = touches both Z-min and Z-max
```

Directional connectivity is calculated as:

```text
Directional connectivity =
directionally through-connected phase voxels / total selected-phase voxels
```

The connectivity module can distinguish:

- Isolated regions
- Boundary-accessible regions
- Boundary-accessible but non-through-connected regions
- X-directionally through-connected regions
- Y-directionally through-connected regions
- Z-directionally through-connected regions
- Regions connected through at least one direction

---

## Phase-fraction analysis

For a selected phase `i`, the volume fraction is calculated as:

```text
Phase fraction = phase voxel count / analysed specimen voxel count
```

The physical phase volume is calculated as:

```text
Physical volume =
phase voxel count × voxel size X × voxel size Y × voxel size Z
```

The phase-fraction module reports:

- Label value
- Phase name
- Voxel count
- Volume fraction
- Volume fraction percentage
- Physical volume in µm³
- Physical volume in mm³
- Directional phase-fraction profiles along X, Y and Z

Only labels explicitly selected by the user are included in the analysed
specimen volume.

---

## Pore-size analysis

CemCT uses a local-thickness method to characterise the local size of a
selected phase.

The local diameter is calculated as:

```text
Local diameter = 2 × local radius × isotropic voxel size
```

The resulting distribution is volume weighted because every selected-phase
voxel contributes one local-diameter value.

The default pore-size classes used for cementitious materials are:

| Local pore diameter | Classification |
|---:|---|
| `< 1 µm` | Small capillary pores |
| `1–10 µm` | Large capillary pores |
| `> 10 µm` | Macropores |

These limits are analytical classifications and should be interpreted
together with the physical voxel size and spatial resolution of the input
XRµCT dataset.

The pore-size module can produce:

- Selected-label binary mask
- Local-diameter map
- Pore-size-class contributions
- Differential pore-size distribution
- Cumulative pore-volume distribution
- Excel tables
- PNG summary figure

---

## Directional diffusion analysis

Finite-difference diffusion is calculated only within the directionally
through-connected region of the selected phase.

For each analysed direction, CemCT reports:

- Directional connectivity
- Effective porosity
- Formation factor
- Relative effective diffusivity
- Diffusion tortuosity
- Solver runtime
- Three-dimensional concentration field

The reported quantities are related by:

```text
Formation factor, F = D0 / Deff
```

```text
Relative effective diffusivity = Deff / D0 = 1 / F
```

```text
Diffusion tortuosity = effective porosity × formation factor
```

Only a phase that connects the two opposing boundaries in the selected
direction can contribute to steady-state directional diffusion.

---

## Random-walk transport analysis

CemCT provides random-walk-based transport analysis within a selected
connected phase.

The method tracks the displacement of random walkers constrained to move
within the accessible transport region.

The random-walk module can report:

- Accessible voxel count
- Accessible phase fraction
- Mean-square displacement
- Directional mean-square displacement
- Fitted directional slopes
- Directional random-walk tortuosity
- Fit ranges
- Coefficients of determination
- Transport mask

Random-walk and finite-difference tortuosity are complementary numerical
descriptions and should not automatically be interpreted as identical
quantities.

---

## Pore-network extraction and permeability

CemCT uses SNOW2-based pore-network extraction and OpenPNM-based
Stokes-flow simulation to estimate directional intrinsic permeability.

The permeability workflow consists of:

```text
Selected phase
→ directional connectivity extraction
→ SNOW2 pore-network extraction
→ OpenPNM network conversion
→ hydraulic model assignment
→ steady-state Stokes-flow simulation
→ intrinsic permeability calculation
```

For each direction, the module can report:

- Selected-phase voxel count
- Through-connected voxel count
- Directional connectivity
- Network pore count
- Network throat count
- Volumetric flow rate
- Pressure drop
- Simulation length
- Cross-sectional area
- Intrinsic permeability in m²
- Intrinsic permeability in Darcy
- Intrinsic permeability in mD
- Three-dimensional pressure field

Permeability results depend on segmentation quality, spatial resolution,
network-extraction parameters and the assumptions used in the pore-network
models.

---

## Software requirements

CemCT is developed in Python and currently targets Python 3.11.

Core dependencies include:

- NumPy
- SciPy
- pandas
- Matplotlib
- tifffile
- scikit-image
- openpyxl

Optional scientific dependencies include:

- PoreSpy
- OpenPNM
- PyAMG
- PyTrax
- pyevtk

The package has primarily been developed and tested on Windows. Additional
cross-platform testing will be performed before formal release.

---

## Installation from source

A dedicated Conda environment is recommended.

Open Anaconda Prompt and run:

```bat
conda create -n cemct python=3.11 -y
conda activate cemct
cd /d path\to\CemCT
python -m pip install -e ".[all]"
```

For development and automated testing, install the scientific and development
dependency groups:

```bat
python -m pip install -e ".[all,dev]"
```

If the `dev` optional dependency group is not defined in `pyproject.toml`,
install the development tools separately:

```bat
python -m pip install pytest build twine
```

---

## Installation verification

Check the installed package version:

```bat
python -c "import cemct; print(cemct.__version__)"
```

Check the installation path:

```bat
python -c "import cemct; print(cemct.__file__)"
```

Run the automated tests:

```bat
python -m pytest -v
```

Display the command-line interface:

```bat
python -m cemct --help
```

On Windows systems where organisational security policies block the generated
`cemct.exe` launcher, use:

```bat
python -m cemct
```

instead of:

```bat
cemct
```

---

## Basic Python usage

### Loading a labelled TIFF volume

```python
from pathlib import Path

from cemct.io import load_labelled_tiff


file_path = Path(
    r"C:\path\to\labelled_volume.tif"
)

volume = load_labelled_tiff(file_path)

print("Volume shape:", volume.shape)
print("Data type:", volume.dtype)
```

### Phase-fraction analysis

```python
from cemct.phase_fraction import calculate_phase_fractions


phase_fraction_results = calculate_phase_fractions(
    volume=volume,
    phases={
        1: "Unreacted particle",
        2: "Reaction product",
        3: "Pore",
        4: "Other phase",
    },
    voxel_size_um=0.7,
)

print(phase_fraction_results)
```

### Connectivity analysis

```python
from cemct.connectivity import analyse_connectivity


connectivity_results = analyse_connectivity(
    volume=volume,
    phase_label=3,
)

print(connectivity_results["summary"])
```

### Pore-size analysis

```python
from cemct.pore_size import calculate_pore_size_distribution


pore_size_results = calculate_pore_size_distribution(
    volume=volume,
    phase_label=3,
    voxel_size_um=0.7,
    number_of_bins=30,
)

print(pore_size_results["summary"])
```

### Directional diffusion analysis

```python
from cemct.diffusion import calculate_directional_diffusion


diffusion_results = calculate_directional_diffusion(
    volume=volume,
    phase_label=3,
    directions=("X", "Y", "Z"),
)

print(diffusion_results["summary"])
```

### Random-walk analysis

```python
from cemct.random_walk import calculate_connected_random_walk


random_walk_results = calculate_connected_random_walk(
    volume=volume,
    phase_label=3,
    connectivity_mode="any",
    number_of_steps=10000,
    number_of_walkers=5000,
    stride=10,
)

print(random_walk_results["summary"])
```

### Intrinsic-permeability analysis

```python
from cemct.permeability import analyse_permeability


permeability_results = analyse_permeability(
    volume=volume,
    phase_label=3,
    voxel_size_um=0.7,
    directions=("X", "Y", "Z"),
    boundary_width=3,
)

print(permeability_results["summary"])
```

---

## Demonstration notebook

The complete demonstration workflow is located at:

```text
examples/CemCT_Demonstration.ipynb
```

The notebook demonstrates:

1. Package and dependency verification
2. TIFF loading and input validation
3. User-defined multi-phase fraction analysis
4. Isolated, accessible and directional connectivity analysis
5. Local-thickness-based pore-size analysis
6. Directional finite-difference diffusion analysis
7. Random-walk-based transport analysis
8. Pore-network extraction and intrinsic-permeability analysis
9. Unified export of numerical, image and metadata results

To validate the notebook, use:

```text
Kernel → Restart Kernel and Run All Cells
```

All cells should execute from beginning to end without requiring hidden
variables from a previous session.

---

## Output files

Depending on the analysis module, CemCT can export:

- Binary TIFF masks
- Float32 scalar-field TIFF volumes
- Excel workbooks
- JSON metadata files
- NPZ numerical archives
- PNG summary figures

Output filenames are derived from the input TIFF stem and the analysis name.

Example output directory:

```text
sample_cemct_results/
```

Example output files:

```text
sample_connectivity_results.xlsx
sample_connectivity_metadata.json
sample_connectivity_through_connected_X.tif
sample_connectivity_through_connected_Y.tif
sample_connectivity_through_connected_Z.tif
sample_connectivity_summary.png
```

---

## Unified export layer

The unified export interface can be called using:

```python
from cemct.export import export_analysis_bundle


exported_paths = export_analysis_bundle(
    input_file=file_path,
    analysis_name="connectivity",
    tables={
        "Connectivity Summary": connectivity_results["summary"],
    },
    masks={
        "through_connected_X": (
            connectivity_results["masks"]["through_connected_X"]
        ),
        "through_connected_Y": (
            connectivity_results["masks"]["through_connected_Y"]
        ),
        "through_connected_Z": (
            connectivity_results["masks"]["through_connected_Z"]
        ),
    },
    metadata={
        "phase_label": 3,
        "voxel_size_um": 0.7,
    },
)

print(exported_paths)
```

---

## Automated testing

Run the complete test suite from the project root:

```bat
conda activate cemct
cd /d path\to\CemCT
python -m pytest -v
```

A release candidate should only be built when all tests pass.

The automated tests currently cover:

- TIFF input validation
- Phase-fraction calculations
- Strict 6-connectivity
- Exclusion of edge-only and corner-only contacts
- Local-diameter conversion
- Pore-size distributions
- Diffusion calculations
- Random-walk analysis
- Pore-network and permeability calculations
- Unified export functions
- Command-line interface behaviour

---

## Scientific scope

CemCT is a post-segmentation analysis package.

The following operations are outside the scope of the current version:

- XRµCT image acquisition
- Tomographic reconstruction
- Ring-artefact correction
- Beam-hardening correction
- Image registration
- Image denoising
- Image segmentation
- Manual correction of segmented labels

The quality and scientific interpretation of the outputs depend on:

- Input image quality
- Segmentation accuracy
- Physical voxel size
- Effective spatial resolution
- Representativeness of the selected volume
- Boundary conditions
- Model parameters
- Numerical convergence

Users are responsible for confirming that the selected analysis methods and
parameters are appropriate for their dataset and scientific question.

---

## Validation status

The current release candidate has been tested using:

- Controlled synthetic labelled 3D datasets
- Known phase-fraction targets
- Analytical and topological unit tests
- Strict face-connectivity test geometries
- Directional integration tests
- End-to-end demonstration-notebook execution
- Clean execution of the six scientific analysis modules
- Automated verification of exported files

Additional validation using experimental cementitious-material datasets,
resolution-sensitivity analysis, parameter-sensitivity analysis and
cross-software comparison is in progress.

---

## Known limitations

The current release candidate has the following limitations:

- Image reconstruction and segmentation are not included.
- Local-thickness analysis assumes an isotropic voxel size.
- The interpretation of sub-resolution pores is limited by the effective
  spatial resolution of the input data.
- Directional transport calculations require a through-connected selected
  phase.
- Random-walk results depend on walker count, step count and fitting range.
- Pore-network extraction depends on SNOW2 parameters.
- Permeability results depend on network geometry and hydraulic-model
  assumptions.
- Very large TIFF volumes may require substantial memory and computation time.
- The package has not yet been systematically validated across all operating
  systems.

---

## Recommended reporting information

Scientific publications using CemCT should report:

- CemCT version
- Input volume dimensions
- NumPy coordinate convention
- Physical voxel dimensions
- Selected phase label
- Connectivity definition
- Pore-size method and binning
- Diffusion solver and convergence tolerance
- Random-walk parameters
- SNOW2 extraction parameters
- Permeability boundary conditions
- Any preprocessing performed before analysis

---

## Project structure

```text
CemCT/
├── pyproject.toml
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── src/
│   └── cemct/
│       ├── __init__.py
│       ├── __main__.py
│       ├── cli.py
│       ├── io.py
│       ├── phase_fraction.py
│       ├── connectivity.py
│       ├── pore_size.py
│       ├── diffusion.py
│       ├── random_walk.py
│       ├── permeability.py
│       └── export.py
├── tests/
└── examples/
    └── CemCT_Demonstration.ipynb
```

---

## Citation

Citation metadata will be supplied through:

```text
CITATION.cff
```

following confirmation of the final software name, version, authorship,
repository URL and DOI.

Until the first formal release, please cite the associated X-SeeO2 Task 8.1
documentation where appropriate.

---

## Development team

- Zixian Su
- Timothy L. Burnett
- Philip J. Withers

### Contact

For development enquiries, contact:

```text
Zixian Su
School of Materials Science and Engineering
Central South University
Changsha, China
```

Formal contact details and corresponding-author information will be added
following confirmation by the development team.

---

## Project acknowledgement

CemCT was developed through work associated with X-SeeO2 Task 8.1 and the
National X-ray Computed Tomography Facility at The University of Manchester.

The formal funding statement, project grant information and institutional
acknowledgements will be added before public release.

---

## License

The intended open-source license is the BSD 3-Clause License.

The final software ownership, copyright statement, licensing terms and
public-release conditions must be confirmed before formal distribution.

Until these matters are confirmed, this release candidate should be treated as
an unreleased working version and should not be redistributed.
