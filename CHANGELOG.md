\# Changelog



All notable changes to CemCT will be documented in this file.



The format is based on Keep a Changelog, and the project intends to follow

Semantic Versioning after formal public release.



\## \[Unreleased]



\### Added



\- Source-layout Python package under `src/cemct`

\- Input validation for integer-labelled 3D TIFF volumes

\- Binary-mask TIFF export using the `0/255` convention

\- User-defined multi-phase fraction analysis

\- Physical phase-volume calculation

\- Directional phase-fraction profiles along X, Y and Z

\- Strict 6-connectivity analysis

\- Isolated pore-region extraction

\- Boundary-accessible pore-region extraction

\- Boundary-accessible but non-through-connected region extraction

\- X-, Y- and Z-directionally through-connected region extraction

\- Directional connectivity calculation

\- Local-thickness-based pore-size analysis

\- Local-diameter TIFF generation

\- Pore-size-class contribution analysis

\- Differential pore-size distribution

\- Cumulative pore-volume analysis

\- Directional finite-difference diffusion analysis

\- Formation-factor calculation

\- Relative-effective-diffusivity calculation

\- Diffusion-tortuosity calculation

\- Three-dimensional concentration-field output

\- Random-walk-based transport analysis

\- Mean-square-displacement analysis

\- Directional random-walk tortuosity calculation

\- SNOW2 pore-network extraction

\- OpenPNM-based steady-state Stokes-flow simulation

\- Directional intrinsic-permeability estimation

\- Three-dimensional pressure-field output

\- Unified TIFF, Excel, JSON, NPZ and PNG export layer

\- Command-line interface

\- Automated unit and integration tests

\- Complete demonstration notebook

\- Controlled synthetic four-phase demonstration dataset

\- Project README and preliminary citation metadata



\### Changed



\- Migrated the original Jupyter-based workflow into reusable Python modules

\- Removed `input()`, `print()` and `plt.show()` from core scientific functions

\- Separated scientific calculations from interactive input and file export

\- Standardised the NumPy coordinate convention as `(Z, Y, X)`

\- Standardised directional naming as `X`, `Y` and `Z`

\- Standardised output filenames using the input-file stem and analysis name

\- Standardised binary output masks using unsigned 8-bit `0/255` values



\### Fixed



\- Corrected binary TIFF handling for `0/1`, Boolean and `0/255` conventions

\- Enforced face-only connectivity in three-dimensional component analysis

\- Excluded edge-only and corner-only voxel contacts from connectivity

\- Corrected local-radius to physical-diameter conversion

\- Improved directional status checking for diffusion calculations

\- Improved compatibility with current PoreSpy, OpenPNM and PyAMG versions

\- Improved metadata serialisation for NumPy scalar and array types

\- Improved non-interactive figure export behaviour



\### Validation



\- Added analytical phase-fraction tests

\- Added face-, edge- and corner-connectivity tests

\- Added local-diameter conversion tests

\- Added pore-size-distribution tests

\- Added directional diffusion tests

\- Added random-walk tests

\- Added permeability tests

\- Added unified export tests

\- Added command-line interface tests

\- Completed end-to-end execution of `CemCT\_Demonstration.ipynb`



\## \[1.0.0] - To be released



\### Planned



\- First formal public release of CemCT

\- Public source-code repository

\- Archived software release

\- Software DOI

\- Final citation metadata

\- Cross-platform clean-install validation

\- Experimental cementitious-material demonstration dataset

\- SoftwareX software publication

