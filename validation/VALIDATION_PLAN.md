\# CemCT Validation Plan



\## 1. Synthetic structural validation



\- Phase fractions with known voxel counts

\- Strict 6-connectivity

\- X-, Y- and Z-directional through-connectivity

\- Isolated and boundary-accessible regions

\- Local-diameter calculation for idealised pores



\## 2. Analytical transport validation



\- Diffusion through a straight channel

\- Formation factor and relative effective diffusivity

\- Random-walk behaviour in an unobstructed domain

\- Intrinsic permeability of an ideal cylindrical or square channel

\- Comparison with analytical solutions



\## 3. Numerical sensitivity



\- Voxel-size sensitivity

\- Image-resolution sensitivity

\- Local-thickness parameter sensitivity

\- SNOW2 parameter sensitivity

\- Random-walk walker-number and step-number sensitivity

\- Solver-tolerance sensitivity



\## 4. Computational performance



Record for each module:



\- Input volume dimensions

\- Number of voxels

\- Runtime

\- Peak memory usage

\- Python version

\- Operating system

\- Processor

\- Available RAM



Recommended volumes:



\- 50 × 50 × 50 voxels

\- 100 × 100 × 100 voxels

\- 200 × 200 × 200 voxels



\## 5. Experimental demonstration



Use at least one segmented experimental XRµCT dataset of a cementitious

material, subject to data ownership and release approval.



Compare selected CemCT outputs with:



\- Existing published measurements

\- Established open-source software

\- Commercial software where an appropriate reference result is available



\## 6. Acceptance criteria



\- Automated unit and topology tests pass

\- Analytical deviations are quantified

\- Numerical trends are physically reasonable

\- Results are reproducible using fixed parameters

\- Demonstration notebook runs from a clean environment

\- All input assumptions and limitations are documented

