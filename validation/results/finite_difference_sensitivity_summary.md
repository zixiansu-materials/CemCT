# Finite-Difference Diffusion Solver Sensitivity

**Validation status:** PASS

## Purpose

This validation evaluates the numerical sensitivity of CemCT's X-direction finite-difference diffusion analysis to solver tolerance and maximum iteration count.

## Validation domain

- Volume shape (Z, Y, X): `(40, 40, 56)`
- Selected phase label: `3`
- Selected phase fraction: `0.438984`
- Direction: `X`
- Random seed: `20260726`
- Connectivity definition: strict three-dimensional 6-connectivity

## Reference configuration

- Solver tolerance: `1.0e-10`
- Maximum iterations: `5000`
- Reference diffusion tortuosity: `2.560663885`
- Reference formation factor: `6.002864501`
- Reference relative effective diffusivity: `0.1665871352`

## Tolerance sensitivity

| Solver tolerance | Diffusion tortuosity | Formation factor | Deff / D0 | Tortuosity difference (%) | Runtime (s) |
|---:|---:|---:|---:|---:|---:|
| 1.0e-04 | 2.2912828 | 5.3713648 | 0.18617242 | 10.52 | 1.3740 |
| 1.0e-06 | 2.5595531 | 6.0002606 | 0.16665943 | 0.043378 | 0.2063 |
| 1.0e-08 | 2.5606541 | 6.0028415 | 0.16658777 | 0.000383239 | 0.3040 |
| 1.0e-10 | 2.5606639 | 6.0028645 | 0.16658714 | 0 | 0.3539 |

## Maximum-iteration sensitivity

| Maximum iterations | Diffusion tortuosity | Formation factor | Deff / D0 | Tortuosity difference (%) | Runtime (s) |
|---:|---:|---:|---:|---:|---:|
| 250 | 2.5606639 | 6.0028645 | 0.16658714 | 0 | 0.3094 |
| 500 | 2.5606639 | 6.0028645 | 0.16658714 | 0 | 0.3676 |
| 1000 | 2.5606639 | 6.0028645 | 0.16658714 | 0 | 0.3393 |
| 5000 | 2.5606639 | 6.0028645 | 0.16658714 | 0 | 0.3539 |

## Numerical consistency

- Successful configurations: `1`
- Configurations within 1% of the reference: `7/8`
- Maximum tortuosity difference from the reference: `10.52%`
- Maximum reciprocal consistency error `abs(F × Deff/D0 - 1)`: `0.000000e+00`

## Interpretation

Configurations with a diffusion-tortuosity difference below 1.0% are classified as numerically converged relative to the highest-accuracy reference configuration.

Runtime measurements are wall-clock measurements and may vary with operating system, processor load and available memory.

## Output files

- Numerical results: `finite_difference_sensitivity_results.csv`
- Summary figure: `finite_difference_sensitivity_analysis.png`
- Validation volume: `finite_difference_sensitivity_domain.tif`
