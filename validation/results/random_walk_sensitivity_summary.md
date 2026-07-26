# CemCT Random-Walk Sensitivity Analysis

## Configuration

- Domain shape: `(64, 64, 64)`
- Selected phase label: `3`
- Expected open-domain tortuosity: `1.000`
- Deterministic seed: `True`
- MSD fit interval: `0.05` to `0.35`

## Walker-count reference case

- Number of walkers: `5000`
- Number of steps: `600`
- Mean tortuosity: `0.984498`
- Mean R-squared: `0.999636`
- Analytical error: `1.550%`
- Change relative to preceding case: `0.568%`
- Status: `PASS`

## Step-count reference case

- Number of walkers: `3000`
- Number of steps: `1000`
- Mean tortuosity: `0.961492`
- Mean R-squared: `0.999623`
- Analytical error: `3.851%`
- Change relative to preceding case: `1.811%`
- Status: `PASS`

## Overall assessment

`PASS`

A review status does not necessarily indicate a software error. It indicates
that the selected numerical parameters have not yet satisfied the predefined
convergence or analytical-tolerance criteria.
