# Analytical pore-size validation

## Objective

Validate the CemCT local-thickness-based pore-diameter calculation using straight square channels of known width.

## Analytical reference

Expected diameter = channel width × voxel size.

## Acceptance criterion

Relative median-diameter error <= 10.0%.

## Results

|   channel_width_voxels |   voxel_size_um |   expected_diameter_um |   measured_median_diameter_um |   measured_mean_diameter_um |   absolute_error_um |   relative_error_percent |   evaluated_voxel_count |   acceptance_limit_percent | status   |
|-----------------------:|----------------:|-----------------------:|------------------------------:|----------------------------:|--------------------:|-------------------------:|------------------------:|---------------------------:|:---------|
|                      6 |               1 |                      6 |                             6 |                     6       |         0           |              0           |                    1152 |                         10 | passed   |
|                     10 |               1 |                     10 |                            10 |                     9.67217 |         0           |              0           |                    3200 |                         10 | passed   |
|                     14 |               1 |                     14 |                            14 |                    13.4512  |         9.53674e-07 |              6.81196e-06 |                    6272 |                         10 | passed   |

## Overall assessment

- All configurations passed: `True`
- Maximum relative error: `0.000007%`
