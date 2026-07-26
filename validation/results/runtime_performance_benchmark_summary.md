# CemCT Runtime and Performance Benchmark

## Purpose

This benchmark evaluates the wall-clock runtime and resident-memory requirements of the six principal CemCT analysis modules.

## Benchmark configuration

- Volume sizes: `[40, 60, 80]`
- Repeats per case: `3`
- Voxel size: `0.7 um`
- Transport direction: `X`
- Random-walk steps: `2000`
- Random-walk walkers: `1000`
- Solver tolerance: `1.0e-10`
- Maximum solver iterations: `5000`

## Computing environment

- Operating system: `Windows-10-10.0.26200-SP0`
- Python: `3.11.15`
- Logical CPU count: `16`
- Total system memory: `31.52 GB`
- CemCT version: `1.0.0`

## Completion status

- Successful benchmark runs: `54`
- Failed benchmark runs: `0`

## Aggregated results

| module                      |   volume_size |   voxel_count |   median_wall_time_seconds |   median_additional_peak_rss_mb |   successful_runs |
|:----------------------------|--------------:|--------------:|---------------------------:|--------------------------------:|------------------:|
| connectivity                |            40 |         64000 |                   0.00688  |                          0.8164 |                 3 |
| connectivity                |            60 |        216000 |                   0.006406 |                          1.145  |                 3 |
| connectivity                |            80 |        512000 |                   0.01464  |                          4.652  |                 3 |
| finite_difference_diffusion |            40 |         64000 |                   0.4227   |                         20.4    |                 3 |
| finite_difference_diffusion |            60 |        216000 |                   0.8337   |                         64.81   |                 3 |
| finite_difference_diffusion |            80 |        512000 |                   1.893    |                        183      |                 3 |
| permeability                |            40 |         64000 |                   0.827    |                          3.449  |                 3 |
| permeability                |            60 |        216000 |                   0.8475   |                         12.16   |                 3 |
| permeability                |            80 |        512000 |                   1.336    |                         35.02   |                 3 |
| phase_fraction              |            40 |         64000 |                   0.003693 |                          0.832  |                 3 |
| phase_fraction              |            60 |        216000 |                   0.003112 |                          0.4805 |                 3 |
| phase_fraction              |            80 |        512000 |                   0.005815 |                          1.086  |                 3 |
| pore_size                   |            40 |         64000 |                   0.6788   |                          1.184  |                 3 |
| pore_size                   |            60 |        216000 |                   0.5564   |                          3.816  |                 3 |
| pore_size                   |            80 |        512000 |                   1.026    |                          9.438  |                 3 |
| random_walk                 |            40 |         64000 |                   0.5813   |                          5.441  |                 3 |
| random_walk                 |            60 |        216000 |                   0.2178   |                          5.812  |                 3 |
| random_walk                 |            80 |        512000 |                   0.2184   |                         12.63   |                 3 |

## Failed calculations

No module failures were recorded.

## Interpretation

Runtime and memory values are implementation- and hardware-specific. They should therefore be interpreted as performance measurements for the computing environment reported above, rather than universal execution times.

Finite-difference diffusion and pore-network permeability were calculated only in the X direction to keep the benchmark computationally practical.

The reported additional peak RSS is the maximum observed increase in process resident memory relative to the baseline immediately before each module call.
