\# CemCT Validation Summary



\## 1. Purpose



This document summarises the verification, analytical validation, sensitivity analysis, and automated software testing completed for CemCT version 1.0.0.



CemCT is a post-segmentation software package for quantitative analysis of binary and integer-labelled three-dimensional XRmicroCT datasets of cementitious materials.



The validation programme evaluates:



\- user-defined phase-fraction calculations;

\- strict face-connected pore classification;

\- local-thickness-based pore-size measurements;

\- finite-difference diffusion calculations;

\- random-walk transport calculations;

\- pore-network-based intrinsic permeability calculations;

\- numerical sensitivity to selected solver and model parameters; and

\- software installation, import, testing, and distribution building.



All three-dimensional NumPy arrays use the following coordinate convention:



&#x20;   Array axis 0 = Z direction

&#x20;   Array axis 1 = Y direction

&#x20;   Array axis 2 = X direction



Connectivity calculations use strict 6-connectivity:



&#x20;   Face connection   = included

&#x20;   Edge connection   = excluded

&#x20;   Corner connection = excluded



\---



\## 2. Validation overview



| Validation item | Reference or expected behaviour | Principal metric | Status |

|---|---|---|---|

| Automated unit tests | All tests pass | Test pass rate | Passed |

| Synthetic structural validation | Prescribed phase fractions and pore connectivity | Fraction and topology agreement | Passed |

| Analytical pore-size validation | Known straight-channel diameter | Diameter error | Passed |

| Analytical diffusion validation | Fully open or straight transport domain | Tortuosity and relative diffusivity | Passed |

| Analytical random-walk validation | Linear MSD and tortuosity close to unity | Tortuosity and R-squared | Passed |

| Analytical permeability validation | Known straight-channel behaviour | Relative permeability error | Passed |

| Random-walk sensitivity analysis | Stable result with increasing walkers and steps | Tortuosity convergence | Passed |

| Finite-difference sensitivity analysis | Stable result with tighter solver tolerance | Formation-factor convergence | Passed |

| Permeability parameter sensitivity | Stable result across justified extraction parameters | Permeability sensitivity | Passed |

| Clean-install validation | Wheel imports outside the source directory | Installed package path | Passed |

| GitHub Actions | Tests, package build, and metadata checks succeed | Workflow status | Passed |



\---



\## 3. Automated software testing



The automated test suite covers:



\- TIFF input validation;

\- binary-mask export;

\- multi-phase fraction calculations;

\- strict face-connected pore connectivity;

\- pore-size analysis;

\- finite-difference diffusion analysis;

\- random-walk analysis;

\- permeability analysis;

\- unified output export; and

\- command-line interface behaviour.



The current local test result is:



&#x20;   40 tests passed



The GitHub Actions workflow additionally performs:



\- package installation;

\- dependency checking;

\- CemCT import verification;

\- automated testing;

\- wheel and source-distribution building; and

\- distribution-metadata checking.



The latest GitHub Actions runs completed successfully.



\---



\## 4. Synthetic structural validation



Controlled synthetic labelled volumes were used to verify phase fractions and pore topology.



The structural validation included:



\- prescribed voxel counts for selected phases;

\- known phase-fraction values;

\- isolated pore regions;

\- boundary-accessible pore regions;

\- X-direction through-connected regions;

\- Y-direction through-connected regions;

\- Z-direction through-connected regions;

\- edge-only contact; and

\- corner-only contact.



Strict 6-connectivity correctly included face-connected voxels while excluding edge-only and corner-only contacts.



The phase-fraction calculations reproduced the prescribed voxel fractions within floating-point precision.



Status: Passed.



\---



\## 5. Analytical pore-size validation



Straight channels with known voxel widths were generated and analysed using the CemCT local-thickness pore-size method.



The analytical diameter was calculated from the prescribed channel geometry and voxel size. The measured local diameter was compared with the expected physical diameter.



Validation cases included multiple channel widths to confirm that the conversion from local radius to physical diameter was applied consistently.



The comparison uses:



&#x20;   Expected diameter = channel width in voxels multiplied by voxel size



and:



&#x20;   Local diameter = 2 multiplied by local radius multiplied by voxel size



The numerical results are stored in:



&#x20;   validation/results/analytical\_pore\_size\_validation\_results.csv



The summary is stored in:



&#x20;   validation/results/analytical\_pore\_size\_validation\_summary.md



The validation figure is stored in:



&#x20;   validation/results/analytical\_pore\_size\_validation.png



Status: Passed.



\---



\## 6. Analytical finite-difference diffusion validation



The finite-difference diffusion implementation was evaluated using open and straight transport domains with known behaviour.



For a straight, unobstructed transport path, the expected behaviour is:



\- a monotonic concentration gradient;

\- relative effective diffusivity close to the analytical reference;

\- formation factor consistent with the domain porosity; and

\- diffusion tortuosity close to the analytical value.



The implementation correctly produced the expected directional transport behaviour.



The validation also confirmed that diffusion is calculated only within directionally through-connected voxels of the selected phase.



Status: Passed.



\---



\## 7. Analytical random-walk validation



Random-walk analysis was evaluated using fully open and straight-channel domains.



The validation used deterministic random seeds to make repeated calculations reproducible.



The principal acceptance criteria were:



\- mean-square displacement increases approximately linearly with time;

\- the fitted MSD slope is positive;

\- the coefficient of determination is close to unity;

\- the fully open domain produces tortuosity close to one; and

\- repeated calculations with the same deterministic seed produce the same result.



The numerical results and diagnostic plots confirmed the expected behaviour.



Status: Passed.



\---



\## 8. Analytical permeability validation



The pore-network permeability implementation was evaluated using controlled straight-channel domains.



The validation examined:



\- successful extraction of the connected pore network;

\- correct identification of inlet and outlet pores;

\- monotonic pressure variation through the network;

\- agreement between inlet and outlet flow-rate magnitudes;

\- correct physical-unit conversion; and

\- agreement with the analytical permeability reference within the defined acceptance tolerance.



Permeability is reported in:



\- square metres;

\- Darcy; and

\- millidarcy.



Status: Passed.



\---



\## 9. Random-walk sensitivity analysis



Random-walk sensitivity was evaluated by varying parameters such as:



\- number of walkers;

\- number of steps;

\- sampling stride;

\- fitting interval; and

\- deterministic seed.



The analysis showed that very small walker populations or short trajectories may produce noisier MSD curves and less stable tortuosity estimates.



Results became stable when sufficient walkers and steps were used. Deterministic seeds enabled direct comparison between parameter configurations.



The selected default parameters represent a balance between convergence and computational cost.



Status: Passed.



\---



\## 10. Finite-difference sensitivity analysis



Finite-difference sensitivity was evaluated by varying:



\- solver tolerance; and

\- maximum number of iterations.



The tested solver tolerances included:



| Solver tolerance | Formation factor | Diffusion tortuosity | Difference from reference |

|---:|---:|---:|---:|

| 1e-4 | 5.371365 | 2.291283 | 10.519972 percent |

| 1e-6 | 6.000261 | 2.559553 | 0.043378 percent |

| 1e-8 | 6.002841 | 2.560654 | 0.000383 percent |

| 1e-10 | 6.002865 | 2.560664 | 0.000000 percent |



The tolerance of 1e-4 was insufficiently converged and produced a noticeable inlet-outlet rate mismatch.



The results at 1e-8 and 1e-10 were effectively converged.



At a solver tolerance of 1e-10, increasing the maximum number of iterations from 250 to 5000 did not materially change the result. This indicates that the solution converged before reaching the iteration limit.



The adopted high-accuracy solver setting is therefore justified.



Generated files include:



&#x20;   validation/results/finite\_difference\_sensitivity\_results.csv

&#x20;   validation/results/finite\_difference\_sensitivity\_analysis.png

&#x20;   validation/results/finite\_difference\_sensitivity\_summary.md

&#x20;   validation/results/finite\_difference\_sensitivity\_domain.tif



Status: Passed.



\---



\## 11. Permeability parameter sensitivity analysis



Permeability sensitivity was evaluated for parameters associated with pore-network extraction and flow simulation.



The investigated parameters included:



\- boundary width;

\- sigma;

\- maximum-filter radius;

\- network extraction stability;

\- number of extracted pores;

\- number of extracted throats; and

\- calculated intrinsic permeability.



The analysis was used to identify whether modest parameter changes caused disproportionate changes in the extracted network or calculated permeability.



The selected default parameters were retained because they provided stable network extraction and physically consistent permeability results for the validation geometry.



Generated numerical tables and figures are retained in the validation results directory.



Status: Passed.



\---



\## 12. Reproducibility and installation validation



CemCT was built as both a wheel and a source distribution.



The wheel was installed into an independent Conda environment. Validation was performed outside the project source directory to prevent accidental imports from the editable source tree.



The clean-install checks confirmed that:



\- CemCT imports from the environment site-packages directory;

\- the installed version is reported correctly;

\- the core phase-fraction calculation runs successfully;

\- package dependencies are consistent;

\- the command-line interface is available; and

\- the built distribution metadata passes Twine validation.



Status: Passed.



\---



\## 13. Validation evidence



Validation evidence is organised under:



&#x20;   validation/



The directory contains:



\- validation scripts;

\- numerical CSV outputs;

\- Markdown summaries;

\- TIFF validation domains; and

\- PNG validation figures.



Large generated data files do not need to be included in the Python package distribution. Selected compact validation tables and figures may be retained in the source repository and cited in the manuscript.



The full validation should be reproducible by running the corresponding scripts from the project root using the CemCT development environment.



\---



\## 14. Current limitations



The completed validation demonstrates software correctness for controlled synthetic and analytical cases. It does not establish universal accuracy for every experimental XRmicroCT dataset.



Important limitations include:



\- results depend on segmentation quality;

\- pore sizes below the image resolution cannot be resolved reliably;

\- local-thickness pore size is a voxel-based local diameter, not an explicit pore-body diameter;

\- pore-network permeability depends on network-extraction and geometrical-model assumptions;

\- random-walk estimates depend on trajectory length, walker count, and fitting interval;

\- finite-difference transport requires a directionally through-connected phase;

\- experimental image artefacts may affect all calculated quantities; and

\- validation using a real segmented cementitious-material dataset is still required.



These limitations will be stated explicitly in the software documentation and SoftwareX manuscript.



\---



\## 15. Overall conclusion



The current validation programme demonstrates that CemCT:



\- reproduces prescribed phase fractions;

\- applies strict face-connected pore classification correctly;

\- recovers known pore dimensions in controlled geometries;

\- produces analytically consistent finite-difference diffusion behaviour;

\- produces linear and reproducible random-walk MSD behaviour;

\- provides stable pore-network permeability calculations;

\- responds predictably to changes in solver and model parameters;

\- passes the complete automated test suite;

\- builds valid Python distributions; and

\- installs and imports successfully in a clean environment.



The software is therefore suitable for release-candidate testing and preparation of the SoftwareX manuscript.



Before formal public release, the remaining validation work is limited to:



1\. runtime and memory benchmarking;

2\. demonstration using at least one real segmented XRmicroCT dataset; and

3\. incorporation of the validation results into the manuscript.



No additional validation branches are planned unless a current validation result identifies a specific unresolved issue.

