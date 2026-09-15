# Experimental CT data for the CemCT manuscript

These files were supplied by Zixian Su for the 14-day sodium-sulfate-activated slag paste example in CemCT_SoftwareX_V3.5.

| File | Content | Type |
|---|---|---|
| raw_CT.tif | Reconstructed grayscale volume, not projection data | uint16 |
| NasSO4_slag_raw_BC_Seg.tif | Segmented label volume | uint8 |

Both arrays have shape (200, 200, 200) in Z, Y, X order. Isotropic voxel size is 0.7 um, giving a 140 x 140 x 140 um field of view. Spatial sampling is not a measurement of imaging resolution.

Labels: 0 = unreacted slag; 1 = reaction products; 2 = resolved capillary pores. Counts are 2,316,476; 4,280,648; and 1,402,876, respectively. All three labels are included in the phase-fraction denominator. Transport analyses select label 2.

The original files contained approximately 0.529233 um spatial calibration. Following the data owner's confirmation, only X/Y resolution, Z spacing and the unit text were corrected to 0.7 um. Every voxel value, dtype and array shape was verified unchanged. No resampling or relabelling was performed. The original filenames are retained. manifest.json records SHA-256 checksums of these corrected copies.

Run examples/reproduce_manuscript.ipynb using examples/manuscript_parameters.json. The grayscale data are included for inspection; calculations use the segmented volume. The final reproduction configuration is explicitly recorded in the parameter file and notebook, including a full-interval (0.0–1.0) through-origin random-walk fit. The grayscale data, voxel counts and numerical comparison outputs make the reproduction inspectable; matching rounded values alone does not prove undocumented historical settings.

## Data ownership and licensing

The experimental data belong to the research team. The data owner reports no additional collaborator requirements. A specific open-data license has not yet been selected. This update makes no new grant of data permissions and does not label the TIFFs or derived experimental results as CC BY or another open-data license.

The software retains its separate BSD-3-Clause license. No dataset DOI has been supplied. Record the selected data license and actual archive identifier when available.
