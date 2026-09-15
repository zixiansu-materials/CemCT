# Changelog

## Unreleased preparation for v1.0.0

The packaged software version remains 1.0.0. This section records preparation work;
it is not a declaration that a public release or DOI deposit has occurred.
Manuscript V3.6 is a separate document revision.

### Included software

- Source-layout package with six analysis modules and a TIFF-inspection CLI.
- Strict face connectivity, local-thickness pore size, finite-difference diffusion,
  random-walk transport and pore-network permeability.
- TIFF, Excel, JSON, NPZ and PNG export interfaces; pressure results are network-node arrays.
- Automated tests, analytical geometries, numerical sensitivity scripts and benchmark records.

### Experimental reproduction added on 2026-09-14

- Grayscale and segmented 200³ TIFF volumes, with author-confirmed 0.7 um calibration.
- Data descriptions, label counts and SHA-256 checksums.
- Reproduction notebook and explicit numerical configuration.
- Full-interval through-origin random-walk fitting (500 saved points, steps 0–4990).
- Completed six-module experimental run and comparison against rounded manuscript values.
- Forty local automated tests passed on the source code included in this snapshot.

### Metadata review on 2026-09-15

- Unified version and release-status wording without changing the software version.
- Replaced citation URL placeholders with the repository address already supplied in the manuscript.
- Corrected the CFF email field to a single address, using the manuscript support contact.
- Unified the software authors as Zixian Su, Timothy L. Burnett, Miguel A. G. Aranda and Philip J. Withers, as confirmed by the user.
- Preserved the BSD-3-Clause license text and the 1.0.0 software version.
- Recorded team-owned experimental data with no additional collaborator requirements; no specific open-data license has been selected.

### Remaining release work

- Data license confirmed as CC BY 4.0; see DATA_LICENSE.md.
- Verify repository visibility and the intended public release tag.
- Create a versioned archive and record the real DOI and release date.
- Add a published article citation when available.
