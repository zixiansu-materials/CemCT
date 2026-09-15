# Release metadata review

## Completed

- Software authors are consistently Zixian Su, Timothy L. Burnett, Miguel A. G. Aranda and Philip J. Withers in README.md, CITATION.cff and pyproject.toml. Zixian Su remains the maintainer.
- Software version remains 1.0.0. Manuscript V3.6 is a document revision only. Public tag and DOI status are not asserted.
- Software license identifier remains BSD-3-Clause, with LICENSE unchanged.
- Citation URL placeholders are replaced by the repository address supplied with the manuscript. The CFF email is a single address matching the manuscript support email.
- The two experimental TIFF volumes in data/experimental/ and the team-generated accompanying results in results/manuscript/ are licensed under Creative Commons Attribution 4.0 International (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/. Software remains licensed under BSD-3-Clause.
- README and CHANGELOG reflect completed experimental reproduction. No source algorithms, TIFFs, configuration, notebook or verified numerical outputs were modified.

## Remaining publication metadata

Verify repository visibility and release tag; add real software/data archive DOI values and article citation when available. No public release or deposit was performed in this review.

## Checks

TOML and CFF YAML parsing, author/version/license consistency, ZIP integrity and byte-for-byte preservation of all pre-existing files outside the edited metadata files were checked. This is not a full external CFF schema validation. Scientific tests were not rerun because source, configuration and results are unchanged.
