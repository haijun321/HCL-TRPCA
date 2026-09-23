# HCL-TRPCA v1.0.2 — submission reproducibility update

Release v1.0.2 synchronizes the public repository with the submission-stage manuscript and Supplementary Material.

## Added

- development-calibrated nonhierarchical attribution-control materials and frozen calibration coefficients;
- recovery-floor and mask-source mechanism-control manuscript artifacts;
- Pavia University rank-budget V2.1 DIAG package;
- Pavia fixed-mask recovery trajectories and publication figure;
- archived paired confirmation statistics and a transparent audit utility;
- final manuscript-oriented figures associated with the added diagnostics.

## Documentation / portability

- added result-to-manuscript mapping;
- removed PILOT/REPRO/SMOKE/duplicate/partial staging artifacts from the public update package;
- sanitized local absolute Pavia data paths. Public code now uses `HCL_TRPCA_DATA_ROOT` or `<repository>/data`.

## Unchanged

- HCL-TRPCA core method;
- frozen primary HCL parameter configuration;
- development and confirmation seed sets;
- primary confirmation observations and reported outcomes.

## Scope note

The exact legacy runner used for the development-stage recovery-mechanism control and the exact legacy confirmation-statistics postprocessor were not present in the staging archive. Their frozen numerical outputs are released. The package includes a separately identified reproduction/audit utility for the archived paired statistics; it is not represented as the original legacy postprocessor.
