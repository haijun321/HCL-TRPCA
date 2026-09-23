# Changelog

All notable changes to the public HCL-TRPCA software archive are documented here.

## v1.0.2 — 2026-09-23

### Added
- Development-calibrated nonhierarchical attribution-control materials.
- Frozen Calibrated-Flat calibration coefficients.
- Recovery-floor and mask-source mechanism-control artifacts.
- Pavia rank-budget V2.1 DIAG results.
- Pavia fixed-mask recovery trajectories.
- Paired confirmation-statistics audit utility.
- Final manuscript-oriented diagnostic figures.

### Documentation / portability
- Added result-to-manuscript mapping.
- Added v1.0.2 portability notes.
- Removed author-machine Pavia data-path dependencies from the public update.
- Public real-data code now resolves data through
  `HCL_TRPCA_DATA_ROOT` or `<repository>/data`.

### Unchanged
- HCL-TRPCA core method.
- Frozen primary HCL configuration.
- Development and confirmation seed sets.
- Primary confirmation observations and reported outcomes.

### Reproducibility note
The exact legacy development runner for the recovery-mechanism control
and the exact legacy paired-statistics postprocessor were not present in
the final staging archive. Their frozen numerical outputs are archived.
The v1.0.2 repository provides a separately identified audit utility for
the archived paired statistics and does not represent it as the original
legacy source.

## v1.0.1 - 2026-09-20

### SBI ground-truth semantics clarification

- Clarified that the distributed SBI `groundtruth/gt*.png` files audited for Board and CAVIAR1 are foreground segmentation masks rather than clean-background reconstruction references.
- Changed `run_sbi_standardized.m` so the presence of `GT` or `BgtUint8` alone can no longer trigger PSNR, SSIM, or AGE evaluation.
- Required both `GroundTruthType = 'clean_background_reference'` and `ValidReconstructionReference = true` before a separately supplied clean background can be used for reconstruction metrics.
- Added explicit annotation, alignment, clean-reference, and metric-status fields to newly generated SBI outputs.
- Replaced the ambiguous legacy `hasGT` field in `results/sbi/SBI_review_small.mat` with semantically explicit fields while preserving all archived ranks, structural diagnostics, and stopping values.
- Added `SBI_REFERENCE_SEMANTICS.md` and updated the README, SBI protocol, result index, reproduction entry point, citation metadata, and checksums.
- Added `audit_sbi_release_semantics.m`, a raw-data-free check of the archived semantic flags, `NaN` metric policy, and CSV/MAT diagnostic consistency.

This patch changes evaluation semantics and documentation only. It does not change the HCL-TRPCA algorithm, optimization procedure, frozen parameters, or archived structural diagnostic values.

## v1.0.0

- Initial public software and reproducibility archive.
