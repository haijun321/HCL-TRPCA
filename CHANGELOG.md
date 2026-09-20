# Changelog

All notable changes to the public HCL-TRPCA software archive are documented here.

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
