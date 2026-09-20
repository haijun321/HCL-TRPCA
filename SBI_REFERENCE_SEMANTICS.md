# SBI reference-semantics clarification

## Scope

The SBI experiments in this repository are used as real-video structural and operational diagnostics of HCL-TRPCA. They are not used as reference-based background-recovery benchmarks.

This note records the semantic distinction introduced in release v1.0.1. It does not modify the HCL-TRPCA algorithm, its frozen parameters, or the archived structural diagnostic values.

## Audit conclusion

The distributed `groundtruth/gt*.png` files examined for the Board and CAVIAR1 sequences have frame/index and spatial correspondence with the associated video data, but their image characteristics and dataset role identify them as foreground segmentation masks. They are not clean-background reference frames.

Alignment and semantics answer different questions:

- spatial and frame-index alignment indicate where an annotation belongs;
- clean-reference semantics determine whether a file can support reconstruction-error metrics.

An aligned foreground mask cannot be substituted for a clean background when computing PSNR, SSIM, MS-SSIM, AGE, pEPs, or pCEPs.

## Sequence-level decision

| Sequence | Annotation status in the archived run | Clean-background reference | Permitted use in this release |
| --- | --- | --- | --- |
| Board | Foreground segmentation ground truth available | Not established | Structural/operational diagnostics; annotation audit only |
| CAVIAR1 | Foreground segmentation ground truth available | Not established | Structural/operational diagnostics; annotation audit only |
| Hall_Monitor | No ground-truth field in the converted MAT used for the archived run | Not established | Structural/operational diagnostics only |
| People_Foliage | No ground-truth field in the converted MAT used for the archived run | Not established | Structural/operational diagnostics only |

Consequently, `results/sbi/sbi_final_summary.csv` reports tensor dimensions, fitted rank, selected fraction, rESS, iteration count, final relative change, tolerance status, and explicit semantic-status fields. Its PSNR, SSIM, and AGE entries are `NaN`; MS-SSIM, pEPs, and pCEPs are not evaluated.

## v1.0.1 code behavior

The SBI runner no longer interprets the existence of `GT` or `BgtUint8` as proof of a valid clean reconstruction reference. Its default behavior is:

```text
HasCleanReference = false
ValidReconstructionReference = false
PSNR = NaN
SSIM = NaN
AGE = NaN
ReconstructionMetricsStatus = not_applicable_no_valid_clean_background_reference
```

A separately supplied clean background is evaluated only when its MAT file explicitly declares both:

```matlab
GroundTruthType = 'clean_background_reference';
ValidReconstructionReference = true;
```

The runner additionally verifies numeric type, finite values, two-dimensional shape, and spatial agreement with the video. These safeguards prevent a segmentation annotation from being evaluated accidentally as a recovered background.

## Archived compact MAT file

Release v1.0.0 used the field name `hasGT` in `results/sbi/SBI_review_small.mat`. That name indicated the presence of a ground-truth-like field but did not distinguish segmentation annotation from clean reconstruction reference.

Release v1.0.1 removes that ambiguous field and records:

- `hasAnnotation`;
- `hasSegmentationGT`;
- `spatialAlignmentVerified`;
- `frameIndexAlignmentVerified`;
- `groundTruthType`;
- `hasCleanReference`;
- `validReconstructionReference`;
- `reconstructionMetricsStatus`.

The ranks, selected fractions, rESS values, update counts, final relative changes, tolerance flags, and compact attribution summaries are unchanged.

## Manuscript interpretation

The appropriate claim supported by the SBI runs is that the hierarchical selection and weighting mechanism remains operational across heterogeneous real-video tensors. The SBI results do not establish superiority in clean-background reconstruction accuracy.
