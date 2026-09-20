# SBI input and diagnostic protocol

Raw SBI sequences and the source-to-MAT conversion utility are not distributed.

## Expected files

- `SBI_Board.mat`
- `SBI_CAVIAR1.mat`
- `SBI_Hall_Monitor.mat`
- `SBI_People_Foliage.mat`

Place these files in `P.sbiRoot`, as configured in `configs/local_paths.m`.

## MAT fields

The preferred observation field is `Yuint8`, a grayscale height-by-width-by-frame tensor of type `uint8`. The runner converts it to double precision and divides by 255. Alternatively, `Y` may contain an already normalized tensor; if its maximum exceeds 1, the runner divides it by 255.

RGB tensors require conversion before this runner; the runner expects the third mode to index frames, not color channels. The release does not establish a new conversion, resize, or frame-selection procedure.

Ground-truth-related fields must distinguish annotation semantics from recovery-reference semantics. Recommended metadata for a segmentation mask is:

```matlab
GroundTruthType = 'foreground_segmentation_mask';
HasSegmentationGT = true;
ValidReconstructionReference = false;
```

The annotation may be stored as `SegmentationGTUint8` or `SegmentationGT`. Legacy converted files may instead contain `BgtUint8` or `GT`; these legacy names do not by themselves establish clean-background semantics.

A separately supplied static clean background may be stored as `CleanBackgroundUint8` or `CleanBackground`. For legacy compatibility, `BgtUint8` or `GT` is accepted only when the file also declares:

```matlab
GroundTruthType = 'clean_background_reference';
ValidReconstructionReference = true;
```

The runner then verifies that the reference is numeric, finite, two-dimensional, and spatially matched to the video. Without both explicit declarations, reference-based recovery metrics remain disabled.

The archived sequence dimensions and fitted ranks are:

| Sequence | Height | Width | Frames | Rank |
| --- | ---: | ---: | ---: | ---: |
| Board | 164 | 200 | 228 | 4 |
| CAVIAR1 | 256 | 384 | 610 | 3 |
| Hall_Monitor | 240 | 352 | 296 | 7 |
| People_Foliage | 240 | 320 | 341 | 15 |

Ranks use the runner's 98% Fourier-energy rule, capped at 15. The frozen HCL settings apply to all sequences.

## Interpretation of outputs

The archived `results/sbi/sbi_final_summary.csv` contains structural and stopping diagnostics plus explicit semantic-status fields. Its PSNR, SSIM, and AGE columns are `NaN` because no valid clean-background reference is available. Across the four archived sequences, the selected fraction is approximately 0.107--0.188 and rESS is approximately 0.813--0.893. Board meets the prescribed relative-step tolerance after 52 updates; the other three sequences reach the 80-update budget. The budget-limited runs are operational diagnostics, not evidence that all sequences converged to a common stationary solution.

The SBI ground-truth images audited for Board and CAVIAR1 are frame-wise foreground segmentation masks rather than clean-background reference images. Spatial or frame-index alignment does not change that semantic distinction. Accordingly, these masks are not used to compute PSNR, SSIM, MS-SSIM, AGE, pEPs, or pCEPs. Hall_Monitor and People_Foliage are likewise not assigned reconstruction-reference metrics in this release.

The v1.0.1 runner records `PSNR`, `SSIM`, and `AGE` as `NaN` unless an independently validated clean background is supplied with the explicit metadata above. If enabled for a user-supplied clean reference, the runner compares that reference with the clipped temporal median of the recovered tensor; AGE is then mean absolute error on [0,1]. Such user-supplied measurements are outside the archived SBI evidence and must still be checked against the intended benchmark convention.

The manuscript and this release therefore use SBI as real-video cross-domain structural and weighting validation, not as a reference-based background-recovery superiority benchmark. The full semantic decision record is in [`SBI_REFERENCE_SEMANTICS.md`](../SBI_REFERENCE_SEMANTICS.md).

Entry point: `reproduce_sbi()`.
