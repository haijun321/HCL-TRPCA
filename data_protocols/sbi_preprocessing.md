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

An optional static background reference may be supplied as `BgtUint8` or `GT`. RGB tensors require conversion before this runner; the runner expects the third mode to index frames, not color channels. The release does not establish a new conversion, resize, or frame-selection procedure.

The archived sequence dimensions and fitted ranks are:

| Sequence | Height | Width | Frames | Rank |
| --- | ---: | ---: | ---: | ---: |
| Board | 164 | 200 | 228 | 4 |
| CAVIAR1 | 256 | 384 | 610 | 3 |
| Hall_Monitor | 240 | 352 | 296 | 7 |
| People_Foliage | 240 | 320 | 341 | 15 |

Ranks use the runner's 98% Fourier-energy rule, capped at 15. The frozen HCL settings apply to all sequences.

## Interpretation of outputs

The archived `results/sbi/sbi_final_summary.csv` contains structural and stopping diagnostics. It does not report reference-based recovery metrics.

The existing runner also calculates optional metrics if a reference field is present: it takes the temporal median of the recovered tensor, clips it to [0,1], and compares it with the supplied background. Its AGE is mean absolute error on [0,1], not an automatically converted 8-bit benchmark score. These computations do not verify reference semantics, pixel registration, grayscale conversion, or official benchmark conventions. They must not be treated as validated benchmark measurements solely because they are returned by the runner.

Entry point: `reproduce_sbi()`.
