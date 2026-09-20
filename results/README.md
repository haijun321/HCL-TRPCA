# Archived experimental results

These files contain saved experimental outputs. They are separate from the working directories populated by the reproduction scripts.

| Location | Evidence | Included entry point |
| --- | --- | --- |
| `synthetic/confirmation/` | Confirmation measurements, summaries, and paired statistical outputs | `reproduce_synthetic('CONFIRMATION')` generates raw measurements and basic NRE summaries; paired-statistics postprocessing is not included |
| `synthetic/gate/` | Attribution-rule comparison on shared HCL states | `reproduce_synthetic('GATE')`, after generating confirmation caches |
| `synthetic/rank/` | Rank-misspecification measurements | `reproduce_synthetic('RANK')` |
| `synthetic/robustness/` | Weak-amplitude and overlapping-support diagnostics | `reproduce_synthetic('WEAK_OVERLAP')` |
| `synthetic/stability/` | Selection-time evidence stability | `run_evidence_stability_final()` |
| `synthetic/runtime/` | Repeated timing measurements and seed-level summaries | `run_runtime_synthetic_final('FINAL')` |
| `pavia/` | Controlled corruption protocol and four-method evaluation | `reproduce_pavia('PROTOCOL')` and `reproduce_pavia('FULL')` |
| `urban/` | Compact Urban diagnostic structure | `reproduce_urban()` produces the full result; compact extraction is not included |
| `sbi/` | Sequence diagnostics and compact attribution summaries | `reproduce_sbi()` produces the full results; compact extraction is not included |

## Field conventions

- `LevelMacroF1` and related historical `Level` fields refer to contamination **granularity**, not severity. Their names are retained for compatibility with the original tables and code.
- Pavia `MPSNR` and `MSSIM` are mean bandwise metrics. `SAMdeg` is in degrees.
- Pavia references are operational references constructed from the supplied data and preprocessing protocol; they are not independently measured noise-free ground truth.
- `ToleranceMet` records a relative-step stopping criterion. It does not certify reference-based recovery accuracy or global optimality.
- SBI segmentation annotations are not clean-background references. The archived SBI CSV therefore records PSNR, SSIM, and AGE as `NaN`; MS-SSIM, pEPs, and pCEPs are not evaluated. Its non-metric values are structural, semantic-status, and stopping diagnostics.
- The summary filename `Fig_S1_confirmation_paired_nre_differences_summary.csv` belongs to the archived statistical analysis. Its `Fig_S1` prefix is not a mapping to the later trajectory figure in the manuscript supplement.

The redundant Pavia MAT snapshot containing local filesystem metadata is not distributed. Its four result/protocol tables are retained as `Pavia_V3_raw.csv`, `Pavia_V3_summary.csv`, `Pavia_V3_protocol_raw.csv`, and `Pavia_V3_protocol_summary.csv`. Scientific protocol settings remain in the runner and protocol document.

The compact Urban MAT file retains its original form. In v1.0.1, the compact SBI MAT file preserves the archived structural values but replaces the ambiguous legacy `hasGT` field with explicit annotation and reconstruction-reference fields. For Board and CAVIAR1, `hasSegmentationGT` is true and `validReconstructionReference` is false; for all four sequences, `hasCleanReference` and `validReconstructionReference` are false. See [`SBI_REFERENCE_SEMANTICS.md`](../SBI_REFERENCE_SEMANTICS.md).

The absence of raw datasets, full recovery tensors, intermediate confirmation caches, or external solvers should not be interpreted as zero-valued or missing experimental results.
