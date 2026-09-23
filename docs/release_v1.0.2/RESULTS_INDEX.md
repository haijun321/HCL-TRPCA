# HCL-TRPCA v1.0.2 result index

This update package contains the additional archived outputs used by the submission-stage manuscript and Supplementary Material.

| Manuscript item | Released material |
|---|---|
| Table 3: HCL vs p-TRPCA paired NRE | `results/statistics/confirmation_paired_differences.csv`, `confirmation_paired_stats.csv` |
| Table 6: Calibrated Flat attribution control | `results/attribution_control/synthetic_gate_*` and `configs/calibrated_flat_parameters.csv` |
| Figure 3: recovery-floor / mask-source control | `figures/HCL_recovery_mechanism.pdf`, `results/recovery_controls/*` |
| Figure 6: Pavia rank-budget diagnosis | `results/pavia_rank_budget/*`, `figures/Pavia_rank_budget_diagnostic.*` |
| Supplement Table S5 | `results/pavia_rank_budget/*` |
| Supplement Figure S1 / Table S6 | `results/pavia_fixedmask_trajectory/*`, `figures/Pavia_fixedmask_trajectory_publication_V2.pdf` |

## Important scope notes

- Pavia rank-budget files in this package are the final V2.1 DIAG outputs only; PILOT, REPRO, SMOKE, duplicate and partial files were removed.
- Fixed-mask trajectory files are the FINAL outputs only; PILOT duplicates were removed.
- The recovery-mechanism control is released as the final manuscript figure plus archived table summaries. The exact legacy development runner was not present in the staging archive and is not represented as independently rerunnable source.
- `scripts/reproduce_confirmation_paired_statistics_from_archived_pairs.m` is a new release reproduction/audit utility. It is **not** the exact legacy `POSTPROCESS_CONFIRMATION_PAIRED_STATS.m` used during manuscript preparation; it validates deterministic quantities against the archived frozen CSV outputs.
- Raw Pavia/Urban/SBI datasets and third-party solver source are not bundled.

### Paired-statistics audit note

The v1.0.2 paired-statistics utility is a release-side audit implementation,
not the exact legacy postprocessor used during manuscript preparation.
Deterministic quantities (paired differences, means, medians, win rates,
Wilcoxon p-values, and Holm-adjusted p-values) reproduce the frozen archive.
Bootstrap confidence-interval endpoints may differ from the archived values
at approximately the 1e-6 level. The frozen
`results/statistics/confirmation_paired_stats.csv` remains the numerical
source used in the manuscript.