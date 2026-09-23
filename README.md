# HCL-TRPCA

MATLAB code and compact experimental results for **HCL-TRPCA: Hierarchical Contamination-Level Attribution for Robust Tensor PCA under Mixed Structured Corruption**.

**Authors:** Haijun Jiang, Li Han, Yue Wang, Jiang Fen, and Ganggang Xu.

HCL-TRPCA separates abnormality detection, entry/block/slice granularity attribution, and numerical recovery weighting. A warm start supplies evidence for selecting a frozen recovery mask. Recovery then uses fixed weights and a prescribed tubal-rank projection. Attribution computed from the final residual is diagnostic and does not feed back into recovery.

Hierarchy is evaluated primarily through granularity attribution. Recovery also depends on suppression strength, rank, and the frozen weighted criterion; the method does not provide uniformly better reconstruction across datasets.

This archive is **release v1.0.2**, synchronized with the submission-stage
HCL-TRPCA manuscript and Supplementary Material. It extends the v1.0.1
public snapshot with attribution-control materials, recovery-mechanism
artifacts, Pavia rank-budget and fixed-mask diagnostics, paired statistical
audit materials, and corresponding publication figures. The HCL-TRPCA core
method, frozen primary configuration, development/confirmation seed split,
and primary reported numerical results are unchanged. The SBI ground-truth
semantics clarification introduced in v1.0.1 is retained. See
[CHANGELOG.md](CHANGELOG.md),
[SBI_REFERENCE_SEMANTICS.md](SBI_REFERENCE_SEMANTICS.md), and
[`docs/release_v1.0.2/RESULTS_INDEX.md`](docs/release_v1.0.2/RESULTS_INDEX.md).

## v1.0.2 submission reproducibility update

Release v1.0.2 synchronizes the public repository with the
submission-stage HCL-TRPCA manuscript and Supplementary Material.

Relative to v1.0.1, this release adds:

- development-calibrated nonhierarchical attribution-control materials
  and frozen calibration coefficients;
- recovery-floor and mask-source mechanism-control artifacts;
- Pavia University rank-budget V2.1 DIAG results;
- Pavia fixed-mask recovery trajectories and the corresponding
  publication figure;
- archived paired confirmation statistics and a release-side
  reproduction/audit utility;
- result-to-manuscript mapping and portability documentation.

The HCL-TRPCA core method, frozen primary configuration, development
seeds (20261001--20261005), confirmation seeds
(20270001--20270030), and primary numerical results are unchanged.

See `docs/release_v1.0.2/RESULTS_INDEX.md` for the mapping between
released numerical files and manuscript tables/figures.

### External data location

Raw third-party datasets are not redistributed. Set the environment
variable `HCL_TRPCA_DATA_ROOT` to the directory containing the required
datasets, or place the expected files under `<repository>/data`.
## Contents

| Directory or file | Contents |
| --- | --- |
| `code/core/` | HCL solver, tubal-rank projection, and metrics |
| `code/synthetic/` | Synthetic generator, confirmation, attribution, robustness, rank, stability, and runtime experiments |
| `code/pavia/` | Controlled Pavia University experiment |
| `code/urban/` | Standardized Urban and uncorrupted Pavia diagnostic runner |
| `code/sbi/` | SBI sequence diagnostic runner |
| `code/baselines/` | Interfaces to external TNN-TRPCA and p-TRPCA implementations |
| `configs/` | Frozen scientific settings and a local data-path template |
| `scripts/` | Setup checks and reproduction entry points |
| `seeds/` | Development, confirmation, robustness, rank, and Pavia seeds |
| `data_protocols/` | Input formats and preprocessing protocols |
| `results/` | Archived numerical results; see the [result index](results/README.md) |
| `figures/` | Four archived experiment figures |
| `environment/` | Recorded MATLAB environment |
| `external/` | External solver interfaces and availability limitations |
| `SBI_REFERENCE_SEMANTICS.md` | SBI ground-truth audit and metric-validity decision |
| `CHANGELOG.md` | Release history |
| `SHA256SUMS` | SHA-256 checksums for the release files |

## Requirements

The recorded experimental environment is MATLAB R2023b, version 23.2.0.2365128, on 64-bit Windows. The full installed-product inventory is in [environment/matlab_version_full.txt](environment/matlab_version_full.txt); it is not a list of required toolboxes.

The core solver uses MATLAB array operations, FFT, SVD, and convolution. Pavia SSIM evaluation uses `ssim` from Image Processing Toolbox; the runner records `NaN` when that function is unavailable. External baseline implementations have their own requirements. Compatibility with other MATLAB releases or GNU Octave is not established by this archive.

## Quick start

Extract the archive, open MATLAB in the directory containing `setup_hcl_trpca.m`, and run:

```matlab
setup_hcl_trpca;
report = preflight_public_release();
R = reproduce_synthetic('SMOKE');
```

The smoke test uses a generated 64-by-64-by-30 tensor and requires neither external datasets nor baseline solvers. It checks the two frozen recovery weights and saves its output under `generated_results/`.

`preflight_public_release` checks configuration loading, entry-point availability, and the compact v1.0.1 SBI semantic artifacts. It is not a numerical validation of every experiment. Missing external solvers are reported separately; they do not prevent the HCL smoke test or the SBI semantic audit from running. The semantic audit can also be invoked directly without raw data:

```matlab
report = audit_sbi_release_semantics();
```

## Synthetic experiments

Development seeds are `20261001:20261005`; independent confirmation seeds are `20270001:20270030`. Scientific settings are defined in [configs/HCL_DSP_FINAL_config.m](configs/HCL_DSP_FINAL_config.m).

| MATLAB command | Scope and dependencies |
| --- | --- |
| `reproduce_synthetic('SMOKE')` | One HCL run; no external solvers |
| `reproduce_synthetic('CONFIRMATION')` | Four scenarios and 30 seeds; requires the configured TNN and p-TRPCA solvers |
| `reproduce_synthetic('GATE')` | Attribution rules sharing saved HCL states; requires confirmation caches |
| `reproduce_synthetic('WEAK_OVERLAP')` | Weak-amplitude and overlapping-support HCL diagnostics |
| `reproduce_synthetic('RANK')` | HCL rank-misspecification experiment |
| `run_evidence_stability_final()` | Selection-time evidence perturbation diagnostic |
| `run_runtime_synthetic_final('FINAL')` | Frozen baseline runtime comparison; requires external solvers |

Run `setup_hcl_trpca` before invoking the lower-level entry points. Full confirmation writes intermediate HCL states for the attribution experiment and can use substantially more disk space than this archive. Archived result tables remain under `results/`; new synthetic outputs are written under `generated_results/`.

## External baseline dependencies

See [external/README.md](external/README.md) before running a complete comparison. The adapters are included; the third-party solver sources are not.

- **p-TRPCA:** obtain the authors' `ptrpca` implementation from the [authors' repository](https://github.com/qguo2010/p-TRPCA) and add its implementation directory to the MATLAB path.
- **TNN-TRPCA:** the reported experiments used the local callable `hcldref.tnn_trpca_admm`, which is not distributed here. Complete baseline reruns require access to that exact implementation. Substituting another TNN solver constitutes a separate comparison and need not reproduce the archived numbers.

The configuration retains both baseline adapters. A complete comparison will stop if a required solver is unavailable; the release does not silently replace it or omit its results.

## Local datasets

Raw datasets are not included. From the repository root, create a local configuration:

```matlab
copyfile(fullfile('configs','local_paths.example.txt'), ...
         fullfile('configs','local_paths.m'));
edit(fullfile('configs','local_paths.m'));
```

Set `P.dataRoot` and, where needed, the individual dataset paths. Alternatively, set the environment variable `HCL_TRPCA_DATA_ROOT`. Keep `configs/local_paths.m` untracked.

| Dataset | Expected input | Protocol |
| --- | --- | --- |
| Pavia University | `PaviaU.mat`, containing the numeric hyperspectral cube | [Pavia](data_protocols/pavia_protocol.md) |
| HYDICE Urban | `Urban_F210.mat`, with `Y`, `nRow`, `nCol`, and `nBand` | [Urban](data_protocols/urban_preprocessing.md) |
| SBI | Four converted sequence MAT files under `MAT/` or the configured `P.sbiRoot` | [SBI](data_protocols/sbi_preprocessing.md) |

Dataset acquisition and SBI source-to-MAT conversion are external to this release. Retain the original providers' dataset attribution and terms.

## Pavia University

After configuring `PaviaU.mat`:

```matlab
reproduce_pavia('INTERFACE');
reproduce_pavia('PROTOCOL');
R = reproduce_pavia('AUDIT');
```

`INTERFACE` checks the HCL interface; the current runner also loads the configured Pavia cube before this check. `PROTOCOL` verifies the generated supports. `AUDIT` evaluates HCL on four scenarios and ten seeds. `FULL` additionally evaluates tSVD-LRA, TNN-TRPCA, and p-TRPCA and requires both external solvers:

```matlab
R = reproduce_pavia('FULL');
```

The controlled experiment uses the prespecified rank `r = 25`. If the separate `hcl_estimate_tubal_rank` helper is available on the path, the runner checks the stored capped rank rule; otherwise it uses the fixed rank. That optional helper is not distributed here. New results are written to timestamped directories under `code/pavia/results_HCL_DSP_FINAL/`.

The `PAVIA` mode in `run_realdata_standardized` is an uncorrupted-data diagnostic and is distinct from the controlled corruption benchmark above.

## Urban and SBI diagnostics

```matlab
R = reproduce_urban();
T = reproduce_sbi();
```

Urban retains all 210 bands and uses rank 25. Without a clean reference, its outputs support structural and weighting diagnostics, not clean-cube reconstruction accuracy.

SBI uses a 98% Fourier-energy rank rule capped at 15. The archived SBI CSV reports ranks, selection fractions, relative effective sample sizes, iteration counts, and stopping diagnostics. The distributed SBI `groundtruth/gt*.png` files audited for Board and CAVIAR1 are foreground segmentation masks rather than clean-background reference images. They are therefore not used to calculate PSNR, SSIM, MS-SSIM, AGE, pEPs, or pCEPs.

The v1.0.1 runner defaults all SBI reconstruction-reference metrics to `NaN`. It enables PSNR, SSIM, and AGE only for a separately supplied clean background that is explicitly tagged with both `GroundTruthType = 'clean_background_reference'` and `ValidReconstructionReference = true`. Field presence alone is insufficient. The manuscript's SBI evidence is consequently interpreted as real-video structural and operational validation, not as a reference-based background-recovery benchmark. See the [SBI protocol](data_protocols/sbi_preprocessing.md) and [reference-semantics note](SBI_REFERENCE_SEMANTICS.md).

## Results and coverage

The [result index](results/README.md) maps archived tables to their experiments and explains historical field names. Numeric CSV values, seeds, and the HCL core are retained from the supplied experimental archive.

This release contains the experiments listed above. It does **not** contain the later Pavia rank-by-budget and fixed-mask trajectory packages, recovery-floor/mask-source controls, or scripts for every manuscript figure and statistical postprocessing step. The included paired-statistics CSV files are archived outputs. They are not a substitute for the absent postprocessing scripts. The four PDFs are archived experiment figures, not a complete set of final manuscript figures.

Full baseline reproduction additionally depends on the unavailable local TNN solver. These coverage limits should be taken into account when citing this software snapshot.

## Citation and license

Citation metadata is provided in [CITATION.cff](CITATION.cff). No publication DOI is assigned in this release metadata.

Original project code is provided under the [BSD 3-Clause License](LICENSE). External software and datasets retain their own terms.
