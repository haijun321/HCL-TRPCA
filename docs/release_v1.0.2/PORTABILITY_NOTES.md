# Public path sanitization note

For the GitHub-ready copy only, three local-machine path references were replaced by portable resolution rules:

1. `configs/HCL_DSP_FINAL_config.m` uses environment variable `HCL_TRPCA_DATA_ROOT` when set, otherwise `<repository>/data`.
2. `code/pavia/run_PaviaU_RANK_BUDGET_DIAG_V2_1.m` no longer contains author-machine `G:` / `H:` Pavia paths.
3. `code/pavia/run_Pavia_FIXEDMASK_TRAJECTORY_DIAG_V1.m` no longer contains author-machine `G:` / `H:` Pavia paths.

These edits change file discovery only; no scientific parameter, seed, rank, update budget, metric or archived result was changed.
