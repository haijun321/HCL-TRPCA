Development-stage recovery-mechanism control.

Experimental unit:
five development seeds.

Mask sources:
ENTRY, count-matched ENTRY-M, HCL.

Recovery floors:
0.01, 0.001, 0.

Common recovery settings:
rank = 5
relaxation = 0.85
maximum recovery updates = 320
relative-step tolerance = 1e-5

The control changes the recovery floor without recomputing scores or
reselecting positions. ENTRY-M matches the number of positions selected
by HCL without using true labels.

The archived files reproduce the manuscript-level summaries and Figure 3.
The exact legacy development runner is not included in the current frozen
package and is therefore not represented as independently rerunnable source.

Files
-----
- `recovery_floor_summary.csv`: machine-readable version of the archived manuscript-level floor summary.
- `mask_source_zero_floor_summary.csv`: machine-readable version of the archived manuscript-level mask-source summary.
- `hcl_recovery_floor.tex` and `weight_source_zero_floor.tex`: archived LaTeX table fragments used during manuscript preparation.
- `../../figures/HCL_recovery_mechanism.pdf`: final manuscript Figure 3 artifact.

The two CSV summaries were transcribed directly from the archived LaTeX tables in this release; no experiment was rerun or retuned to create them.
