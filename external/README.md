# External baseline dependencies

The HCL-TRPCA repository contains the proposed method, experimental
protocols, project-specific baseline adapters, and compact result files
used in the manuscript.

Some baseline solver implementations are not redistributed because
their redistribution rights could not be established with sufficient
certainty from the tested local snapshots.

## TNN-TRPCA

Project adapter:

`code/baselines/my_TNN_adapter.m`

The manuscript experiments used a verified local ADMM reference
implementation of convex TNN-TRPCA.

This implementation is not the official implementation distributed by
the original TNN-TRPCA authors.

The local solver source is not redistributed in this public release
because no explicit license, README, Git metadata, or other
redistribution metadata was found in the tested local snapshot.

The adapter is retained to document the exact interface and experimental
configuration used in the manuscript.

The local stopping configuration used in the reported experiments was:

- tolerance: `1e-6`
- maximum iterations: `500`

In the audited local implementation, the clean reference was used only
for diagnostic NRE tracing. It did not affect optimization updates,
parameter selection, stopping, or iterate selection.

Compact numerical results obtained with this baseline are included under
`results/`.

Alternative public implementations of convex TNN-TRPCA may be used for
independent comparison, but they are not guaranteed to reproduce the
exact numerical values of the local reference implementation used in
the manuscript.

## p-TRPCA

Reference:

T. Yan and Q. Guo,
"Tensor robust principal component analysis via dual lp quasi-norm
sparse constraints,"
Digital Signal Processing, vol. 150, Art. no. 104520, 2024.

Official code repository:

[qguo2010/p-TRPCA](https://github.com/qguo2010/p-TRPCA).

Required MATLAB callable:

`ptrpca`

The manuscript experiments used the authors' distributed implementation
through:

`code/baselines/my_pTRPCA_adapter.m`

The third-party p-TRPCA implementation is not redistributed in this
repository. Users should obtain it from the authors' repository and add
the corresponding implementation directory to the MATLAB path.

The frozen p-TRPCA parameters used in the manuscript are:

- `w = [1, 1.1, 1.5]`
- `p1 = 1.0`
- `p2 = 0.8`
- `mu = 1e-4`
- `tol = 1e-8`
- `rho = 1.1`
- `maxIter = 500`

## Reproduction scope

The HCL-TRPCA method and standalone synthetic smoke test do not require
the external p-TRPCA implementation.

Reproducing the complete baseline comparison requires the corresponding
external baseline dependencies described above.

The exact compact result tables used in the manuscript are included
under `results/`.