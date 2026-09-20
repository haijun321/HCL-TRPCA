# Synthetic protocol

- Tensor size: 64 × 64 × 30.
- True tubal rank: 5.
- Additive Gaussian noise: 35 dB SNR.
- Gross corruption magnitudes: Uniform[6,10] with independent random signs.
- Moderate scenarios:
  - EB: 5% entry + 10% block
  - ES: 5% entry + 10% slice
  - BS: 10% block + 10% slice
  - EBS: 5% entry + 10% block + 10% slice
- Supports are generated as disjoint labels in the controlled benchmark.
- Development seeds: 20261001–20261005.
- Independent confirmation seeds: 20270001–20270030.
- Weak EBS uses amplitude range [3,5].
- Rank-misspecification grid: [3,4,5,6,8].

See `code/synthetic/` for the exact final generator and experiment entry points.
