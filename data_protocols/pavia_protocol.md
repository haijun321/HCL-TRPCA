# Pavia University controlled protocol

Raw Pavia University data are not redistributed by this repository.

Expected source file: `PaviaU.mat`.

Final protocol:
- source size: 610 × 340 × 103;
- crop rows 178:433 and columns 43:298 -> 256 × 256 × 103;
- normalize using the supplied [0,8000] intensity range;
- frozen rank: r = 25;
- scenarios:
  - Entry5
  - Block10
  - Slice10
  - Mixed E5+B10+S10
- generation order: slice -> block -> entry;
- supports are mutually exclusive;
- block sizes: [8,12,16,24];
- corruption magnitude: Uniform[6,10] × std(reference), independent random signs;
- seeds: 20280001–20280010;
- HCL recovery floor: epsRec = 1e-3.

Final runner:
`code/pavia/run_PaviaU_CONTROLLED_FINAL_V3_2.m`.
