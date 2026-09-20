# Urban preprocessing protocol

Raw HYDICE Urban data are not redistributed by this repository.

Expected source file: `Urban_F210.mat`.

Final preprocessing:
- stored matrix layout: 210 × 94249 (bands × pixels);
- metadata: 307 × 307 × 210;
- reconstruct to 307 × 307 × 210;
- retain all 210 bands;
- global min-max normalization before cropping;
- crop rows and columns 26:281 -> 256 × 256 × 210;
- rank: 25;
- HCL recovery floor: epsRec = 1e-3;
- no clean reference is assumed, so Urban is used for structural diagnostics rather than clean-cube NRE/PSNR claims.
