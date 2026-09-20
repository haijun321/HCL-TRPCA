function R = reproduce_synthetic(mode)
%REPRODUCE_SYNTHETIC Reproduce the frozen synthetic experiments.
if nargin<1, mode='SMOKE'; end
setup_hcl_trpca();
R = run_HCL_DSP_FINAL(mode);
end
