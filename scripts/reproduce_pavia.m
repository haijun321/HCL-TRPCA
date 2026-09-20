function R = reproduce_pavia(mode)
%REPRODUCE_PAVIA Reproduce the final controlled PaviaU experiment.
if nargin<1, mode='INTERFACE'; end
setup_hcl_trpca();
R = run_PaviaU_CONTROLLED_FINAL_V3_2(mode);
end
