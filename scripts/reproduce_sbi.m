function T = reproduce_sbi()
%REPRODUCE_SBI Reproduce the final SBI runs.
setup_hcl_trpca();
T = run_sbi_standardized();
end
