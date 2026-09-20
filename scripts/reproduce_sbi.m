function T = reproduce_sbi()
%REPRODUCE_SBI Reproduce the SBI structural and operational diagnostics.
% Reference-based metrics remain disabled unless a separately validated
% clean background is explicitly declared; see data_protocols/sbi_preprocessing.md.
setup_hcl_trpca();
T = run_sbi_standardized();
end
