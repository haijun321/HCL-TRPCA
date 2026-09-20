function R = reproduce_urban()
%REPRODUCE_URBAN Reproduce the final Urban run.
setup_hcl_trpca();
R = run_realdata_standardized('URBAN');
end
