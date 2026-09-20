function out = my_TNN_adapter(Y,Xstar,truth,C)
% MY_TNN_ADAPTER
% Adapter for the verified V33 local reference implementation of
% convex TNN-TRPCA via ADMM.
%
% IMPORTANT:
% This is NOT the original authors' official implementation.
%
% Interface expected by HCL_DSP_FINAL:
%   out = adapter(Y,Xstar,truth,C)
% Required:
%   out.Xhat

% ------------------------------------------------------------
% Use the same V33 stopping configuration.
% ------------------------------------------------------------
opt = struct();
opt.tol     = 1e-6;
opt.maxIter = 500;

% ------------------------------------------------------------
% Run the verified V33 TNN-TRPCA implementation.
% Xstar is supplied only for diagnostic NRE tracing and is NOT
% used in the optimization or stopping rule.
% ------------------------------------------------------------
out = hcldref.tnn_trpca_admm(Y,opt,Xstar);

% ------------------------------------------------------------
% Continuous anomaly score for optional AP evaluation.
% ------------------------------------------------------------
if isfield(out,'S')
    out.score = abs(out.S);
else
    out.score = [];
end

% Record provenance.
out.adapterName = 'V33 local reference TNN-TRPCA-ADMM';
out.isOfficialAuthorsCode = false;

end