function m = macro_f1_present_strict(pred,truth)
%MACRO_F1_PRESENT_STRICT
% Exact copy of the official hcl_metrics.m Macro-F1 definition.

classes = unique(double(truth(:)));
classes(classes==0) = [];

vals = nan(numel(classes),1);

for ii = 1:numel(classes)
    c = classes(ii);
    [~,~,vals(ii)] = binary_prf_strict(pred==c,truth==c);
end

m = mean(vals,'omitnan');
end

function [P,R,F1] = binary_prf_strict(pred,truth)
pred  = logical(pred);
truth = logical(truth);

tp = nnz(pred & truth);
fp = nnz(pred & ~truth);
fn = nnz(~pred & truth);

P = tp/max(tp+fp,1);
R = tp/max(tp+fn,1);
F1 = 2*tp/max(2*tp+fp+fn,1);
end
