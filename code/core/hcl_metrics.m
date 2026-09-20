function M = hcl_metrics(Xhat,D,hclOut)
% HCL_METRICS  Common synthetic metrics.

M = struct();
M.NRE = norm(Xhat(:)-D.Xstar(:)) / norm(D.Xstar(:));

if nargin < 3 || isempty(hclOut)
    return;
end

trueMask = D.truth.unionMask;
detMask = hclOut.final.po >= 0.50;
recMask = hclOut.mask;

[M.DetPrecision,M.DetRecall,M.DetF1] = binary_prf(detMask,trueMask);
[M.MaskPrecision,M.MaskRecall,M.MaskF1] = binary_prf(recMask,trueMask);
M.AP = average_precision(hclOut.final.po,trueMask);
M.rESS = hclOut.recoveryESS;
M.SelectedFraction = hclOut.selectedFraction;
M.Iterations = hclOut.iterations;
M.FinalRelChange = hclOut.finalRelChange;
M.ToleranceMet = hclOut.toleranceMet;

if ~isempty(D.truth.level)
    predLevel = zeros(size(detMask),'uint8');
    A = cat(4,hclOut.final.piEntry,hclOut.final.piBlock,hclOut.final.piSlice);
    [~,lab] = max(A,[],4);  % ties resolved entry->block->slice
    predLevel(detMask)=uint8(lab(detMask));
    M.LevelMacroF1 = macro_f1_present(predLevel,D.truth.level);
else
    M.LevelMacroF1 = NaN;
end
end

% ========================================================================
function [P,R,F1] = binary_prf(pred,truth)
pred=logical(pred); truth=logical(truth);
tp=nnz(pred & truth); fp=nnz(pred & ~truth); fn=nnz(~pred & truth);
P = tp/max(tp+fp,1);
R = tp/max(tp+fn,1);
F1 = 2*tp/max(2*tp+fp+fn,1);
end

function m = macro_f1_present(pred,truth)
classes = unique(double(truth(:)));
classes(classes==0)=[];
vals = nan(numel(classes),1);
for ii=1:numel(classes)
    c=classes(ii);
    [~,~,vals(ii)] = binary_prf(pred==c,truth==c);
end
m = mean(vals,'omitnan');
end

function ap = average_precision(score,truth)
truth=logical(truth(:)); score=score(:);
npos=nnz(truth);
if npos==0, ap=NaN; return; end
[~,ord]=sort(score,'descend');
y=truth(ord);
tp=cumsum(y);
prec=tp./(1:numel(y))';
ap=sum(prec(y))/npos;
end
