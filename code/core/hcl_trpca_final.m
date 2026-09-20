function out = hcl_trpca_final(Y, r, P)
% HCL_TRPCA_FINAL
% Implements the manuscript HCL-TRPCA frozen-mask algorithm.
%
% Outputs:
%   Xhat, mask, W0, recoveryESS, iterations, finalRelChange
%   selection: ae,ab,as,po,qe,qb,qs,piClean,piEntry,piBlock,piSlice
%   final    : same fields recomputed from final residual
%
% NOTE: selection uses epsSel; recovery uses epsRec. They are distinct.

validateattributes(Y,{'double','single'},{'real','finite','nonempty'});
[I,J,K] = size(Y);
M = numel(Y);

%% 1) Initialize
X = hcl_project_tubal_rank(Y,r);

%% 2) Robust warm start
Tw = min(P.warmIterations,P.maxIterations);
for t = 1:Tw
    E = compute_evidence(Y-X,P);
    Wwarm = max(P.epsWarm, ...
        tukey_weight(E.ze,P.tukeyEntry) .* ...
        tukey_weight(E.zb,P.tukeyBlock) .* ...
        tukey_weight(E.zs,P.tukeySlice));
    Xnew = relaxed_update(Y,X,Wwarm,r,P.gamma);
    X = Xnew;
end

warmX = X;

%% 3) Selection-time attribution
Sel = compute_evidence(Y-X,P);
Sel = add_attribution(Sel,P);

%% 4) Select and freeze O0
[mask,n0,nstar,etaSel] = select_mask(Sel.po,P);

%% 5) Frozen recovery weights
W0 = ones(size(Y),'like',Y);
W0(mask) = P.epsRec;
recoveryESS = normalized_ess(W0);

%% 6) Fixed-weight recovery
finalRel = Inf;
iter = Tw;

while iter < P.maxIterations
    Xnew = relaxed_update(Y,X,W0,r,P.gamma);
    finalRel = norm(Xnew(:)-X(:)) / max(norm(X(:)),eps(class(X)));
    X = Xnew;
    iter = iter + 1;
    if finalRel < P.relTol
        break;
    end
end

%% 7) Final diagnostic attribution
Fin = compute_evidence(Y-X,P);
Fin = add_attribution(Fin,P);

out = struct();
out.Xhat = X;
out.warmX = warmX;
out.mask = mask;
out.W0 = W0;
out.recoveryESS = recoveryESS;
out.selectionESS = etaSel;
out.n0 = n0;
out.nstar = nstar;
out.selectedFraction = nnz(mask)/M;
out.iterations = iter;
out.finalRelChange = finalRel;
out.toleranceMet = (finalRel < P.relTol);
out.selection = Sel;
out.final = Fin;
out.rank = r;
out.epsRec = P.epsRec;
end

% ========================================================================
function Xnew = relaxed_update(Y,X,W,r,gamma)
V = W .* Y + (1-W) .* X;
Z = hcl_project_tubal_rank(V,r);
Xnew = gamma .* Z + (1-gamma) .* X;
end

% ========================================================================
function E = compute_evidence(R,P)
[I,J,K] = size(R);

% Warm-start RMS statistic B, full b^2 denominator with zero padding.
b = P.b;
kerB = ones(b,b);
B = zeros(size(R),'like',R);
for k = 1:K
    B(:,:,k) = sqrt(max(conv2(R(:,:,k).^2,kerB,'same')/(b^2),0));
end

% Slice RMS
S = zeros(K,1);
for k = 1:K
    tmp = R(:,:,k);
    S(k) = sqrt(mean(tmp(:).^2));
end

% Robust normalized residual quantities
mR = median(R(:));
sR = robust_scale(R(:),P.scaleEps);
ze = abs(R-mR) ./ max(sR,P.scaleEps);

mB = median(B(:));
sB = robust_scale(B(:),P.scaleEps);
zb = max(B-mB,0) ./ max(sB,P.scaleEps);

mS = median(S);
sS = robust_scale(S,P.scaleEps);
zsVec = max(S-mS,0) ./ max(sS,P.scaleEps);
zs = zeros(size(R),'like',R);
for k = 1:K
    zs(:,:,k) = zsVec(k);
end

% Entry evidence
ae = sigmoid_clip((ze-P.entryCenter)./P.entryTemp);

% Structural evidence
h = min([5,I,J]);
kerH = ones(h,h);
ab = zeros(size(R),'like',R);
Cevb = zeros(size(R),'like',R);

validCount = conv2(ones(I,J),kerH,'same');
for k = 1:K
    Cevb(:,:,k) = conv2(ae(:,:,k),kerH,'same') ./ validCount;
    ab(:,:,k) = sigmoid_clip((Cevb(:,:,k)-P.blockCenter)./P.blockTemp);
end

sliceMean = zeros(K,1);
as = zeros(size(R),'like',R);
for k = 1:K
    tmp = ae(:,:,k);
    sliceMean(k) = mean(tmp(:));
    as(:,:,k) = sigmoid_clip((sliceMean(k)-P.sliceCenter)./P.sliceTemp);
end

E = struct('ze',ze,'zb',zb,'zs',zs, ...
           'ae',ae,'ab',ab,'as',as, ...
           'B',B,'S',S,'Cevb',Cevb,'sliceMean',sliceMean);
end

% ========================================================================
function A = add_attribution(E,P)
ae = E.ae; ab = E.ab; as = E.as;

qs = as;
qb = (1-as).*ab;
qe = (1-as).*(1-ab);

beff = min(1, P.betaBlock .* (ae.^P.nu) .* ab);
seff = min(1, P.betaSlice .* (ae.^P.nu) .* as);

po = 1 - (1-ae).*(1-beff).*(1-seff);
po = min(1,max(0,po));

piC = 1-po;
piE = po.*qe;
piB = po.*qb;
piS = po.*qs;

den = piC+piE+piB+piS;
den(den==0) = 1;

A = E;
A.po = po;
A.qe = qe; A.qb = qb; A.qs = qs;
A.piClean = piC./den;
A.piEntry = piE./den;
A.piBlock = piB./den;
A.piSlice = piS./den;
end

% ========================================================================
function [mask,n0,nstar,eta] = select_mask(po,P)
M = numel(po);
fdhat = mean(po(:) >= P.tauDetect);

countHigh = nnz(po(:) >= P.tauHigh);
n0 = min(round(M*min(fdhat,P.fMax)), countHigh);
n0 = max(0,min(M,n0));

[~,ord] = sort(po(:),'descend');

nstar = n0;
while nstar > 0
    eta = two_level_ess(nstar,M,P.epsSel);
    if eta >= P.etaMin
        break;
    end
    nstar = nstar - 1;
end

if nstar == 0
    eta = 1;
end

mask = false(size(po));
if nstar > 0
    mask(ord(1:nstar)) = true;
end
end

function eta = two_level_ess(n,M,epsFloor)
f = n/M;
eta = (1-f+epsFloor*f)^2 / (1-f+(epsFloor^2)*f);
end

% ========================================================================
function w = tukey_weight(z,c)
u = z./c;
w = (1-u.^2).^2;
w(abs(z) >= c) = 0;
w = max(0,w);
end

function y = sigmoid_clip(x)
x = min(40,max(-40,x));
y = 1 ./ (1+exp(-x));
end

function s = robust_scale(x,epsS)
x = x(isfinite(x));
if isempty(x)
    s = 1;
    return;
end
m = median(x);
d = median(abs(x-m)) / 0.6744897501960817;
v = sqrt(mean((x-m).^2));
if d >= epsS
    s = d;
elseif v >= epsS
    s = v;
else
    s = 1;
end
end

function e = normalized_ess(W)
w = W(:);
e = (sum(w)^2) / (numel(w)*sum(w.^2));
end
