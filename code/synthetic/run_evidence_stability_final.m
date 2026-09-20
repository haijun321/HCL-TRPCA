function R = run_evidence_stability_final()
% RUN_EVIDENCE_STABILITY_FINAL
% Proposition-1 diagnostic for the frozen HCL-TRPCA protocol.
%
% IMPORTANT:
%   1) This diagnostic uses H.selection, because H.mask is selected from the
%      selection-time score p_o^0. It must NOT use H.final.
%   2) It does not modify hcl_trpca_final.m or any frozen HCL parameter.
%   3) It reports both fixed-n top-n mask agreement and full-rule agreement.
%
% Expected sanity check:
%   delta = 0:
%       ActualPoInf             = 0
%       FixedNMaskAgreement     = 1
%       FullRuleMaskAgreement   = 1
%       FullRuleSelectedCount   = H.nstar
%
% The experiment is a diagnostic of Proposition 1, not a new tuning stage.

C = HCL_DSP_FINAL_config();

outDir = fullfile(C.outputRoot,'EVIDENCE_STABILITY_FINAL');
if ~exist(outDir,'dir'), mkdir(outDir); end

deltaGrid = [0 1e-4 1e-3 1e-2 5e-2];
seeds = C.confirmationSeeds(1:min(10,numel(C.confirmationSeeds)));

rows = {};
rid = 0;

for iseed = 1:numel(seeds)
    seed = seeds(iseed);

    D = hcl_synthetic_case(seed,'EBS',C,'DISJOINT');
    H = hcl_trpca_final(D.Y,D.rank,C.hcl);

    % -------------------------------------------------------------
    % Selection-time evidence:
    % Frozen mask H.mask was generated from the SELECTION-TIME score.
    % Therefore the stability diagnostic must perturb H.selection.
    % -------------------------------------------------------------
    E = H.selection;
    basePo   = E.po;
    baseMask = H.mask;
    baseN    = nnz(baseMask);

    % Boundary gap of the original selection-time score for top-n stability.
    boundaryGap = topn_boundary_gap(basePo,baseN);

    for id = 1:numel(deltaGrid)
        delta = deltaGrid(id);

        % Deterministic perturbations for reproducibility.
        pertSeed = mod(double(seed) + 100000*double(id), 2^32-1);
        rng(pertSeed,'twister');

        if delta == 0
            ae2 = E.ae;
            ab2 = E.ab;
            as2 = E.as;
        else
            ae2 = clip01(E.ae + delta*(2*rand(size(E.ae))-1));
            ab2 = clip01(E.ab + delta*(2*rand(size(E.ab))-1));
            as2 = clip01(E.as + delta*(2*rand(size(E.as))-1));
        end

        actualEvidenceDelta = max([ ...
            max(abs(ae2(:)-E.ae(:))), ...
            max(abs(ab2(:)-E.ab(:))), ...
            max(abs(as2(:)-E.as(:))) ]);

        po2 = recompute_po(ae2,ab2,as2,C.hcl);

        actualPoInf = max(abs(po2(:)-basePo(:)));

        % Proposition-1 score bound, using the ACTUAL evidence perturbation.
        B = min(1, ...
            2.4*actualEvidenceDelta + 1.4*sqrt(actualEvidenceDelta));

        boundSatisfied = (actualPoInf <= B + 1e-12);

        % ---------------------------------------------------------
        % A) Fixed-n top-n mask agreement
        % ---------------------------------------------------------
        fixedNMask = topn_mask(po2,baseN);
        fixedNAgreement = mean(fixedNMask(:) == baseMask(:));

        if isnan(boundaryGap)
            topNConditionSatisfied = NaN;
        else
            topNConditionSatisfied = double(boundaryGap > 2*B);
        end

        % ---------------------------------------------------------
        % B) Full frozen-selection-rule agreement
        %    Reapply Eqs. (10)--(13) with unchanged calibration.
        % ---------------------------------------------------------
        [fullMask,n0new,nstarNew,etaNew] = ...
            select_mask_full_rule(po2,C.hcl);

        fullRuleAgreement = mean(fullMask(:) == baseMask(:));

        % ---------------------------------------------------------
        % Mandatory delta=0 sanity checks
        % ---------------------------------------------------------
        if delta == 0
            assert(actualEvidenceDelta < 1e-14, ...
                'delta=0 failed: evidence changed.');
            assert(actualPoInf < 1e-14, ...
                'delta=0 failed: p_o changed.');
            assert(abs(fixedNAgreement-1) < 1e-14, ...
                'delta=0 failed: fixed-n mask mismatch.');
            assert(abs(fullRuleAgreement-1) < 1e-14, ...
                'delta=0 failed: full-rule mask mismatch.');
            assert(nstarNew == H.nstar, ...
                'delta=0 failed: selected count differs from H.nstar.');
        end

        rid = rid + 1;
        rows(rid,:) = { ...
            seed, delta, actualEvidenceDelta, actualPoInf, B, ...
            boundSatisfied, boundaryGap, topNConditionSatisfied, ...
            fixedNAgreement, fullRuleAgreement, ...
            baseN, n0new, nstarNew, etaNew};
    end
end

T = cell2table(rows,'VariableNames', { ...
    'Seed','RequestedDelta','ActualEvidenceDelta','ActualPoInf', ...
    'TheoreticalBound','BoundSatisfied','BoundaryGap', ...
    'TopNConditionSatisfied','FixedNMaskAgreement', ...
    'FullRuleMaskAgreement','BaseSelectedCount','FullRuleN0', ...
    'FullRuleSelectedCount','FullRuleESS'});

writetable(T,fullfile(outDir,'evidence_stability_final_raw.csv'));

S = groupsummary(T,'RequestedDelta',{'mean','std'}, { ...
    'ActualEvidenceDelta','ActualPoInf','TheoreticalBound', ...
    'FixedNMaskAgreement','FullRuleMaskAgreement', ...
    'FullRuleSelectedCount'});

writetable(S,fullfile(outDir,'evidence_stability_final_summary.csv'));

% Overall audit
A = table();
A.NumRuns = height(T);
A.NumBoundViolations = sum(~T.BoundSatisfied);
A.NumDeltaZeroFailures = sum( ...
    T.RequestedDelta==0 & ...
    (abs(T.FixedNMaskAgreement-1)>1e-14 | ...
     abs(T.FullRuleMaskAgreement-1)>1e-14 | ...
     T.ActualPoInf>1e-14));
writetable(A,fullfile(outDir,'evidence_stability_final_audit.csv'));

R = struct('raw',T,'summary',S,'audit',A);
save(fullfile(outDir,'evidence_stability_final_results.mat'),'R');

fprintf('\n=== EVIDENCE STABILITY FINAL ===\n');
fprintf('Runs                : %d\n',A.NumRuns);
fprintf('Bound violations    : %d\n',A.NumBoundViolations);
fprintf('delta=0 failures    : %d\n',A.NumDeltaZeroFailures);

if A.NumBoundViolations==0 && A.NumDeltaZeroFailures==0
    fprintf('STATUS              : PASS\n');
else
    fprintf('STATUS              : CHECK REQUIRED\n');
end

end

% ========================================================================
function po = recompute_po(ae,ab,as,P)
beff = min(1,P.betaBlock.*(ae.^P.nu).*ab);
seff = min(1,P.betaSlice.*(ae.^P.nu).*as);
po = 1-(1-ae).*(1-beff).*(1-seff);
po = clip01(po);
end

function x = clip01(x)
x = min(1,max(0,x));
end

function mask = topn_mask(score,n)
mask = false(size(score));
if n<=0, return; end
[~,ord] = sort(score(:),'descend');
n = min(n,numel(ord));
mask(ord(1:n)) = true;
end

function gap = topn_boundary_gap(score,n)
M = numel(score);
if n<=0 || n>=M
    gap = NaN;
    return;
end
s = sort(score(:),'descend');
gap = s(n)-s(n+1);
end

function [mask,n0,nstar,eta] = select_mask_full_rule(po,P)
M = numel(po);

fdhat = mean(po(:) >= P.tauDetect);

countHigh = nnz(po(:) >= P.tauHigh);
n0 = min(round(M*min(fdhat,P.fMax)),countHigh);
n0 = max(0,min(M,n0));

[~,ord] = sort(po(:),'descend');

nstar = n0;
while nstar > 0
    eta = two_level_ess(nstar,M,P.epsSel);
    if eta >= P.etaMin
        break;
    end
    nstar = nstar-1;
end

if nstar==0
    eta = 1;
end

mask = false(size(po));
if nstar>0
    mask(ord(1:nstar)) = true;
end
end

function eta = two_level_ess(n,M,epsFloor)
f = n/M;
eta = (1-f+epsFloor*f)^2 / ...
      (1-f+(epsFloor^2)*f);
end
