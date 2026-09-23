function R = reproduce_confirmation_paired_statistics_from_archived_pairs(repoRoot)
% REPRODUCE_CONFIRMATION_PAIRED_STATISTICS_FROM_ARCHIVED_PAIRS
% Release-side reproduction/audit utility for the archived HCL-TRPCA
% confirmation paired rows.
%
% IMPORTANT
%   This is a clean v1.0.2 release utility reconstructed from the frozen
%   statistical protocol and archived paired rows. It is NOT represented as
%   the exact legacy POSTPROCESS_CONFIRMATION_PAIRED_STATS.m source file.
%
% Protocol frozen in the archived audit:
%   - scenarios: EB, ES, BS, EBS
%   - n = 30 paired tensors per scenario
%   - DeltaNRE = NRE_HCL - NRE_pTRPCA
%   - 10,000 paired percentile-bootstrap resamples
%   - bootstrap seed = 2026091601
%   - two-sided Wilcoxon signed-rank test
%   - Holm correction over the four scenario-wise tests
%
% Usage:
%   R = reproduce_confirmation_paired_statistics_from_archived_pairs();
%   R = reproduce_confirmation_paired_statistics_from_archived_pairs(repoRoot);

if nargin < 1 || isempty(repoRoot)
    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(thisFile));
end

statsDir = fullfile(repoRoot,'results','statistics');
pairFile = fullfile(statsDir,'confirmation_paired_differences.csv');
refFile  = fullfile(statsDir,'confirmation_paired_stats.csv');
assert(isfile(pairFile),'Missing paired file: %s',pairFile);
assert(isfile(refFile),'Missing frozen summary: %s',refFile);

T = readtable(pairFile,'TextType','string');
Ref = readtable(refFile,'TextType','string');
scenarios = ["EB","ES","BS","EBS"];
B = 10000;
rng(2026091601,'twister');

rows = cell(numel(scenarios),15);
rawP = nan(numel(scenarios),1);

for k=1:numel(scenarios)
    s = scenarios(k);
    A = T(upper(T.Scenario)==s,:);
    assert(height(A)==30,'Expected 30 pairs for %s; found %d.',s,height(A));
    d = A.DeltaNRE(:);
    n = numel(d);

    idx = randi(n,n,B);
    bootMeans = mean(d(idx),1);
    ci = prctile(bootMeans,[2.5 97.5]);

    if exist('signrank','file')==2
        rawP(k)=signrank(d,0,'tail','both');
    else
        warning('signrank not available; Wilcoxon field set to NaN.');
    end

    rows{k,1}=s;
    rows{k,2}=n;
    rows{k,3}=mean(A.HCL_NRE);
    rows{k,4}=std(A.HCL_NRE,0);
    rows{k,5}=mean(A.pTRPCA_NRE);
    rows{k,6}=std(A.pTRPCA_NRE,0);
    rows{k,7}=mean(d<0);
    rows{k,8}=mean(d>0);
    rows{k,9}=mean(d==0);
    rows{k,10}=mean(d);
    rows{k,11}=median(d);
    rows{k,12}=ci(1);
    rows{k,13}=ci(2);
    rows{k,14}=rawP(k);
    rows{k,15}=NaN;
end

if all(isfinite(rawP))
    hp = holm_adjust(rawP);
    for k=1:numel(scenarios), rows{k,15}=hp(k); end
end

Stats = cell2table(rows,'VariableNames',{...
    'Scenario','N','Mean_HCL_NRE','SD_HCL_NRE','Mean_pTRPCA_NRE',...
    'SD_pTRPCA_NRE','WinRate','LossRate','TieRate','Mean_DeltaNRE',...
    'Median_DeltaNRE','Bootstrap95CI_Lower','Bootstrap95CI_Upper',...
    'WilcoxonP','HolmAdjustedP'});

% Deterministic quantities must match the frozen archive essentially exactly.
fieldsExact = {'N','Mean_HCL_NRE','SD_HCL_NRE','Mean_pTRPCA_NRE',...
    'SD_pTRPCA_NRE','WinRate','LossRate','TieRate','Mean_DeltaNRE',...
    'Median_DeltaNRE'};
tol = 5e-12;
for k=1:numel(scenarios)
    s=scenarios(k);
    j=find(upper(Ref.Scenario)==s,1);
    assert(~isempty(j),'Missing frozen reference row for %s.',s);
    for q=1:numel(fieldsExact)
        f=fieldsExact{q};
        a=Stats.(f)(k); b=Ref.(f)(j);
        assert(abs(double(a)-double(b))<tol, ...
            'Mismatch in %s/%s: reproduced %.15g vs frozen %.15g',s,f,a,b);
    end
end

fprintf('\n============================================================\n');
fprintf('ARCHIVED CONFIRMATION PAIRED-STATISTICS REPRODUCTION\n');
fprintf('============================================================\n');
disp(Stats);

% The exact legacy RNG implementation is not claimed here. Show any CI/test
% deviations relative to the frozen archive rather than overwriting it.
Compare = table(strings(numel(scenarios),1),nan(numel(scenarios),1),...
    nan(numel(scenarios),1),nan(numel(scenarios),1),nan(numel(scenarios),1),...
    'VariableNames',{'Scenario','CILowerDiff','CIUpperDiff','WilcoxonPDiff','HolmPDiff'});
for k=1:numel(scenarios)
    s=scenarios(k); j=find(upper(Ref.Scenario)==s,1);
    Compare.Scenario(k)=s;
    Compare.CILowerDiff(k)=Stats.Bootstrap95CI_Lower(k)-Ref.Bootstrap95CI_Lower(j);
    Compare.CIUpperDiff(k)=Stats.Bootstrap95CI_Upper(k)-Ref.Bootstrap95CI_Upper(j);
    Compare.WilcoxonPDiff(k)=Stats.WilcoxonP(k)-Ref.WilcoxonP(j);
    Compare.HolmPDiff(k)=Stats.HolmAdjustedP(k)-Ref.HolmAdjustedP(j);
end
fprintf('\nDifferences from frozen archived inference fields:\n');
disp(Compare);

R = struct('stats',Stats,'frozenReference',Ref,'comparison',Compare);
end

function pAdj = holm_adjust(p)
p=p(:); m=numel(p);
[ps,ord]=sort(p,'ascend');
adj=zeros(m,1); running=0;
for i=1:m
    running=max(running,(m-i+1)*ps(i));
    adj(i)=min(1,running);
end
pAdj=nan(m,1); pAdj(ord)=adj;
end
