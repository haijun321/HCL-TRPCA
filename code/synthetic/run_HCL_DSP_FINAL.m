function R = run_HCL_DSP_FINAL(mode)
% RUN_HCL_DSP_FINAL
%
% Recommended order:
%   run_HCL_DSP_FINAL('SMOKE')
%   run_HCL_DSP_FINAL('CONFIRMATION')
%   run_HCL_DSP_FINAL('GATE')
%   run_HCL_DSP_FINAL('WEAK_OVERLAP')
%   run_HCL_DSP_FINAL('RANK')
%
% Output:
%   generated_results/<experiment>/

if nargin<1, mode='SMOKE'; end
C = HCL_DSP_FINAL_config();
assert_final_protocol(C);

if ~exist(C.outputRoot,'dir'), mkdir(C.outputRoot); end
save(fullfile(C.outputRoot,'FINAL_PROTOCOL.mat'),'C');

switch upper(mode)
    case 'SMOKE'
        R = exp_smoke(C);
    case 'CONFIRMATION'
        R = exp_confirmation(C);
    case 'GATE'
        R = exp_gate(C);
    case 'WEAK_OVERLAP'
        R = exp_weak_overlap(C);
    case 'RANK'
        R = exp_rank(C);
    case 'ALL_SYNTHETIC'
        R.smoke = exp_smoke(C);
        R.confirmation = exp_confirmation(C);
        R.gate = exp_gate(C);
        R.weak_overlap = exp_weak_overlap(C);
        R.rank = exp_rank(C);
    otherwise
        error('Unknown mode: %s',mode);
end
end

% =========================================================================
function R = exp_smoke(C)
fprintf('\n=== HCL DSP FINAL : SMOKE ===\n');
seed=C.confirmationSeeds(1);
D=hcl_synthetic_case(seed,'EBS',C,'DISJOINT');

H=hcl_trpca_final(D.Y,D.rank,C.hcl);
M=hcl_metrics(H.Xhat,D,H);

fprintf('seed=%d | NRE=%.6g | DetF1=%.4f | MaskF1=%.4f | LevelF1=%.4f\n', ...
    seed,M.NRE,M.DetF1,M.MaskF1,M.LevelMacroF1);
fprintf('selected=%.4f | rESS=%.4f | it=%d | rel=%.3e\n', ...
    H.selectedFraction,H.recoveryESS,H.iterations,H.finalRelChange);

% Algorithm-identity checks
assert(all(abs(H.W0(H.mask)-C.hcl.epsRec)<1e-12), ...
    'Frozen mask does not use epsRec.');
assert(all(abs(H.W0(~H.mask)-1)<1e-12), ...
    'Unselected positions must have recovery weight 1.');

R=struct('D',D,'H',H,'M',M);
save(fullfile(C.outputRoot,'SMOKE_result.mat'),'R','-v7.3');
end

% =========================================================================
function R = exp_confirmation(C)
fprintf('\n=== HCL DSP FINAL : INDEPENDENT CONFIRMATION ===\n');
outDir=fullfile(C.outputRoot,'E01_CONFIRMATION');
if ~exist(outDir,'dir'),mkdir(outDir);end

scenarios={'EB','ES','BS','EBS'};
rows={};
rowId=0;

for iseed=1:numel(C.confirmationSeeds)
    seed=C.confirmationSeeds(iseed);
    for is=1:numel(scenarios)
        sc=scenarios{is};
        fprintf('[CONF] seed %d | %s\n',seed,sc);

        D=hcl_synthetic_case(seed,sc,C,'DISJOINT');
        hashY = cheap_hash(D.Y);
        hashX = cheap_hash(D.Xstar);

        % 1) tSVD-LRA
        X0=hcl_project_tubal_rank(D.Y,D.rank);
        rowId=rowId+1;
        rows(rowId,:)=make_row(seed,sc,'tSVD-LRA', ...
            norm(X0(:)-D.Xstar(:))/norm(D.Xstar(:)), ...
            NaN,NaN,NaN,NaN,NaN,hashY,hashX,C);

        % 2) External TNN if registered
        if ~isempty(C.baselines.TNN)
            O=C.baselines.TNN(D.Y,D.Xstar,D.truth,C);
            assert(isfield(O,'Xhat'),'TNN adapter must return out.Xhat');
            nre=norm(O.Xhat(:)-D.Xstar(:))/norm(D.Xstar(:));
            rowId=rowId+1;
            rows(rowId,:)=make_row(seed,sc,'TNN-TRPCA',nre, ...
                NaN,NaN,NaN,NaN,NaN,hashY,hashX,C);
        end

        % 3) External p-TRPCA if registered
        if ~isempty(C.baselines.pTRPCA)
            O=C.baselines.pTRPCA(D.Y,D.Xstar,D.truth,C);
            assert(isfield(O,'Xhat'),'pTRPCA adapter must return out.Xhat');
            nre=norm(O.Xhat(:)-D.Xstar(:))/norm(D.Xstar(:));
            rowId=rowId+1;
            rows(rowId,:)=make_row(seed,sc,'p-TRPCA',nre, ...
                NaN,NaN,NaN,NaN,NaN,hashY,hashX,C);
        end

        % 4) External Entry-Robust if registered
        if ~isempty(C.baselines.EntryRobust)
            O=C.baselines.EntryRobust(D.Y,D.Xstar,D.truth,C);
            assert(isfield(O,'Xhat'),'EntryRobust adapter must return out.Xhat');
            nre=norm(O.Xhat(:)-D.Xstar(:))/norm(D.Xstar(:));
            rowId=rowId+1;
            rows(rowId,:)=make_row(seed,sc,'Entry-Robust',nre, ...
                NaN,NaN,NaN,NaN,NaN,hashY,hashX,C);
        end

        % 5) Final HCL-TRPCA
        H=hcl_trpca_final(D.Y,D.rank,C.hcl);
        M=hcl_metrics(H.Xhat,D,H);
        rowId=rowId+1;
        rows(rowId,:)=make_row(seed,sc,'HCL-TRPCA', ...
            M.NRE,M.DetF1,M.MaskF1,M.LevelMacroF1, ...
            H.recoveryESS,H.selectedFraction,hashY,hashX,C);

        % 6) True-mask diagnostic oracle with SAME fixed recovery operator
        Xoracle=oracle_fixed_recovery(D.Y,D.truth.unionMask,D.rank,C.hcl);
        nre=norm(Xoracle(:)-D.Xstar(:))/norm(D.Xstar(:));
        rowId=rowId+1;
        rows(rowId,:)=make_row(seed,sc,'Oracle',nre, ...
            NaN,NaN,NaN,NaN,mean(D.truth.unionMask(:)),hashY,hashX,C);

        % Save HCL state for zero-cost gate ablation
        if C.flags.saveFullSyntheticState
            runDir=fullfile(outDir,sprintf('seed_%d',seed),sc);
            if ~exist(runDir,'dir'),mkdir(runDir);end
            save(fullfile(runDir,'HCL_state.mat'),'D','H','M','-v7.3');
        end
    end
end

T=cell2table(rows,'VariableNames', ...
    {'Seed','Scenario','Method','NRE','DetectionF1','RecoveryMaskF1', ...
     'LevelMacroF1','rESS','SelectedFraction','InputHash','ReferenceHash', ...
     'EpsRec'});
writetable(T,fullfile(outDir,'confirmation_raw.csv'));

S=summarize_nre(T);
writetable(S,fullfile(outDir,'confirmation_summary.csv'));

R=struct('raw',T,'summary',S);
save(fullfile(outDir,'confirmation_results.mat'),'R','-v7.3');
end

% =========================================================================
function R = exp_gate(C)
fprintf('\n=== HCL DSP FINAL : GATE ABLATION ===\n');
confDir=fullfile(C.outputRoot,'E01_CONFIRMATION');
outDir=fullfile(C.outputRoot,'E02_GATE');
if ~exist(outDir,'dir'),mkdir(outDir);end

scenarios={'EB','ES','BS','EBS'};
allT=table();

for iseed=1:numel(C.confirmationSeeds)
    seed=C.confirmationSeeds(iseed);
    for is=1:numel(scenarios)
        sc=scenarios{is};
        f=fullfile(confDir,sprintf('seed_%d',seed),sc,'HCL_state.mat');
        if ~exist(f,'file')
            error('Missing confirmation cache: %s. Run CONFIRMATION first.',f);
        end
        L=load(f,'D','H');
        T=hcl_gate_ablation(L.H,L.D);
        T.Seed=repmat(seed,height(T),1);
        T.Scenario=repmat(string(sc),height(T),1);
        allT=[allT;T]; %#ok<AGROW>
    end
end

writetable(allT,fullfile(outDir,'gate_raw.csv'));

G=groupsummary(allT,{'Scenario','Variant'}, ...
    {'mean','std'},{'LevelMacroF1','EntryF1','BlockF1','SliceF1'});
writetable(G,fullfile(outDir,'gate_summary.csv'));

R=struct('raw',allT,'summary',G);
save(fullfile(outDir,'gate_results.mat'),'R');
end

% =========================================================================
function R = exp_weak_overlap(C)
fprintf('\n=== HCL DSP FINAL : WEAK + OVERLAP ===\n');
outDir=fullfile(C.outputRoot,'E03_WEAK_OVERLAP');
if ~exist(outDir,'dir'),mkdir(outDir);end

rows={}; rid=0;
for seed=C.robustnessSeeds
    % Weak disjoint EBS
    D=hcl_synthetic_case(seed,'EBS',C,'WEAK');
    H=hcl_trpca_final(D.Y,D.rank,C.hcl);
    M=hcl_metrics(H.Xhat,D,H);
    rid=rid+1;
    rows(rid,:)={seed,'WEAK_EBS',M.NRE,M.DetF1,M.MaskF1, ...
                 M.LevelMacroF1,NaN,H.selectedFraction};

    % Overlapping EBS
    D=hcl_synthetic_case(seed,'EBS',C,'OVERLAP');
    H=hcl_trpca_final(D.Y,D.rank,C.hcl);
    M=hcl_metrics(H.Xhat,D,H);
    aam=active_attribution_mass(H.final,D.truth.activeLevels,D.truth.unionMask);
    rid=rid+1;
    rows(rid,:)={seed,'OVERLAP_EBS',M.NRE,M.DetF1,M.MaskF1, ...
                 NaN,aam,H.selectedFraction};
end

T=cell2table(rows,'VariableNames', ...
    {'Seed','Scenario','NRE','DetectionF1','RecoveryMaskF1', ...
     'LevelMacroF1','ActiveAttributionMass','SelectedFraction'});
writetable(T,fullfile(outDir,'weak_overlap_raw.csv'));

G=groupsummary(T,'Scenario',{'mean','std'}, ...
    {'NRE','DetectionF1','RecoveryMaskF1','LevelMacroF1', ...
     'ActiveAttributionMass','SelectedFraction'});
writetable(G,fullfile(outDir,'weak_overlap_summary.csv'));

R=struct('raw',T,'summary',G);
save(fullfile(outDir,'weak_overlap_results.mat'),'R');
end

% =========================================================================
function R = exp_rank(C)
fprintf('\n=== HCL DSP FINAL : RANK MISSPECIFICATION ===\n');
outDir=fullfile(C.outputRoot,'E04_RANK');
if ~exist(outDir,'dir'),mkdir(outDir);end

rows={}; rid=0;
for seed=C.failureSeeds
    D=hcl_synthetic_case(seed,'EBS',C,'DISJOINT');
    for r=C.rankGrid
        H=hcl_trpca_final(D.Y,r,C.hcl);
        M=hcl_metrics(H.Xhat,D,H);
        rid=rid+1;
        rows(rid,:)={seed,r,M.NRE,M.DetF1,M.MaskF1,M.LevelMacroF1, ...
                     H.iterations,H.finalRelChange,H.toleranceMet};
    end
end

T=cell2table(rows,'VariableNames', ...
    {'Seed','FittedRank','NRE','DetectionF1','RecoveryMaskF1', ...
     'LevelMacroF1','Iterations','FinalRelChange','ToleranceMet'});
writetable(T,fullfile(outDir,'rank_raw.csv'));

G=groupsummary(T,'FittedRank',{'mean','std'}, ...
    {'NRE','DetectionF1','RecoveryMaskF1','LevelMacroF1','Iterations'});
writetable(G,fullfile(outDir,'rank_summary.csv'));

R=struct('raw',T,'summary',G);
save(fullfile(outDir,'rank_results.mat'),'R');
end

% =========================================================================
function X = oracle_fixed_recovery(Y,trueMask,r,P)
X=hcl_project_tubal_rank(Y,r);
W=ones(size(Y));
W(trueMask)=P.epsRec;
for t=1:P.maxIterations
    V=W.*Y+(1-W).*X;
    Z=hcl_project_tubal_rank(V,r);
    Xnew=P.gamma*Z+(1-P.gamma)*X;
    rel=norm(Xnew(:)-X(:))/max(norm(X(:)),eps);
    X=Xnew;
    if t>=P.warmIterations && rel<P.relTol, break; end
end
end

function aam=active_attribution_mass(F,active,unionMask)
Q=cat(4,F.qe,F.qb,F.qs);
s=sum(Q.*double(active),4);
v=s(unionMask);
aam=mean(v);
end

function row=make_row(seed,sc,method,nre,df1,mf1,lf1,ress,sel,hY,hX,C)
row={seed,sc,method,nre,df1,mf1,lf1,ress,sel,hY,hX,C.hcl.epsRec};
end

function S=summarize_nre(T)
methods=unique(T.Method,'stable');
scens=unique(T.Scenario,'stable');
rows={}; k=0;
for i=1:numel(methods)
    for j=1:numel(scens)
        idx=strcmp(T.Method,methods{i}) & strcmp(T.Scenario,scens{j});
        if any(idx)
            x=T.NRE(idx);
            k=k+1;
            rows(k,:)={methods{i},scens{j},mean(x,'omitnan'),std(x,'omitnan'),sum(~isnan(x))};
        end
    end
end
S=cell2table(rows,'VariableNames',{'Method','Scenario','MeanNRE','SDNRE','N'});
end

function h=cheap_hash(X)
% Pairing fingerprint, not a cryptographic hash.
v=double(X(:));
h=sprintf('%.12e_%.12e_%.12e_%d',sum(v),sum(v.^2),sum(abs(v)),numel(v));
end

function assert_final_protocol(C)
assert(abs(C.hcl.epsWarm-1e-2)<1e-14,'epsWarm changed.');
assert(abs(C.hcl.epsSel -1e-2)<1e-14,'epsSel changed.');
assert(abs(C.hcl.epsRec -1e-3)<1e-14,'epsRec changed.');
assert(C.hcl.b==8,'b changed.');
assert(abs(C.hcl.tauHigh-0.90)<1e-14,'tauHigh changed.');
assert(abs(C.hcl.fMax-0.32)<1e-14,'fMax changed.');
assert(abs(C.hcl.etaMin-0.68)<1e-14,'etaMin changed.');
assert(abs(C.hcl.gamma-0.85)<1e-14,'gamma changed.');
assert(abs(C.hcl.relTol-1e-5)<1e-14,'relTol changed.');
assert(isempty(intersect(C.confirmationSeeds,C.developmentSeeds)), ...
    'Development and confirmation seeds overlap.');
end
