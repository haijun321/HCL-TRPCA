function R = run_HCL_ATTRIBUTION_REVIEW_UPGRADE_V1(mode)
%RUN_HCL_ATTRIBUTION_REVIEW_UPGRADE_V1
% Robust post-hoc attribution audit and reviewer-control analysis.
%
% Usage:
%   R = run_HCL_ATTRIBUTION_REVIEW_UPGRADE_V1('AUDIT');
%   R = run_HCL_ATTRIBUTION_REVIEW_UPGRADE_V1('ALL');
%
% Put this file in the HCL_DSP_FINAL_package root and run it there.  The
% script does not rerun HCL-TRPCA; it recursively reads existing MAT files.

if nargin < 1 || isempty(mode), mode = 'AUDIT'; end
mode = upper(char(string(mode)));
if ~ismember(mode,{'AUDIT','ALL'})
    error('Mode must be AUDIT or ALL.');
end

CFG = make_config();
if ~exist(CFG.outDir,'dir'), mkdir(CFG.outDir); end

fprintf('\n============================================================\n');
fprintf('HCL ATTRIBUTION REVIEW UPGRADE V1 FIXED\n');
fprintf('MODE : %s\nROOT : %s\n',mode,CFG.rootDir);
fprintf('============================================================\n\n');

fprintf('[1/7] Discovering candidate MAT files...\n');
FILES = discover_mat_files(CFG);
fprintf('Candidate MAT files found: %d\n',height(FILES));
if isempty(FILES), error('No MAT files found below: %s',CFG.rootDir); end

fprintf('[2/7] Auditing attribution states...\n');
[AUDIT,CASES] = audit_candidate_files(FILES,CFG);
writetable(AUDIT,fullfile(CFG.outDir,'AUDIT_inventory.csv'));

fprintf('\nAUDIT table variables:\n');
disp(AUDIT.Properties.VariableNames');
fprintf('\nUsable state counts:\n');
for n = {'DEV','CONFIRM','OVERLAP','PAVIA'}
    show_dataset_count(AUDIT,n{1});
end
print_audit_details(AUDIT);

R = struct('config',CFG,'audit',AUDIT,'cases',CASES);
save(fullfile(CFG.outDir,'AUDIT_workspace.mat'),'R','-v7.3');

if strcmp(mode,'AUDIT')
    fprintf('\nAUDIT ONLY complete. Inspect:\n%s\n', ...
        fullfile(CFG.outDir,'AUDIT_inventory.csv'));
    return;
end

fprintf('\n[3/7] Fitting calibrated nonhierarchical control on DEV...\n');
DEV = CASES(strcmp({CASES.dataset},'DEV') & [CASES.usableClass]);
if isempty(DEV)
    error(['No usable DEV cases. Inspect AUDIT_inventory.csv. ', ...
           'ALL requires ae/ab/as/po and a 0/1/2/3 GT label map.']);
end
if numel(DEV) ~= 20
    error('Protocol violation: expected 20 usable DEV states, found %d.',numel(DEV));
end
CAL = fit_calibrated_flat(DEV,CFG);
R.calibration = CAL;
writetable(struct2table(CAL,'AsArray',true), ...
    fullfile(CFG.outDir,'calibrated_flat_parameters.csv'));

fprintf('[4/7] Evaluating synthetic confirmation states...\n');
CONF = CASES(strcmp({CASES.dataset},'CONFIRM') & [CASES.usableClass]);
if isempty(CONF)
    warning('No usable CONFIRM cases; synthetic outputs are skipped.');
    SYN = struct('raw',table(),'summary',table(),'paired',table());
elseif numel(CONF) ~= 120
    error('Protocol violation: expected 120 usable CONFIRM states, found %d.',numel(CONF));
else
    SYN = evaluate_synthetic(CONF,CAL,CFG);
    writetable(SYN.raw,fullfile(CFG.outDir,'synthetic_gate_comparison_raw.csv'));
    writetable(SYN.summary,fullfile(CFG.outDir,'synthetic_gate_comparison_summary.csv'));
    writetable(SYN.paired,fullfile(CFG.outDir,'synthetic_gate_paired_HCL_vs_CalFlat.csv'));
    make_gate_figure(SYN,CFG);
end
R.synthetic = SYN;

fprintf('[5/7] Evaluating Pavia attribution states...\n');
PC = CASES(strcmp({CASES.dataset},'PAVIA') & [CASES.usableClass]);
if numel(PC) < 40
    warning('Pavia formal attribution control skipped: only %d/40 usable detailed states.',numel(PC));
    PAV = struct('raw',table(),'summary',table(),'confusion',struct());
else
    PAV = evaluate_pavia(PC,CAL,CFG);
    writetable(PAV.raw,fullfile(CFG.outDir,'pavia_attribution_metrics_raw.csv'));
    writetable(PAV.summary,fullfile(CFG.outDir,'pavia_attribution_metrics_summary.csv'));
    write_pavia_confusions(PAV,CFG);
    make_pavia_confusion_figure(PAV,CFG);
end
R.pavia = PAV;

fprintf('[6/7] Evaluating overlap active-attribution mass...\n');
OC = CASES(strcmp({CASES.dataset},'OVERLAP') & [CASES.usableOverlap]);
if isempty(OC)
    warning('No usable OVERLAP cases; overlap outputs are skipped.');
    OV = struct('raw',table(),'summary',table());
else
    OV = evaluate_overlap(OC,CFG);
    writetable(OV.raw,fullfile(CFG.outDir,'overlap_aam_stratified_raw.csv'));
    writetable(OV.summary,fullfile(CFG.outDir,'overlap_aam_stratified_summary.csv'));
    make_overlap_figure(OV,CFG);
end
R.overlap = OV;

fprintf('[7/7] Saving consolidated results...\n');
save(fullfile(CFG.outDir,'HCL_ATTRIBUTION_REVIEW_UPGRADE_V1_FIXED_RESULTS.mat'), ...
    'R','-v7.3');
fprintf('\nALL complete. Outputs: %s\n',CFG.outDir);
end

%% Configuration
function CFG = make_config()
CFG = struct;
thisFile = mfilename('fullpath');
CFG.rootDir = fileparts(thisFile);
CFG.outDir = fullfile(CFG.rootDir,'ATTRIBUTION_REVIEW_UPGRADE_V1');
CFG.devSeeds = 20261001:20261005;

CFG.confirmSeeds = 20270001:20270030;

CFG.overlapSeeds = [];

CFG.paviaSeeds = 20280001:20280010;

CFG.syntheticScenarios = {'EB','ES','BS','EBS'};
CFG.detectThreshold = 0.50;
CFG.calibrationMaxPerClass = 200000;
CFG.bootstrapB = 5000;
CFG.bootstrapSeed = 93021;
end

%% Discovery and audit
function T = discover_mat_files(CFG)
D = dir(fullfile(CFG.rootDir,'**','*.mat'));
keep = ~[D.isdir];
D = D(keep);
paths = strings(numel(D),1);
for i=1:numel(D), paths(i)=string(fullfile(D(i).folder,D(i).name)); end
bad = contains(paths,string(CFG.outDir),'IgnoreCase',true);
paths = paths(~bad);
T = table(cellstr(paths),'VariableNames',{'File'});
end

function [A,CASES] = audit_candidate_files(FILES,CFG)
n = height(FILES);
File=string(FILES.File); Dataset=strings(n,1); Scenario=strings(n,1);
Seed=nan(n,1); HasAE=false(n,1); HasAB=false(n,1); HasAS=false(n,1);
HasPO=false(n,1); HasQE=false(n,1); HasQB=false(n,1); HasQS=false(n,1);
HasGTLabel=false(n,1); HasGTEntry=false(n,1); HasGTBlock=false(n,1);
HasGTSlice=false(n,1); Usable=false(n,1); UsableClass=false(n,1);
UsableOverlap=false(n,1); Reason=strings(n,1);
emptyCase = blank_case(); CASES = repmat(emptyCase,n,1);

for i=1:n
    f=char(File(i)); c=emptyCase(); c.file=f;
    try
        S=load(f);
    catch ME
        Reason(i)="LOAD_FAILED: "+string(ME.identifier); CASES(i)=c; continue;
    end
    Seed(i)=parse_seed(f,S,CFG); Dataset(i)=classify_dataset(Seed(i),f,CFG);
    Scenario(i)=classify_scenario(f,S,Dataset(i));
    c.seed=Seed(i); c.dataset=char(Dataset(i)); c.scenario=char(Scenario(i));

    c.ae=force_double(recursive_find(S,alias_ae()));
    c.ab=force_double(recursive_find(S,alias_ab()));
    c.as=force_double(recursive_find(S,alias_as()));
    c.po=force_double(recursive_find(S,alias_po()));
    c.qe=force_double(recursive_find(S,alias_qe()));
    c.qb=force_double(recursive_find(S,alias_qb()));
    c.qs=force_double(recursive_find(S,alias_qs()));
    L=recursive_find(S,alias_gtlabel());
    c.gtEntry=logical_or_empty(recursive_find(S,alias_gtentry()));
    c.gtBlock=logical_or_empty(recursive_find(S,alias_gtblock()));
    c.gtSlice=logical_or_empty(recursive_find(S,alias_gtslice()));

    HasAE(i)=~isempty(c.ae); HasAB(i)=~isempty(c.ab); HasAS(i)=~isempty(c.as);
    HasPO(i)=~isempty(c.po); HasQE(i)=~isempty(c.qe); HasQB(i)=~isempty(c.qb);
    HasQS(i)=~isempty(c.qs); HasGTEntry(i)=~isempty(c.gtEntry);
    HasGTBlock(i)=~isempty(c.gtBlock); HasGTSlice(i)=~isempty(c.gtSlice);
    if ~isempty(L)
        try, c.gtLabel=normalize_label_map(L); catch, c.gtLabel=[]; end
    elseif all([HasGTEntry(i),HasGTBlock(i),HasGTSlice(i)])
        c.gtLabel=label_from_masks(c.gtEntry,c.gtBlock,c.gtSlice);
    end
    HasGTLabel(i)=~isempty(c.gtLabel);

    base = all([HasAE(i),HasAB(i),HasAS(i),HasPO(i)]);
    qok = all([HasQE(i),HasQB(i),HasQS(i)]);
    shapes = compatible_shapes(c);
    UsableClass(i)=base && qok && HasGTLabel(i) && shapes;
    UsableOverlap(i)=qok && all([HasGTEntry(i),HasGTBlock(i),HasGTSlice(i)]) && shapes;
    Usable(i)=UsableClass(i) || UsableOverlap(i);
    c.usableClass=UsableClass(i); c.usableOverlap=UsableOverlap(i);
    if Usable(i), Reason(i)="OK"; else, Reason(i)=missing_reason(i); end
    CASES(i)=c;
end

% Explicit names prevent version-dependent table-name inference.
A = table(File,Dataset,Scenario,Seed,HasAE,HasAB,HasAS,HasPO, ...
    HasQE,HasQB,HasQS,HasGTLabel,HasGTEntry,HasGTBlock,HasGTSlice, ...
    UsableClass,UsableOverlap,Usable,Reason, ...
    'VariableNames',{'File','Dataset','Scenario','Seed','HasAE','HasAB', ...
    'HasAS','HasPO','HasQE','HasQB','HasQS','HasGTLabel','HasGTEntry', ...
    'HasGTBlock','HasGTSlice','UsableClass','UsableOverlap','Usable','Reason'});

    function r=missing_reason(k)
        miss={};
        flags=[HasAE(k),HasAB(k),HasAS(k),HasPO(k),HasQE(k),HasQB(k),HasQS(k),HasGTLabel(k)];
        names={'ae','ab','as','po','qe','qb','qs','GTLabel'};
        for z=1:numel(flags), if ~flags(z), miss{end+1}=names{z}; end, end %#ok<AGROW>
        if ~shapes, miss{end+1}='shapeMismatch'; end
        r="MISSING: "+string(strjoin(miss,','));
    end
end

function c=blank_case()
c=struct('file','','dataset','OTHER','scenario','','seed',NaN, ...
 'ae',[],'ab',[],'as',[],'po',[],'qe',[],'qb',[],'qs',[], ...
 'gtLabel',[],'gtEntry',[],'gtBlock',[],'gtSlice',[], ...
 'usableClass',false,'usableOverlap',false);
end

function tf=compatible_shapes(c)
arr={c.ae,c.ab,c.as,c.po,c.qe,c.qb,c.qs,c.gtLabel,c.gtEntry,c.gtBlock,c.gtSlice};
sz=[]; tf=true;
for i=1:numel(arr)
    if isempty(arr{i}), continue; end
    if isempty(sz), sz=size(arr{i}); elseif ~isequal(size(arr{i}),sz), tf=false; return; end
end
end

function show_dataset_count(A,name)
required={'Dataset','Usable'};
for k=1:numel(required)
    if ~ismember(required{k},A.Properties.VariableNames)
        error('AUDIT table is missing required variable: %s',required{k});
    end
end
idx=string(A.Dataset)==string(name) & logical(A.Usable);
fprintf('  %-8s : %d usable files\n',char(string(name)),sum(idx));
end

function print_audit_details(A)
fprintf('\n============================================================\n');
fprintf('AUDIT DETAILS BY DATASET\n');
fprintf('============================================================\n');
for d=["DEV","CONFIRM","OVERLAP","PAVIA"]
 idx=string(A.Dataset)==d; fprintf('\n[%s]\n',d);
 fprintf('Total candidate files : %d\n',sum(idx));
 fprintf('Usable files          : %d\n',sum(idx&A.Usable));
 for v={'HasAE','HasAB','HasAS','HasPO','HasQE','HasQB','HasQS','HasGTLabel','HasGTEntry','HasGTBlock','HasGTSlice'}
   fprintf('%-22s: %d\n',v{1},sum(idx&A.(v{1})));
 end
end
end

%% Calibration
function CAL = fit_calibrated_flat(CASES,CFG)
% Calibration uses FINAL HCL evidence only.
% DEV labels are used exclusively here.

Xall = [];
yall = [];

fprintf('\nStrict DEV final-state loading...\n');

for i = 1:numel(CASES)

    X = load_synthetic_final_state( ...
        CASES(i).file);

    idx = ...
        X.gtLabel > 0 & ...
        isfinite(X.ae) & ...
        isfinite(X.ab) & ...
        isfinite(X.as);

    if ~any(idx(:))
        continue;
    end

    Xi = [ ...
        X.ae(idx), ...
        X.ab(idx), ...
        X.as(idx)];

    yi = double(X.gtLabel(idx));

    Xall = [Xall;Xi]; %#ok<AGROW>
    yall = [yall;yi]; %#ok<AGROW>

end


if isempty(Xall)

    error('DEV contains no labeled abnormal positions.');

end


%% Balanced calibration subsampling

rng(CFG.bootstrapSeed,'twister');

keep = [];

for k = 1:3

    q = find(yall==k);

    if numel(q) > CFG.calibrationMaxPerClass

        q = q(randperm( ...
            numel(q), ...
            CFG.calibrationMaxPerClass));

    end

    keep = [keep;q]; %#ok<AGROW>

end

Xall = Xall(keep,:);
yall = yall(keep);


%% Log transform and standardization

Z = log(max(Xall,eps));

mu = mean(Z,1);

sg = std(Z,0,1);

sg(sg<eps) = 1;

Z = (Z-mu)./sg;


%% Calibrated flat optimization

theta0 = [ ...
    ones(1,3), ...
    zeros(1,3)];

opt = optimset( ...
    'Display','off', ...
    'MaxIter',2000, ...
    'MaxFunEvals',10000, ...
    'TolX',1e-8, ...
    'TolFun',1e-8);

theta = fminsearch( ...
    @(t) calibration_loss(t,Z,yall), ...
    theta0, ...
    opt);


CAL = struct( ...
    'slopeEntry',theta(1), ...
    'slopeBlock',theta(2), ...
    'slopeSlice',theta(3), ...
    'biasEntry',theta(4), ...
    'biasBlock',theta(5), ...
    'biasSlice',theta(6), ...
    'muEntry',mu(1), ...
    'muBlock',mu(2), ...
    'muSlice',mu(3), ...
    'sdEntry',sg(1), ...
    'sdBlock',sg(2), ...
    'sdSlice',sg(3), ...
    'NFit',numel(yall));

fprintf('Strict DEV calibration complete.\n');
fprintf('NFit = %d\n',CAL.NFit);

end

function L=calibration_loss(t,Z,y)
score=Z.*t(1:3)+t(4:6); score=score-max(score,[],2);
P=exp(score); P=P./sum(P,2); ind=sub2ind(size(P),(1:numel(y))',y);
L=-mean(log(max(P(ind),realmin)))+1e-4*sum((t(1:3)-1).^2)+1e-5*sum(t(4:6).^2);
end

function [qe,qb,qs]=apply_calibrated_flat(ae,ab,as,CAL)
sz=size(ae); Z=[log(max(ae(:),eps)),log(max(ab(:),eps)),log(max(as(:),eps))];
mu=[CAL.muEntry,CAL.muBlock,CAL.muSlice]; sd=[CAL.sdEntry,CAL.sdBlock,CAL.sdSlice];
Z=(Z-mu)./sd; t=[CAL.slopeEntry,CAL.slopeBlock,CAL.slopeSlice];
b=[CAL.biasEntry,CAL.biasBlock,CAL.biasSlice]; score=Z.*t+b; score=score-max(score,[],2);
P=exp(score); P=P./sum(P,2); qe=reshape(P(:,1),sz); qb=reshape(P(:,2),sz); qs=reshape(P(:,3),sz);
end

function OUT=evaluate_synthetic(CASES,CAL,CFG)

raw = evaluate_synthetic_strict_cases( ...
    CASES,CAL,CFG);

summary = summarize_metrics(raw);

rows = {};

sc = unique(raw.Scenario,'stable');

for s = 1:numel(sc)

    H = raw( ...
        strcmp(raw.Scenario,sc{s}) & ...
        strcmp(raw.Method,'FullHCL'),:);

    C = raw( ...
        strcmp(raw.Scenario,sc{s}) & ...
        strcmp(raw.Method,'CalibratedFlat'),:);

    seeds = intersect(H.Seed,C.Seed);

    d = [];

    for k = 1:numel(seeds)

        h = H.GranularityMacroF1( ...
            H.Seed==seeds(k));

        c = C.GranularityMacroF1( ...
            C.Seed==seeds(k));

        if ...
            ~isempty(h) && ...
            ~isempty(c) && ...
            isfinite(h(1)) && ...
            isfinite(c(1))

            d(end+1,1) = ...
                h(1)-c(1); %#ok<AGROW>

        end
    end

    [lo,hi] = bootstrap_ci( ...
        d, ...
        CFG.bootstrapB, ...
        CFG.bootstrapSeed+s);

    rows(end+1,:) = { ...
        sc{s}, ...
        numel(d), ...
        mean(d,'omitnan'), ...
        median(d,'omitnan'), ...
        mean(d>0), ...
        lo, ...
        hi}; %#ok<AGROW>

end


paired = cell2table( ...
    rows, ...
    'VariableNames',{ ...
    'Scenario', ...
    'NPaired', ...
    'MeanDeltaF1_HCLminusCalFlat', ...
    'MedianDeltaF1_HCLminusCalFlat', ...
    'HCLWinRate', ...
    'BootstrapCI_Low', ...
    'BootstrapCI_High'});


OUT = struct( ...
    'raw',raw, ...
    'summary',summary, ...
    'paired',paired);

end
function raw=evaluate_classification_cases(CASES,CAL,CFG)
rows={}; methods={'DirectFlat','CalibratedFlat','FullHCL'};
for i=1:numel(CASES)
 C=CASES(i);
 for m=1:numel(methods)
  switch methods{m}
   case 'DirectFlat', den=C.ae+C.ab+C.as+eps; qe=C.ae./den; qb=C.ab./den; qs=C.as./den;
   case 'CalibratedFlat', [qe,qb,qs]=apply_calibrated_flat(C.ae,C.ab,C.as,CAL);
   otherwise, qe=C.qe; qb=C.qb; qs=C.qs;
  end
  pred=attribution_label(C.po,qe,qb,qs,CFG.detectThreshold); M=classification_metrics(C.gtLabel,pred);
  rows(end+1,:)={C.scenario,C.seed,methods{m},M.DetectionF1,M.GranularityMacroF1,M.CleanFPR, ...
   M.PrecEntry,M.RecEntry,M.F1Entry,M.PrecBlock,M.RecBlock,M.F1Block,M.PrecSlice,M.RecSlice,M.F1Slice}; %#ok<AGROW>
 end
end
raw=cell2table(rows,'VariableNames',{'Scenario','Seed','Method','DetectionF1','GranularityMacroF1','CleanFPR','PrecisionEntry','RecallEntry','F1Entry','PrecisionBlock','RecallBlock','F1Block','PrecisionSlice','RecallSlice','F1Slice'});
end

function S=summarize_metrics(raw)
sc=unique(raw.Scenario,'stable'); me=unique(raw.Method,'stable'); rows={};
for i=1:numel(sc), for j=1:numel(me)
 idx=strcmp(raw.Scenario,sc{i})&strcmp(raw.Method,me{j});
 rows(end+1,:)={sc{i},me{j},sum(idx),mean(raw.DetectionF1(idx),'omitnan'),std(raw.DetectionF1(idx),0,'omitnan'),mean(raw.GranularityMacroF1(idx),'omitnan'),std(raw.GranularityMacroF1(idx),0,'omitnan'),mean(raw.CleanFPR(idx),'omitnan'),mean(raw.PrecisionEntry(idx),'omitnan'),mean(raw.RecallEntry(idx),'omitnan'),mean(raw.PrecisionBlock(idx),'omitnan'),mean(raw.RecallBlock(idx),'omitnan'),mean(raw.PrecisionSlice(idx),'omitnan'),mean(raw.RecallSlice(idx),'omitnan')}; %#ok<AGROW>
end, end
S=cell2table(rows,'VariableNames',{'Scenario','Method','N','DetectionF1_Mean','DetectionF1_SD','GranMacroF1_Mean','GranMacroF1_SD','CleanFPR_Mean','PrecisionEntry_Mean','RecallEntry_Mean','PrecisionBlock_Mean','RecallBlock_Mean','PrecisionSlice_Mean','RecallSlice_Mean'});
end

function [lo,hi]=bootstrap_ci(d,B,seed)
d=d(isfinite(d)); if isempty(d), lo=NaN; hi=NaN; return; end
rng(seed,'twister'); n=numel(d); v=zeros(B,1);
for b=1:B, v(b)=mean(d(randi(n,n,1))); end
q=prctile(v,[2.5 97.5]); lo=q(1); hi=q(2);
end

%% Pavia
function OUT=evaluate_pavia(CASES,CAL,CFG)
raw=evaluate_classification_cases(CASES,CAL,CFG); summary=summarize_metrics(raw); CONF=struct;
methods={'DirectFlat','CalibratedFlat','FullHCL'};
for i=1:numel(CASES), C=CASES(i);
 for m=1:numel(methods)
  switch methods{m}
   case 'DirectFlat', den=C.ae+C.ab+C.as+eps; qe=C.ae./den; qb=C.ab./den; qs=C.as./den;
   case 'CalibratedFlat', [qe,qb,qs]=apply_calibrated_flat(C.ae,C.ab,C.as,CAL);
   otherwise, qe=C.qe; qb=C.qb; qs=C.qs;
  end
  pred=attribution_label(C.po,qe,qb,qs,CFG.detectThreshold);
  key=matlab.lang.makeValidName(sprintf('%s_%s',C.scenario,methods{m}));
  if ~isfield(CONF,key), CONF.(key)=zeros(4); end
  CONF.(key)=CONF.(key)+confusion4(C.gtLabel,pred);
 end
end
OUT=struct('raw',raw,'summary',summary,'confusion',CONF);
end

function write_pavia_confusions(PAV,CFG)
names=fieldnames(PAV.confusion); cls={'Clean','Entry','Block','Slice'};
for i=1:numel(names)
 C=PAV.confusion.(names{i}); N=C./max(sum(C,2),1);
 writetable(array2table(C,'VariableNames',cls,'RowNames',cls),fullfile(CFG.outDir,['pavia_confusion_counts_' names{i} '.csv']),'WriteRowNames',true);
 writetable(array2table(N,'VariableNames',cls,'RowNames',cls),fullfile(CFG.outDir,['pavia_confusion_rownorm_' names{i} '.csv']),'WriteRowNames',true);
end
end

%% Overlap
function OUT=evaluate_overlap(CASES,~)
rows={};
for i=1:numel(CASES), C=CASES(i);
 active=double(C.gtEntry)+double(C.gtBlock)+double(C.gtSlice);
 local=C.qe.*double(C.gtEntry)+C.qb.*double(C.gtBlock)+C.qs.*double(C.gtSlice);
 omega=active>0; total=sum(omega(:)); if total==0, continue; end
 for k=1:3
  idx=active==k; nk=sum(idx(:));
  if nk==0, aam=NaN; frac=0; else, aam=mean(local(idx),'omitnan'); frac=nk/total; end
  rows(end+1,:)={C.seed,k,nk,frac,aam,k/3}; %#ok<AGROW>
 end
 rows(end+1,:)={C.seed,0,total,1,mean(local(omega),'omitnan'),mean(active(omega)/3,'omitnan')}; %#ok<AGROW>
end
raw=cell2table(rows,'VariableNames',{'Seed','ActiveMechanisms','NPositions','FractionOfCorruptedUnion','HCL_AAM','Uniform_AAM'});
sr={}; for k=[1 2 3 0]
 idx=raw.ActiveMechanisms==k;
 sr(end+1,:)={k,sum(idx),mean(raw.FractionOfCorruptedUnion(idx),'omitnan'),std(raw.FractionOfCorruptedUnion(idx),0,'omitnan'),mean(raw.HCL_AAM(idx),'omitnan'),std(raw.HCL_AAM(idx),0,'omitnan'),mean(raw.Uniform_AAM(idx),'omitnan'),std(raw.Uniform_AAM(idx),0,'omitnan')}; %#ok<AGROW>
end
summary=cell2table(sr,'VariableNames',{'ActiveMechanisms','NSeeds','FractionMean','FractionSD','HCL_AAM_Mean','HCL_AAM_SD','Uniform_AAM_Mean','Uniform_AAM_SD'});
OUT=struct('raw',raw,'summary',summary);
end

%% Classification
function pred=attribution_label(po,qe,qb,qs,tau)
pred=zeros(size(po),'uint8'); detected=po>=tau; Q=[qe(:),qb(:),qs(:)]; [~,lab]=max(Q,[],2);
lab=reshape(uint8(lab),size(po)); pred(detected)=lab(detected);
end

function M=classification_metrics(gt,pred)
gt=uint8(gt); pred=uint8(pred); C=confusion4(gt,pred);
t=gt~=0; p=pred~=0; M.DetectionF1=safe_f1(sum(t(:)&p(:)),sum(~t(:)&p(:)),sum(t(:)&~p(:)));
n=sum(gt(:)==0); if n==0, M.CleanFPR=NaN; else, M.CleanFPR=sum(gt(:)==0&pred(:)~=0)/n; end
P=nan(3,1); R=P; F=P;
for c=1:3
 tp=sum(gt(:)==c&pred(:)==c); fp=sum(gt(:)~=c&pred(:)==c); fn=sum(gt(:)==c&pred(:)~=c);
 P(c)=safe_div(tp,tp+fp); R(c)=safe_div(tp,tp+fn); F(c)=safe_f1(tp,fp,fn);
end
present=arrayfun(@(c)any(gt(:)==c),1:3); M.GranularityMacroF1=mean(F(present),'omitnan');
M.PrecEntry=P(1);M.PrecBlock=P(2);M.PrecSlice=P(3);M.RecEntry=R(1);M.RecBlock=R(2);M.RecSlice=R(3);M.F1Entry=F(1);M.F1Block=F(2);M.F1Slice=F(3);M.Confusion=C;
end

function C=confusion4(gt,pred)
C=zeros(4); for t=0:3, for p=0:3, C(t+1,p+1)=sum(gt(:)==t&pred(:)==p); end, end
end
function f=safe_f1(tp,fp,fn), d=2*tp+fp+fn; if d==0,f=NaN;else,f=2*tp/d;end, end
function x=safe_div(a,b), if b==0,x=NaN;else,x=a/b;end, end

%% Figures
function make_gate_figure(SYN,CFG)
S=SYN.summary; sc=CFG.syntheticScenarios; me={'DirectFlat','CalibratedFlat','FullHCL'}; Y=nan(numel(sc),3); E=Y;
for i=1:numel(sc),for j=1:3, q=strcmp(S.Scenario,sc{i})&strcmp(S.Method,me{j}); if any(q),Y(i,j)=S.GranMacroF1_Mean(q);E(i,j)=S.GranMacroF1_SD(q);end,end,end
fig=figure('Color','w','Position',[100 100 900 520]); bar(Y,'grouped'); hold on; grouped_errors(Y,E); ylim([0 1.05]);
set(gca,'XTick',1:numel(sc),'XTickLabel',sc,'FontName','Times New Roman','FontSize',12); ylabel('Granularity Macro-F1'); xlabel('Corruption mixture'); legend({'Direct Flat','Calibrated Flat','Full HCL'},'Location','northoutside','Orientation','horizontal','Box','off'); grid on; box on; export_both(fig,CFG,'gate_calibrated_flat_comparison');
end

function make_pavia_confusion_figure(PAV,CFG)
kH=matlab.lang.makeValidName('Mixed_FullHCL'); kC=matlab.lang.makeValidName('Mixed_CalibratedFlat'); if ~isfield(PAV.confusion,kH)||~isfield(PAV.confusion,kC),return;end
NH=PAV.confusion.(kH)./max(sum(PAV.confusion.(kH),2),1); NC=PAV.confusion.(kC)./max(sum(PAV.confusion.(kC),2),1); lab={'Clean','Entry','Block','Slice'};
fig=figure('Color','w','Position',[100 100 1000 430]); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(NC,[0 1]); axis image; matrix_axes(lab); title('(a) Calibrated Flat'); colorbar; add_matrix_text(NC);
nexttile; imagesc(NH,[0 1]); axis image; matrix_axes(lab); title('(b) Full HCL'); colorbar; add_matrix_text(NH); export_both(fig,CFG,'pavia_mixed_confusion_comparison');
end

function make_overlap_figure(OV,CFG)
S=OV.summary(ismember(OV.summary.ActiveMechanisms,[1 2 3]),:); [~,o]=sort(S.ActiveMechanisms);S=S(o,:);Y=[S.HCL_AAM_Mean,S.Uniform_AAM_Mean];E=[S.HCL_AAM_SD,S.Uniform_AAM_SD];
fig=figure('Color','w','Position',[100 100 800 500]);bar(Y,'grouped');hold on;grouped_errors(Y,E);ylim([0 1.05]);set(gca,'XTick',1:3,'XTickLabel',{'1','2','3'},'FontName','Times New Roman','FontSize',12);xlabel('Number of active corruption mechanisms');ylabel('Active-attribution mass');legend({'Full HCL','Uniform allocation'},'Location','northoutside','Orientation','horizontal','Box','off');grid on;box on;export_both(fig,CFG,'overlap_aam_stratified');
end

function grouped_errors(Y,E)
ng=size(Y,1); nb=size(Y,2); gw=min(.8,nb/(nb+1.5));
for i=1:nb, x=(1:ng)-gw/2+(2*i-1)*gw/(2*nb); errorbar(x,Y(:,i),E(:,i),'k.','LineWidth',1); end
end
function matrix_axes(lab),set(gca,'XTick',1:4,'XTickLabel',lab,'YTick',1:4,'YTickLabel',lab,'FontName','Times New Roman','FontSize',11);xlabel('Predicted');ylabel('True');end
function add_matrix_text(M),for r=1:size(M,1),for c=1:size(M,2),text(c,r,sprintf('%.3f',M(r,c)),'HorizontalAlignment','center','FontName','Times New Roman','FontSize',10);end,end,end
function export_both(fig,CFG,name)
try,exportgraphics(fig,fullfile(CFG.outDir,[name '.pdf']),'ContentType','vector');exportgraphics(fig,fullfile(CFG.outDir,[name '.png']),'Resolution',600);catch,saveas(fig,fullfile(CFG.outDir,[name '.png']));end;close(fig);
end

%% Aliases and recursive extraction
function a=alias_ae(),a={'ae','a_e','ae_final','final_ae','entryevidence','evidenceentry','entry_evidence'};end
function a=alias_ab(),a={'ab','a_b','ab_final','final_ab','blockevidence','evidenceblock','block_evidence'};end
function a=alias_as(),a={'as','a_s','as_final','final_as','sliceevidence','evidenceslice','slice_evidence'};end
function a=alias_po(),a={'po','p_o','po_final','final_po','pofinal','finalpo','abnormalityscore','finalabnormalityscore'};end
function a=alias_qe(),a={'qe','q_e','qe_final','final_qe','entryshare','qentry','pi_e','pie'};end
function a=alias_qb(),a={'qb','q_b','qb_final','final_qb','blockshare','qblock','pi_b','pib'};end
function a=alias_qs(),a={'qs','q_s','qs_final','final_qs','sliceshare','qslice','pi_s','pis'};end
function a=alias_gtlabel(),a={'gtlabel','gt_label','truelabel','true_label','referencelabel','reference_label','granularitygt','granularity_gt','truegranularity','true_granularity','injectedgranularity','injected_granularity'};end
function a=alias_gtentry(),a={'gtentry','gt_entry','maskentry','entrymask','trueentrymask','entrysupport','maske','supporte','ee'};end
function a=alias_gtblock(),a={'gtblock','gt_block','maskblock','blockmask','trueblockmask','blocksupport','maskb','supportb','eb'};end
function a=alias_gtslice(),a={'gtslice','gt_slice','maskslice','slicemask','trueslicemask','slicesupport','masks','supports','es'};end

function value=recursive_find(S,aliases)
value=[]; if ~isstruct(S),return;end; target=cellfun(@normalize_one,aliases,'UniformOutput',false); fn=fieldnames(S);
for i=1:numel(fn)
 if any(strcmp(normalize_one(fn{i}),target)),v=S.(fn{i});if isnumeric(v)||islogical(v)||ischar(v)||isstring(v),value=v;return;end,end
end
for i=1:numel(fn),v=S.(fn{i});if isstruct(v)&&numel(v)==1,value=recursive_find(v,aliases);if ~isempty(value),return;end,end,end
end
function x=normalize_one(x),x=regexprep(lower(char(x)),'[^a-z0-9]','');end

%% Dataset/scenario/seed and label normalization
function seed=parse_seed(file,S,CFG)
seed=NaN; v=recursive_find(S,{'seed','randomseed','rngseed'}); if isnumeric(v)&&isscalar(v),seed=double(v);end
if isfinite(seed),return;end
tok=regexp(char(file),'20\d{6}','match'); candidates=str2double(tok);
valid=[CFG.devSeeds,CFG.confirmSeeds,CFG.overlapSeeds,CFG.paviaSeeds]; hit=candidates(ismember(candidates,valid)); if ~isempty(hit),seed=hit(end);elseif ~isempty(candidates),seed=candidates(end);end
end

function dataset=classify_dataset(seed,file,CFG)
dataset="OTHER"; if ismember(seed,CFG.devSeeds),dataset="DEV";elseif ismember(seed,CFG.confirmSeeds),dataset="CONFIRM";elseif ismember(seed,CFG.overlapSeeds),dataset="OVERLAP";elseif ismember(seed,CFG.paviaSeeds),dataset="PAVIA";else
 f=upper(string(file));if contains(f,"OVERLAP"),dataset="OVERLAP";elseif contains(f,"PAVIA"),dataset="PAVIA";elseif contains(f,"CONFIRM"),dataset="CONFIRM";elseif contains(f,"DEVELOP"),dataset="DEV";end
end
end

function scenario=classify_scenario(file,S,dataset)
scenario="";v=recursive_find(S,{'scenario','case','casename','setting'});if ischar(v)||isstring(v),txt=upper(string(v));else,txt=upper(string(file));end
switch char(dataset)
 case {'DEV','CONFIRM'},if contains(txt,"EBS"),scenario="EBS";elseif contains(txt,"ES"),scenario="ES";elseif contains(txt,"BS"),scenario="BS";elseif contains(txt,"EB"),scenario="EB";end
 case 'OVERLAP',scenario="OverlapEBS";
 case 'PAVIA',if contains(txt,"MIXED"),scenario="Mixed";elseif contains(txt,"ENTRY"),scenario="Entry";elseif contains(txt,"BLOCK"),scenario="Block";elseif contains(txt,"SLICE"),scenario="Slice";end
end
end

function L=normalize_label_map(L)
L=double(L);u=unique(L(:));if all(ismember(u,[0 1 2 3])),L=uint8(L);else,error('GT label map must use 0/1/2/3.');end
end
function L=label_from_masks(e,b,s)
% Exclusive labels only. For overlaps, retain the finest-first deterministic
% convention solely for classification; overlap AAM always uses all masks.
L=zeros(size(e),'uint8');L(logical(s))=3;L(logical(b))=2;L(logical(e))=1;
end
function x=force_double(x),if isempty(x),return;elseif isnumeric(x)||islogical(x),x=double(x);else,x=[];end,end
function x=logical_or_empty(x),if isempty(x),return;elseif isnumeric(x)||islogical(x),x=logical(x);else,x=[];end,end
function raw = evaluate_synthetic_strict_cases(CASES,CAL,CFG)
% Synthetic reviewer-control evaluator using FINAL HCL state only.

rows = {};

methods = { ...
    'DirectFlat', ...
    'CalibratedFlat', ...
    'FullHCL'};

tol = 1e-12;


for i = 1:numel(CASES)

    C = CASES(i);

    % ---------------------------------------------------------
    % Strictly reload the original HCL state.
    % Never use recursively extracted attribution fields here.
    % ---------------------------------------------------------

    X = load_synthetic_final_state(C.file);

    gt = X.gtLabel;

    detMask = X.detMask;


    % ---------------------------------------------------------
    % Mandatory Full-HCL fidelity check
    % ---------------------------------------------------------

    f1OfficialCheck = ...
        macro_f1_present_strict( ...
        X.predHCL, ...
        gt);

    assert( ...
        abs(f1OfficialCheck-X.officialGranF1) < tol, ...
        ['FULL HCL GRANULARITY CONSISTENCY FAILURE\n' ...
         'File: %s\n' ...
         'Strict = %.16g\n' ...
         'Official = %.16g'], ...
        char(X.file), ...
        f1OfficialCheck, ...
        X.officialGranF1);


    for m = 1:numel(methods)

        switch methods{m}

            % =================================================
            % Direct Flat
            % =================================================

            case 'DirectFlat'

                den = ...
                    X.ae + ...
                    X.ab + ...
                    X.as + eps;

                qe = X.ae ./ den;
                qb = X.ab ./ den;
                qs = X.as ./ den;

                pred = ...
                    attribution_label_from_support( ...
                    detMask, ...
                    qe,qb,qs);


            % =================================================
            % Development-calibrated flat
            % =================================================

            case 'CalibratedFlat'

                [qe,qb,qs] = ...
                    apply_calibrated_flat( ...
                    X.ae, ...
                    X.ab, ...
                    X.as, ...
                    CAL);

                pred = ...
                    attribution_label_from_support( ...
                    detMask, ...
                    qe,qb,qs);


            % =================================================
            % Full HCL
            % =================================================

            case 'FullHCL'

                pred = X.predHCL;

        end


        % -----------------------------------------------------
        % Identical metric definition
        % -----------------------------------------------------

        M = classification_metrics( ...
            gt, ...
            pred);


        % -----------------------------------------------------
        % HARD Full-HCL assertion
        % -----------------------------------------------------

        if strcmp(methods{m},'FullHCL')

            assert( ...
                abs(M.DetectionF1-X.officialDetF1) < tol, ...
                ['FULL HCL DETECTION CONSISTENCY FAILURE\n' ...
                 'File: %s\n' ...
                 'Reviewer = %.16g\n' ...
                 'Official = %.16g'], ...
                char(X.file), ...
                M.DetectionF1, ...
                X.officialDetF1);

            assert( ...
                abs(M.GranularityMacroF1- ...
                    X.officialGranF1) < tol, ...
                ['FULL HCL GRANULARITY CONSISTENCY FAILURE\n' ...
                 'File: %s\n' ...
                 'Reviewer = %.16g\n' ...
                 'Official = %.16g'], ...
                char(X.file), ...
                M.GranularityMacroF1, ...
                X.officialGranF1);

        end


        rows(end+1,:) = { ...
            C.scenario, ...
            C.seed, ...
            methods{m}, ...
            M.DetectionF1, ...
            M.GranularityMacroF1, ...
            M.CleanFPR, ...
            M.PrecEntry, ...
            M.RecEntry, ...
            M.F1Entry, ...
            M.PrecBlock, ...
            M.RecBlock, ...
            M.F1Block, ...
            M.PrecSlice, ...
            M.RecSlice, ...
            M.F1Slice}; %#ok<AGROW>

    end

end


raw = cell2table( ...
    rows, ...
    'VariableNames',{ ...
    'Scenario', ...
    'Seed', ...
    'Method', ...
    'DetectionF1', ...
    'GranularityMacroF1', ...
    'CleanFPR', ...
    'PrecisionEntry', ...
    'RecallEntry', ...
    'F1Entry', ...
    'PrecisionBlock', ...
    'RecallBlock', ...
    'F1Block', ...
    'PrecisionSlice', ...
    'RecallSlice', ...
    'F1Slice'});


fprintf('\n');
fprintf('============================================================\n');
fprintf('STRICT SYNTHETIC CONSISTENCY: PASS\n');
fprintf('All Full-HCL cases reproduce official saved metrics.\n');
fprintf('Number of confirmation states checked = %d\n', ...
    numel(CASES));
fprintf('============================================================\n');

end
function pred = attribution_label_from_support( ...
    detMask,qe,qb,qs)

A = cat(4,qe,qb,qs);

[~,lab] = max(A,[],4);

pred = zeros(size(detMask),'uint8');

pred(detMask) = ...
    uint8(lab(detMask));

end