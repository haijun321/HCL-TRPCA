function R = make_HCL_DSP_final_figures_v3(opts)
%MAKE_HCL_DSP_FINAL_FIGURES_V3  Create the two final synthetic figures for DSP.
%
% Outputs (under <projectRoot>\figures):
%   gate_ablation_final.pdf          vector PDF for Overleaf
%   gate_ablation_final.png          300 dpi
%   gate_ablation_final.tif          300 dpi
%   gate_ablation_final.fig          editable MATLAB figure
%
%   rank_stability_final.pdf         vector PDF for Overleaf
%   rank_stability_final.png         300 dpi
%   rank_stability_final.tif         300 dpi
%   rank_stability_final.fig         editable MATLAB figure
%
% It also writes the exact plot data used:
%   gate_ablation_final_plotdata.csv
%   rank_final_plotdata.csv
%   stability_final_plotdata.csv
%
% Usage:
%   R = make_HCL_DSP_final_figures_v3();
%
% Optional explicit source paths:
%   opts.gateCsv      = '...\gate_raw.csv';
%   opts.rankCsv      = '...\rank_raw.csv';
%   opts.stabilityCsv = '...\stability_raw.csv';
%   R = make_HCL_DSP_final_figures_v3(opts);
%
% Notes:
% 1) Gate and stability figures MUST come from saved final CSV results.
%    If automatic discovery fails, pass their paths explicitly.
% 2) For the rank panel only, the validated final summary values are used as
%    a fallback if the rank CSV cannot be located.
% 3) No algorithm is rerun here; this function only reads completed results.
% 4) V3 applies publication-oriented layout changes: compact Gate legend,
%    no redundant Gate title, taller rank/stability canvas, smaller panel
%    titles, no oversized true-rank label, and unclipped perturbation axes.

if nargin < 1 || isempty(opts), opts = struct(); end
if ~isfield(opts,'gateCsv'), opts.gateCsv = ''; end
if ~isfield(opts,'rankCsv'), opts.rankCsv = ''; end
if ~isfield(opts,'stabilityCsv'), opts.stabilityCsv = ''; end
if ~isfield(opts,'closeFigures'), opts.closeFigures = false; end

% Publication-style controls
if ~isfield(opts,'fontName'),           opts.fontName = 'Times New Roman'; end
if ~isfield(opts,'axisFontSize'),       opts.axisFontSize = 8.0; end
if ~isfield(opts,'panelTitleFontSize'), opts.panelTitleFontSize = 8.8; end
if ~isfield(opts,'legendFontSize'),     opts.legendFontSize = 7.2; end
if ~isfield(opts,'labelFontSize'),      opts.labelFontSize = 8.5; end
if ~isfield(opts,'lineWidth'),          opts.lineWidth = 1.0; end

% -------------------------------------------------------------------------
% Project paths
% -------------------------------------------------------------------------
cfgPath = which('HCL_DSP_FINAL_config');
assert(~isempty(cfgPath), ['HCL_DSP_FINAL_config.m is not on the MATLAB path. ' ...
    'Add the final package root first.']);
projectRoot = fileparts(cfgPath);
C = HCL_DSP_FINAL_config();
outputRoot = C.outputRoot;
figDir = fullfile(projectRoot,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n============================================================\n');
fprintf(' HCL-DSP FINAL PUBLICATION FIGURES\n');
fprintf('============================================================\n');
fprintf('Project root : %s\n',projectRoot);
fprintf('Results root : %s\n',outputRoot);
fprintf('Figure dir   : %s\n\n',figDir);

% -------------------------------------------------------------------------
% Locate source CSV files
% -------------------------------------------------------------------------
gateCsv = resolve_csv(opts.gateCsv, outputRoot, 'gate', @is_gate_table);
rankCsv = resolve_csv(opts.rankCsv, outputRoot, 'rank', @is_rank_table, true);
stabCsv = resolve_csv(opts.stabilityCsv, outputRoot, 'stability', @is_stability_table);

fprintf('Gate source      : %s\n',gateCsv);
if isempty(rankCsv)
    fprintf('Rank source      : <not found; validated frozen summary fallback>\n');
else
    fprintf('Rank source      : %s\n',rankCsv);
end
fprintf('Stability source : %s\n\n',stabCsv);

% -------------------------------------------------------------------------
% Read and summarize final data
% -------------------------------------------------------------------------
G = build_gate_plotdata(readtable(gateCsv));
if isempty(rankCsv)
    K = rank_fallback_final();
else
    K = build_rank_plotdata(readtable(rankCsv));
end
S = build_stability_plotdata(readtable(stabCsv));

% Save exact plot data for audit/reproducibility
writetable(G,fullfile(figDir,'gate_ablation_final_plotdata.csv'));
writetable(K,fullfile(figDir,'rank_final_plotdata.csv'));
writetable(S,fullfile(figDir,'stability_final_plotdata.csv'));

% -------------------------------------------------------------------------
% Figure 1: Gate ablation
% Publication design:
%   - no redundant in-figure title
%   - compact one-line legend
%   - restrained fonts and bar width
%   - vector PDF + 300 dpi raster + editable FIG
% -------------------------------------------------------------------------
f1 = figure('Color','w','Units','inches','Position',[1 1 7.15 2.95], ...
    'Renderer','painters','Name','Gate ablation final');

tl1 = tiledlayout(f1,1,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl1,1);
hold(ax,'on');

scenarioOrder = ["EB","ES","BS","EBS"];
variantOrder = ["EntryOnly","Flat","NoBlock","NoSlice","FullHCL"];
variantLabel = {'Entry only','Flat','No block','No slice','Full HCL'};

mu = nan(numel(scenarioOrder),numel(variantOrder));
sd = nan(size(mu));
for i=1:numel(scenarioOrder)
    for j=1:numel(variantOrder)
        idx = G.Scenario==scenarioOrder(i) & G.Variant==variantOrder(j);
        assert(sum(idx)==1,'Missing/duplicate gate cell: %s / %s.', ...
            scenarioOrder(i),variantOrder(j));
        mu(i,j) = G.MeanLevelMacroF1(idx);
        sd(i,j) = G.SDLevelMacroF1(idx);
    end
end

b = bar(ax,1:numel(scenarioOrder),mu,'grouped','BarWidth',0.76);
for j=1:numel(b)
    b(j).EdgeColor = 'none';
end

for j=1:numel(b)
    x = b(j).XEndPoints;
    errorbar(ax,x,mu(:,j),sd(:,j),'k','LineStyle','none', ...
        'LineWidth',0.65,'CapSize',2.5);
end

xlim(ax,[0.50 4.50]);
ylim(ax,[0 1.035]);
xticks(ax,1:4);
xticklabels(ax,cellstr(scenarioOrder));

xlabel(ax,'Corruption mixture', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);
ylabel(ax,'Level Macro-F1', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);

grid(ax,'on');
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.GridAlpha = 0.13;

% Compact publication legend above the axes.
lgd1 = legend(ax,b,variantLabel, ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'NumColumns',5);
lgd1.Layout.Tile = 'north';
lgd1.FontName = opts.fontName;
lgd1.FontSize = opts.legendFontSize;

% No plot title: the paper caption carries the figure title and interpretation.
style_axis(ax,opts);

export_publication_set(f1,figDir,'gate_ablation_final');

% -------------------------------------------------------------------------
% Figure 2: Rank + stability, three panels
% Publication design:
%   - taller canvas to prevent bottom clipping
%   - smaller panel titles
%   - no oversized "True rank" text
%   - concise perturbation labels
%   - zooming of panel (c) disclosed in the manuscript caption, not title
% -------------------------------------------------------------------------
f2 = figure('Color','w','Units','inches','Position',[1 1 7.15 3.45], ...
    'Renderer','painters','Name','Rank stability final');

tl = tiledlayout(f2,1,3,'TileSpacing','compact','Padding','loose');

% ---------------------------- (a) Rank ----------------------------
ax1 = nexttile(tl,1);
hold(ax1,'on');

errorbar(ax1,K.Rank,K.MeanNRE,K.SDNRE,'-o', ...
    'LineWidth',1.05,'MarkerSize',4.2,'CapSize',3.5);

% True rank is indicated only by a light dashed reference line.
xline(ax1,5,'--','LineWidth',0.85);

xticks(ax1,K.Rank);
xlim(ax1,[min(K.Rank)-0.35 max(K.Rank)+0.35]);

ymax = max(K.MeanNRE + K.SDNRE);
ylim(ax1,[0 max(0.56,1.07*ymax)]);

xlabel(ax1,'Fitted rank', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);
ylabel(ax1,'NRE', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);

title(ax1,'(a) Rank misspecification', ...
    'FontName',opts.fontName, ...
    'FontSize',opts.panelTitleFontSize, ...
    'FontWeight','normal');

grid(ax1,'on');
ax1.GridAlpha = 0.13;
style_axis(ax1,opts);

% ------------------------- (b) Stability bound -------------------------
ax2 = nexttile(tl,2);
hold(ax2,'on');

xp = 1:height(S);

if all(isfinite(S.SDObservedScorePerturbation))
    hObs = errorbar(ax2,xp,S.MeanObservedScorePerturbation, ...
        S.SDObservedScorePerturbation,'-o', ...
        'LineWidth',1.00,'MarkerSize',4.0,'CapSize',3);
else
    hObs = plot(ax2,xp,S.MeanObservedScorePerturbation,'-o', ...
        'LineWidth',1.00,'MarkerSize',4.0);
end

hBnd = plot(ax2,xp,S.TheoreticalBound,'--s', ...
    'LineWidth',1.00,'MarkerSize',3.8);

xticks(ax2,xp);
xticklabels(ax2,delta_labels(S.Delta));
xtickangle(ax2,25);
xlim(ax2,[0.65 height(S)+0.35]);

ymax2 = max([S.MeanObservedScorePerturbation(:); S.TheoreticalBound(:)]);
ylim(ax2,[0 1.08*ymax2]);

xlabel(ax2,'Perturbation, \delta', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);
ylabel(ax2,'Score perturbation', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);

title(ax2,'(b) Stability bound', ...
    'FontName',opts.fontName, ...
    'FontSize',opts.panelTitleFontSize, ...
    'FontWeight','normal');

lgd2 = legend(ax2,[hObs hBnd],{'Observed','B(\delta)'}, ...
    'Location','northwest','Box','off');
lgd2.FontName = opts.fontName;
lgd2.FontSize = opts.legendFontSize;

grid(ax2,'on');
ax2.GridAlpha = 0.13;
style_axis(ax2,opts);

% ------------------------- (c) Mask agreement -------------------------
ax3 = nexttile(tl,3);
hold(ax3,'on');

fixedMean = 100*S.MeanFixedNMaskAgreement;
fullMean  = 100*S.MeanFullRuleMaskAgreement;
fixedSD   = 100*S.SDFixedNMaskAgreement;
fullSD    = 100*S.SDFullRuleMaskAgreement;

if all(isfinite(fixedSD))
    hF = errorbar(ax3,xp,fixedMean,fixedSD,'-o', ...
        'LineWidth',1.00,'MarkerSize',4.0,'CapSize',3);
else
    hF = plot(ax3,xp,fixedMean,'-o', ...
        'LineWidth',1.00,'MarkerSize',4.0);
end

if all(isfinite(fullSD))
    hR = errorbar(ax3,xp,fullMean,fullSD,'--s', ...
        'LineWidth',1.00,'MarkerSize',3.8,'CapSize',3);
else
    hR = plot(ax3,xp,fullMean,'--s', ...
        'LineWidth',1.00,'MarkerSize',3.8);
end

xticks(ax3,xp);
xticklabels(ax3,delta_labels(S.Delta));
xtickangle(ax3,25);
xlim(ax3,[0.65 height(S)+0.35]);

% Expanded y-axis for small differences.  The manuscript caption should
% explicitly state that this axis is expanded.
vals = [fixedMean(:)-fixedSD(:); fullMean(:)-fullSD(:)];
vals = vals(isfinite(vals));
if isempty(vals)
    low = 99.93;
else
    low = floor((min(vals)-0.002)*100)/100;
    low = max(99.90,min(low,99.95));
end
ylim(ax3,[low 100.002]);

xlabel(ax3,'Perturbation, \delta', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);
ylabel(ax3,'Mask agreement (%)', ...
    'FontName',opts.fontName,'FontSize',opts.labelFontSize);

title(ax3,'(c) Recovery-mask agreement', ...
    'FontName',opts.fontName, ...
    'FontSize',opts.panelTitleFontSize, ...
    'FontWeight','normal');

lgd3 = legend(ax3,[hF hR],{'Fixed-n','Full rule'}, ...
    'Location','southwest','Box','off');
lgd3.FontName = opts.fontName;
lgd3.FontSize = opts.legendFontSize;

grid(ax3,'on');
ax3.GridAlpha = 0.13;
style_axis(ax3,opts);

export_publication_set(f2,figDir,'rank_stability_final');

% Return paths and plot data
R = struct();
R.projectRoot = projectRoot;
R.outputRoot = outputRoot;
R.figureDir = figDir;
R.sources = struct('gateCsv',gateCsv,'rankCsv',rankCsv,'stabilityCsv',stabCsv);
R.gate = G;
R.rank = K;
R.stability = S;
R.files = struct( ...
    'gatePDF',fullfile(figDir,'gate_ablation_final.pdf'), ...
    'gatePNG',fullfile(figDir,'gate_ablation_final.png'), ...
    'gateTIF',fullfile(figDir,'gate_ablation_final.tif'), ...
    'gateFIG',fullfile(figDir,'gate_ablation_final.fig'), ...
    'rankStabilityPDF',fullfile(figDir,'rank_stability_final.pdf'), ...
    'rankStabilityPNG',fullfile(figDir,'rank_stability_final.png'), ...
    'rankStabilityTIF',fullfile(figDir,'rank_stability_final.tif'), ...
    'rankStabilityFIG',fullfile(figDir,'rank_stability_final.fig'));

fprintf('\nGenerated publication files:\n');
fprintf('  %s\n',R.files.gatePDF);
fprintf('  %s\n',R.files.gatePNG);
fprintf('  %s\n',R.files.gateTIF);
fprintf('  %s\n',R.files.gateFIG);
fprintf('  %s\n',R.files.rankStabilityPDF);
fprintf('  %s\n',R.files.rankStabilityPNG);
fprintf('  %s\n',R.files.rankStabilityTIF);
fprintf('  %s\n',R.files.rankStabilityFIG);

if opts.closeFigures
    close(f1); close(f2);
end
end

% =========================================================================
% Data discovery
% =========================================================================
function pathOut = resolve_csv(explicitPath,root,keyword,validator,allowEmpty)
if nargin < 5, allowEmpty = false; end
pathOut = '';
if ~isempty(explicitPath)
    assert(exist(explicitPath,'file')==2,'CSV not found: %s',explicitPath);
    T = readtable(explicitPath);
    assert(validator(T),'Explicit %s CSV does not contain expected columns.',keyword);
    pathOut = explicitPath;
    return;
end

D = dir(fullfile(root,'**','*.csv'));
if isempty(D)
    if allowEmpty, return; end
    error('No CSV files found under %s.',root);
end

names = lower(string({D.name}));
hit = contains(names,lower(keyword));
isRaw = contains(names,'raw');
isSummary = contains(names,'summary');
% Prefer final RAW data, then summary files, then any other matching CSV.
ord = [find(hit & isRaw) find(hit & ~isRaw & isSummary) ...
       find(hit & ~isRaw & ~isSummary) find(~hit & isRaw) find(~hit & ~isRaw)];
ord = unique(ord,'stable');

for ii = ord
    p = fullfile(D(ii).folder,D(ii).name);
    try
        T = readtable(p,'VariableNamingRule','preserve');
        if validator(T)
            pathOut = p;
            return;
        end
    catch
    end
end

if allowEmpty, return; end
error(['Could not automatically locate the final %s CSV under:\n%s\n' ...
    'Pass opts.%sCsv explicitly.'],keyword,root,keyword);
end

function tf = is_gate_table(T)
v = lower(string(T.Properties.VariableNames));
tf = any(contains(v,'scenario')) && ...
     (any(contains(v,'variant')) || any(contains(v,'method'))) && ...
     any(contains(v,'macro') & contains(v,'f1'));
end

function tf = is_rank_table(T)
v = lower(string(T.Properties.VariableNames));
tf = any(contains(v,'rank')) && any(strcmp(v,'nre') | contains(v,'meannre'));
end

function tf = is_stability_table(T)
% Accept both the original generic schema and the FINAL stability schema:
% RequestedDelta, ActualEvidenceDelta, ActualPoInf, TheoreticalBound,
% FixedNMaskAgreement, FullRuleMaskAgreement.
v = lower(regexprep(string(T.Properties.VariableNames),'[^a-zA-Z0-9]',''));
hasDelta = any(v=="delta" | v=="requesteddelta" | ...
               contains(v,'evidencedelta') | contains(v,'perturb'));
hasAgree = any(contains(v,'fixednmaskagreement') | ...
               contains(v,'fullrulemaskagreement') | ...
               (contains(v,'agreement') & (contains(v,'fixed') | contains(v,'full'))));
hasScore = any(v=="actualpoinf" | v=="meanactualpoinf" | ...
               contains(v,'scoreperturb') | contains(v,'scoredelta') | ...
               contains(v,'actualscore'));
tf = hasDelta && hasAgree && hasScore;
end

% =========================================================================
% Gate table
% =========================================================================
function G = build_gate_plotdata(T)
scenarioVar = pick_var(T,["Scenario"],["scenario"]);
variantVar = pick_var(T,["Variant","Method"],["variant"]);
if isempty(variantVar)
    variantVar = pick_var(T,["Variant","Method"],["method"]);
end
metricVar = pick_var(T,["LevelMacroF1"],["macro","f1"]);

meanVar = pick_var(T,["MeanLevelMacroF1","MeanMacroF1"],["mean","macro","f1"]);
sdVar = pick_var(T,["SDLevelMacroF1","StdLevelMacroF1","SDMacroF1","StdMacroF1"],["macro","f1","sd"]);
if isempty(sdVar), sdVar = pick_var(T,[],["macro","f1","std"]); end

scenario = upper(string(T.(scenarioVar)));
variant = strings(height(T),1);
for i=1:height(T), variant(i)=canon_variant(T.(variantVar)(i)); end

scenarioOrder = ["EB","ES","BS","EBS"];
variantOrder = ["EntryOnly","Flat","NoBlock","NoSlice","FullHCL"];
rows = cell(numel(scenarioOrder)*numel(variantOrder),6);
r=0;

for i=1:numel(scenarioOrder)
    for j=1:numel(variantOrder)
        idx = scenario==scenarioOrder(i) & variant==variantOrder(j);
        assert(any(idx),'Gate source missing %s / %s.',scenarioOrder(i),variantOrder(j));
        r=r+1;
        if ~isempty(meanVar)
            m = mean(double(T.(meanVar)(idx)),'omitnan');
            if ~isempty(sdVar)
                s = mean(double(T.(sdVar)(idx)),'omitnan');
            else
                s = NaN;
            end
            n = sum(idx);
        else
            x = double(T.(metricVar)(idx));
            m = mean(x,'omitnan');
            s = std(x,'omitnan');
            n = sum(isfinite(x));
        end
        rows(r,:) = {scenarioOrder(i),variantOrder(j),m,s,n,true};
    end
end

G = cell2table(rows,'VariableNames', ...
    {'Scenario','Variant','MeanLevelMacroF1','SDLevelMacroF1','N','FromFinalCSV'});
G.Scenario = string(G.Scenario);
G.Variant = string(G.Variant);
end

function s = canon_variant(x)
t = lower(regexprep(char(string(x)),'[^a-zA-Z0-9]',''));
if contains(t,'noblock')
    s="NoBlock";
elseif contains(t,'noslice')
    s="NoSlice";
elseif contains(t,'flat')
    s="Flat";
elseif contains(t,'entry') && (contains(t,'only') || contains(t,'entryonly'))
    s="EntryOnly";
elseif contains(t,'full') || contains(t,'fullhcl') || strcmp(t,'hcl')
    s="FullHCL";
else
    error('Unrecognized gate variant label: %s',string(x));
end
end

% =========================================================================
% Rank table
% =========================================================================
function K = build_rank_plotdata(T)
rankVar = pick_var(T,["FittedRank","Rank","RecoveryRank"],["rank"]);
meanVar = pick_var(T,["MeanNRE"],["mean","nre"]);
sdVar = pick_var(T,["SDNRE","StdNRE"],["nre","sd"]);
if isempty(sdVar), sdVar = pick_var(T,[],["nre","std"]); end
nreVar = pick_var(T,["NRE"],["nre"]);

% If a method column exists, retain HCL rows only when possible
methodVar = pick_var(T,["Method"],["method"]);
idxMethod = true(height(T),1);
if ~isempty(methodVar)
    ms = lower(string(T.(methodVar)));
    h = contains(ms,'hcl');
    if any(h), idxMethod = h; end
end

rv = double(T.(rankVar));
ranks = unique(rv(idxMethod & isfinite(rv)))';
wanted = [3 4 5 6 8];
if all(ismember(wanted,ranks)), ranks=wanted; end

rows = cell(numel(ranks),5);
for i=1:numel(ranks)
    idx = idxMethod & rv==ranks(i);
    if ~isempty(meanVar)
        m=mean(double(T.(meanVar)(idx)),'omitnan');
        if ~isempty(sdVar), s=mean(double(T.(sdVar)(idx)),'omitnan'); else, s=NaN; end
        n=sum(idx);
    else
        x=double(T.(nreVar)(idx));
        m=mean(x,'omitnan'); s=std(x,'omitnan'); n=sum(isfinite(x));
    end
    rows(i,:)={ranks(i),m,s,n,false};
end
K=cell2table(rows,'VariableNames',{'Rank','MeanNRE','SDNRE','N','Fallback'});
end

function K = rank_fallback_final()
K = table([3;4;5;6;8], ...
    [0.50786;0.33203;0.008736;0.08706;0.29788], ...
    [0.00514;0.00623;0.000035;0.04200;0.14250], ...
    [10;10;10;10;10], ...
    true(5,1), ...
    'VariableNames',{'Rank','MeanNRE','SDNRE','N','Fallback'});
warning(['Rank CSV was not found. Panel (a) uses the already validated final ' ...
    'summary values recorded in the final experiment report.']);
end

% =========================================================================
% Stability table
% =========================================================================
function S = build_stability_plotdata(T)
% Supports FINAL raw:
% Seed, RequestedDelta, ActualEvidenceDelta, ActualPoInf, TheoreticalBound,
% BoundSatisfied, BoundaryGap, TopNConditionSatisfied,
% FixedNMaskAgreement, FullRuleMaskAgreement, ...
%
% Supports FINAL summary:
% RequestedDelta, GroupCount, mean_ActualEvidenceDelta,
% mean_ActualPoInf, std_ActualPoInf, mean_TheoreticalBound,
% mean_FixedNMaskAgreement, std_FixedNMaskAgreement,
% mean_FullRuleMaskAgreement, std_FullRuleMaskAgreement, ...

deltaVar = pick_var(T, ...
    ["RequestedDelta","Delta","EvidenceDelta","Perturbation"],["delta"]);
if isempty(deltaVar), deltaVar=pick_var(T,[],["perturb"]); end

boundVar = pick_var(T, ...
    ["TheoreticalBound","mean_TheoreticalBound"],["theoretical","bound"]);

% RAW observed score perturbation.  In the final stability runner this is
% ||p_tilde-p||_inf and is stored as ActualPoInf.
obsVar = pick_var(T, ...
    ["ActualPoInf","ObservedScorePerturbation","ActualScorePerturbation", ...
     "ScorePerturbation","ScoreDeltaInf","ObservedDeltaP"],["actual","po","inf"]);
if isempty(obsVar), obsVar=pick_var(T,[],["score","perturb"]); end
if isempty(obsVar), obsVar=pick_var(T,[],["score","delta"]); end

fixedVar = pick_var(T, ...
    ["FixedNMaskAgreement","FixedNAgreement","MaskAgreementFixedN"], ...
    ["fixed","agreement"]);
fullVar = pick_var(T, ...
    ["FullRuleMaskAgreement","FullRuleAgreement","MaskAgreementFull"], ...
    ["full","agreement"]);

% SUMMARY columns
meanObsVar = pick_var(T, ...
    ["mean_ActualPoInf","MeanActualPoInf","MeanObservedScorePerturbation", ...
     "MeanActualScorePerturbation","MeanScorePerturbation"], ...
    ["mean","actual","po","inf"]);
if isempty(meanObsVar), meanObsVar=pick_var(T,[],["mean","score","perturb"]); end

sdObsVar = pick_var(T, ...
    ["std_ActualPoInf","SDActualPoInf","StdActualPoInf", ...
     "SDObservedScorePerturbation","StdObservedScorePerturbation"], ...
    ["std","actual","po","inf"]);
if isempty(sdObsVar), sdObsVar=pick_var(T,[],["sd","actual","po","inf"]); end
if isempty(sdObsVar), sdObsVar=pick_var(T,[],["score","perturb","std"]); end

meanBoundVar = pick_var(T, ...
    ["mean_TheoreticalBound","MeanTheoreticalBound"],["mean","theoretical","bound"]);
sdBoundVar = pick_var(T, ...
    ["std_TheoreticalBound","SDTheoreticalBound","StdTheoreticalBound"], ...
    ["std","theoretical","bound"]);

meanFixedVar = pick_var(T, ...
    ["mean_FixedNMaskAgreement","MeanFixedNMaskAgreement","MeanFixedNAgreement"], ...
    ["mean","fixed","agreement"]);
sdFixedVar = pick_var(T, ...
    ["std_FixedNMaskAgreement","SDFixedNMaskAgreement","StdFixedNMaskAgreement"], ...
    ["std","fixed","agreement"]);
if isempty(sdFixedVar), sdFixedVar=pick_var(T,[],["sd","fixed","agreement"]); end

meanFullVar = pick_var(T, ...
    ["mean_FullRuleMaskAgreement","MeanFullRuleMaskAgreement","MeanFullRuleAgreement"], ...
    ["mean","full","agreement"]);
sdFullVar = pick_var(T, ...
    ["std_FullRuleMaskAgreement","SDFullRuleMaskAgreement","StdFullRuleMaskAgreement"], ...
    ["std","full","agreement"]);
if isempty(sdFullVar), sdFullVar=pick_var(T,[],["sd","full","agreement"]); end

countVar = pick_var(T,["GroupCount","N"],["group","count"]);

hasRaw = ~isempty(deltaVar) && ~isempty(obsVar) && ...
         ~isempty(fixedVar) && ~isempty(fullVar);
hasSummary = ~isempty(deltaVar) && ~isempty(meanObsVar) && ...
             ~isempty(meanFixedVar) && ~isempty(meanFullVar);

assert(hasRaw || hasSummary, ...
    ['Could not identify stability columns. Found variables:\n%s'], ...
    strjoin(string(T.Properties.VariableNames),', '));

dv = double(T.(deltaVar));
deltas = sort(unique(dv(isfinite(dv))))';
rows = cell(numel(deltas),11);

for i=1:numel(deltas)
    idx = abs(dv-deltas(i)) <= max(1e-14,eps(max(1,abs(deltas(i))))*10);

    if hasSummary
        mo = mean(double(T.(meanObsVar)(idx)),'omitnan');
        mf = mean(double(T.(meanFixedVar)(idx)),'omitnan');
        mr = mean(double(T.(meanFullVar)(idx)),'omitnan');

        if ~isempty(sdObsVar)
            so = mean(double(T.(sdObsVar)(idx)),'omitnan');
        else
            so = NaN;
        end
        if ~isempty(sdFixedVar)
            sf = mean(double(T.(sdFixedVar)(idx)),'omitnan');
        else
            sf = NaN;
        end
        if ~isempty(sdFullVar)
            sr = mean(double(T.(sdFullVar)(idx)),'omitnan');
        else
            sr = NaN;
        end

        if ~isempty(meanBoundVar)
            bnd = mean(double(T.(meanBoundVar)(idx)),'omitnan');
        elseif ~isempty(boundVar)
            bnd = mean(double(T.(boundVar)(idx)),'omitnan');
        else
            bnd = min(1,2.4*deltas(i)+1.4*sqrt(max(deltas(i),0)));
        end

        if ~isempty(sdBoundVar)
            sbnd = mean(double(T.(sdBoundVar)(idx)),'omitnan');
        else
            sbnd = NaN;
        end

        if ~isempty(countVar)
            n = round(mean(double(T.(countVar)(idx)),'omitnan'));
        else
            n = sum(idx);
        end
    else
        xo = double(T.(obsVar)(idx));
        xf = double(T.(fixedVar)(idx));
        xr = double(T.(fullVar)(idx));

        mo = mean(xo,'omitnan'); so = std(xo,'omitnan');
        mf = mean(xf,'omitnan'); sf = std(xf,'omitnan');
        mr = mean(xr,'omitnan'); sr = std(xr,'omitnan');

        if ~isempty(boundVar)
            xb = double(T.(boundVar)(idx));
            bnd = mean(xb,'omitnan');
            sbnd = std(xb,'omitnan');
        else
            bnd = min(1,2.4*deltas(i)+1.4*sqrt(max(deltas(i),0)));
            sbnd = NaN;
        end
        n = sum(isfinite(xo));
    end

    rows(i,:) = {deltas(i),mo,so,bnd,sbnd,mf,sf,mr,sr,n,true};
end

S = cell2table(rows,'VariableNames', ...
    {'Delta','MeanObservedScorePerturbation','SDObservedScorePerturbation', ...
     'TheoreticalBound','SDTheoreticalBound', ...
     'MeanFixedNMaskAgreement','SDFixedNMaskAgreement', ...
     'MeanFullRuleMaskAgreement','SDFullRuleMaskAgreement', ...
     'N','FromFinalCSV'});

% Final integrity checks for the publication plot.
assert(height(S)==5, ...
    'Expected five stability perturbation levels, found %d.',height(S));
assert(all(S.MeanObservedScorePerturbation <= S.TheoreticalBound + 1e-12), ...
    'At least one mean observed perturbation exceeds the theoretical bound.');
end

% =========================================================================
% Utilities
% =========================================================================
function name = pick_var(T,preferred,tokens)
vars = string(T.Properties.VariableNames);
nv = lower(regexprep(vars,'[^a-zA-Z0-9]',''));
name = '';

for p = string(preferred)
    np = lower(regexprep(p,'[^a-zA-Z0-9]',''));
    k = find(nv==np,1);
    if ~isempty(k), name=char(vars(k)); return; end
end

if ~isempty(tokens)
    mask=true(size(vars));
    for t=string(tokens)
        nt=lower(regexprep(t,'[^a-zA-Z0-9]',''));
        mask = mask & contains(nv,nt);
    end
    k=find(mask,1);
    if ~isempty(k), name=char(vars(k)); return; end
end
end

function labels = delta_labels(d)
labels = strings(numel(d),1);
for i=1:numel(d)
    x=d(i);
    if x==0
        labels(i)="0";
    elseif abs(x-1e-4)<1e-12
        labels(i)="10^{-4}";
    elseif abs(x-1e-3)<1e-12
        labels(i)="10^{-3}";
    elseif abs(x-1e-2)<1e-12
        labels(i)="10^{-2}";
    elseif abs(x-5e-2)<1e-12
        labels(i)="5\times10^{-2}";
    else
        labels(i)=string(sprintf('%.2g',x));
    end
end
labels=cellstr(labels);
end

function style_axis(ax,opts)
set(ax, ...
    'FontName',opts.fontName, ...
    'FontSize',opts.axisFontSize, ...
    'LineWidth',0.75, ...
    'TickDir','out', ...
    'TickLength',[0.015 0.015], ...
    'Box','off', ...
    'Layer','top');
end

function export_publication_set(fig,figDir,baseName)
pdfFile = fullfile(figDir,[baseName '.pdf']);
pngFile = fullfile(figDir,[baseName '.png']);
tifFile = fullfile(figDir,[baseName '.tif']);
figFile = fullfile(figDir,[baseName '.fig']);

drawnow;
savefig(fig,figFile);
exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,pngFile,'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,tifFile,'Resolution',300,'BackgroundColor','white');

fprintf('[OK] %s -> PDF(vector), PNG(300 dpi), TIF(300 dpi), FIG\n',baseName);
end
