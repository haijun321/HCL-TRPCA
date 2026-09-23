function R = run_PaviaU_RANK_BUDGET_DIAG_V2_1(mode,resumeRoot)
%RUN_PAVIAU_RANK_BUDGET_DIAG_V2_1
% Paper-oriented rank x iteration-budget diagnosis for controlled PaviaU.
%
% V2.1 POSTPROCESSING FIX
% -------------------------
% REPRO contains only budget=80. Therefore no 160-vs-320 plateau rows exist.
% V2 attempted cell2table(rows,...) with rows={} in summarize_plateau(),
% causing a VariableNames/variable-count error AFTER all 10 expensive HCL
% runs had already completed. V2.1 returns a correctly typed empty plateau
% table when a mode does not contain both budgets 160 and 320.
%
% This fix changes NO experiment, NO HCL setting, NO seed, and NO metric.
% Existing V2 raw results can be resumed without rerunning completed HCL jobs.
%
% IMPORTANT SCIENTIFIC ROLE
% -------------------------
% This diagnostic does NOT replace the frozen primary experiment
% run_PaviaU_CONTROLLED_FINAL_V3_2.m.
%
% The primary Pavia experiment remains:
%   rank   = 25
%   budget = 80 updates
%   tol    = 1e-5
%
% The present file explains the large reconstruction gap observed in that
% frozen primary experiment by separating two factors only:
%   (1) projection rank / model order;
%   (2) maximum recovery-iteration budget.
%
% No HCL threshold, evidence calibration, recovery floor, warm-start rule,
% corruption generator, amplitude, ROI, normalization, or metric is retuned.
%
% PREVIOUSLY COMPLETED REFERENCE-BASED RANK AUDIT
% -----------------------------------------------
% The exact Pavia reference/preprocessing used here was independently audited
% BEFORE this V2 rank-budget experiment. The median sampled-Fourier-slice
% energy ranks (mode-3 FFT, stride 4, 26 sampled slices) were:
%
%   90% energy : r90 = 30
%   95% energy : r95 = 61
%   98% energy : r98 = 105
%   99% energy : r99 = 134
%
% These values are post-hoc REFERENCE-BASED diagnostic landmarks. They are
% NOT deployable rank estimates and are NOT "true tubal ranks".
%
% REQUIRED RUN ORDER
% ------------------
%   1) SMOKE
%      one Mixed seed; r=25; budgets [80 160 320]
%      -> checks whether budget alone changes recovery while the frozen mask
%         remains invariant.
%
%   2) REPRO
%      ten Mixed seeds; r=25; budget=80
%      -> must reproduce the manuscript-level frozen primary Pavia result.
%
%   3) PILOT
%      three Mixed seeds;
%      ranks [25 30 61 105 134];
%      budgets [80 160 320]
%      -> boundary diagnosis, including r99=134 as an intentionally extreme
%         high-rank pilot condition.
%
%   4) DIAG
%      ten Mixed seeds;
%      ranks [25 30 61 105];
%      budgets [80 160 320]
%      -> final quantitative rank-budget diagnosis used for the paper.
%
% Why r99=134 is PILOT-only:
%   r99 is >50% of the maximum matrix rank (256) and is used only as an
%   extreme boundary check. The final 10-seed diagnostic grid is frozen
%   BEFORE the pilot at [25 30 61 105]; this prevents post-hoc selection
%   based on the pilot NRE.
%
% MODES
% -----
% SMOKE : seed 20280001, rank 25, budgets 80/160/320.
% REPRO : seeds 20280001:20280010, rank 25, budget 80.
% PILOT : seeds 20280001:20280003, ranks 25/30/61/105/134, all budgets.
% DIAG  : seeds 20280001:20280010, ranks 25/30/61/105, all budgets.
%
% OPTIONAL RESUME
% ---------------
%   R = run_PaviaU_RANK_BUDGET_DIAG_V2_1('DIAG', existingOutputFolder)
% reloads Pavia_RBD_raw.csv (or the PARTIAL file) and skips completed
% Seed-Rank-Budget rows.
%
% CORE OUTPUTS
% ------------
% Pavia_RBD_rank_landmarks.csv
% Pavia_RBD_protocol.csv
% Pavia_RBD_raw.csv
% Pavia_RBD_summary.csv
% Pavia_RBD_publication_summary.csv
% Pavia_RBD_reproduction.csv
% Pavia_RBD_effects.csv
% Pavia_RBD_plateau.csv
% Pavia_RBD_gate_report.txt
% Pavia_RBD_results_small.mat
%
% FIGURES
% -------
% figures/Pavia_rank_budget_diagnostic.*
% figures/Pavia_rank_structural_stability.*
% figures/Pavia_rank_budget_effects.*
%
% INTERPRETATION RULE
% -------------------
% This is a mechanism/failure diagnosis, not a parameter-tuning benchmark.
% Do NOT take the best (rank,budget) observed on these same seeds and present
% it as an independent performance result. If the primary Pavia protocol is
% changed after diagnosis, use fresh confirmation seeds.
%
% MATLAB target: R2023b-compatible syntax.
% -------------------------------------------------------------------------

if nargin < 1 || isempty(mode)
    mode = 'SMOKE';
end
if nargin < 2
    resumeRoot = '';
end

mode = upper(char(string(mode)));
validModes = {'SMOKE','REPRO','PILOT','DIAG'};
assert(ismember(mode,validModes), ...
    'Mode must be SMOKE, REPRO, PILOT, or DIAG.');

%% ------------------------------------------------------------------------
% 0. Project/config/data
% -------------------------------------------------------------------------
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));

assert(exist('HCL_DSP_FINAL_config','file')==2, ...
    'HCL_DSP_FINAL_config.m is not on the MATLAB path.');
C = HCL_DSP_FINAL_config();

paviaFile = resolve_pavia_file(C,projectRoot);

%% ------------------------------------------------------------------------
% 1. Frozen primary protocol + diagnostic-only factors
% -------------------------------------------------------------------------
P = struct();
P.version            = 'PAVIA_RANK_BUDGET_DIAG_V2_1';

% Same Pavia reference/preprocessing as V3.2.
P.roiRows            = 178:433;
P.roiCols            = 43:298;
P.normalizationMax   = 8000;
P.rankStride         = 4;

% Same controlled Mixed corruption as V3.2.
P.blockSizes         = [8 12 16 24];
P.ampSigmaMin        = 6;
P.ampSigmaMax        = 10;
P.mixedScenario      = mk_scenario('Mixed',0.05,0.10,0.10);
P.mixedScenarioIndex = 4; % preserves V3.2 generatorSeed = seed + 1000*4

% Same HCL frozen scientific settings.
P.detectThreshold    = 0.50;
P.epsRecExpected     = 1e-3;
P.relTolExpected     = 1e-5;
P.primaryRank        = 25;
P.primaryBudget      = 80;
P.maxUpdatesExpected = P.primaryBudget; % used only by interface preflight

% ---------------------------------------------------------------------
% Previously audited reference-based rank landmarks. These values are
% frozen BEFORE SMOKE/PILOT/DIAG and are never selected from recovery NRE.
% ----------------------------------------------------------------------
P.rank90             = 30;
P.rank95             = 61;
P.rank98             = 105;
P.rank99             = 134;
P.rankAuditStride    = 4;
P.rankAuditSlices    = 26;
P.rankAuditTargets   = [0.90 0.95 0.98 0.99];

% Diagnostic grids.
P.rankGridPilot      = [25 30 61 105 134]; % r99 is boundary-only pilot
P.rankGridDiag       = [25 30 61 105];     % final 10-seed paper diagnosis
P.budgetGrid         = [80 160 320];

% Seeds: same diagnostic set as the current controlled Pavia experiment.
P.seedsAll           = 20280001:20280010;
P.seedsPilot         = 20280001:20280003;
P.seedSmoke          = 20280001;

% Reference identity guard. This is NOT a fitted parameter.
P.referenceSigmaExpected = 0.11118511;
P.referenceSigmaAbsTol   = 5e-8;

% Current manuscript reference values for the reproduction gate only.
% These are NOT optimization targets.
P.paperMixedNREMean      = 0.37319;
P.paperMixedNREStd       = 0.00477;
P.paperMixedDetectionF1  = 0.9766;
P.paperMixedMaskF1       = 0.9896;
P.paperMixedLevelF1      = 0.9557;

P.reproGreenAbsTol       = 0.0020;
P.reproFailAbsTol        = 0.0100;
P.reproStructuralTol     = 0.0150;

% Gate thresholds for diagnostics, not algorithm parameters.
P.maskAgreementWarn      = 1 - 1e-12;
P.plateauRelativeTol     = 0.02; % 2% relative mean-NRE change 160 -> 320

% Paired bootstrap configuration for the final 10-seed effect table.
P.bootstrapReplicates    = 10000;
P.bootstrapSeed          = 2026091802;

%% ------------------------------------------------------------------------
% 2. Output folder
% -------------------------------------------------------------------------
if isempty(resumeRoot)
    outBase = fullfile(projectRoot,'results_HCL_DSP_FINAL');
    if exist(outBase,'dir')~=7, mkdir(outBase); end
    stamp = datestr(now,'yyyymmddTHHMMSS');
    outRoot = fullfile(outBase, ...
        sprintf('PAVIA_RANK_BUDGET_DIAG_V2_1_%s_%s',mode,stamp));
    mkdir(outRoot);
else
    outRoot = char(resumeRoot);
    assert(exist(outRoot,'dir')==7, ...
        'resumeRoot does not exist: %s',outRoot);
end

figRoot = fullfile(outRoot,'figures');
if exist(figRoot,'dir')~=7, mkdir(figRoot); end

fprintf('\n============================================================\n');
fprintf(' PAVIAU RANK-BUDGET DIAG V2.1 | %s\n',mode);
fprintf('============================================================\n');
fprintf('Project root : %s\n',projectRoot);
fprintf('Pavia file   : %s\n',paviaFile);
fprintf('Output root  : %s\n',outRoot);
fprintf('Primary      : rank=%d | budget=%d | tol=%.1e\n', ...
    P.primaryRank,P.primaryBudget,P.relTolExpected);
fprintf('Pilot ranks  : %s\n',mat2str(P.rankGridPilot));
fprintf('DIAG ranks   : %s\n',mat2str(P.rankGridDiag));
fprintf('Budgets      : %s\n',mat2str(P.budgetGrid));
fprintf('Rank audit   : r90=%d | r95=%d | r98=%d | r99=%d\n', ...
    P.rank90,P.rank95,P.rank98,P.rank99);
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
% 3. Load the EXACT same Pavia reference as V3.2
% -------------------------------------------------------------------------
[Xfull,sourceVar] = load_largest_3d_numeric(paviaFile);

assert(ndims(Xfull)==3,'PaviaU source must be numeric 3-D.');
assert(size(Xfull,1)>=max(P.roiRows) && ...
       size(Xfull,2)>=max(P.roiCols) && size(Xfull,3)==103, ...
       'Unexpected PaviaU dimensions.');

rawMin = min(Xfull(:));
rawMax = max(Xfull(:));

Xfull = double(Xfull) ./ P.normalizationMax;
Xfull = min(max(Xfull,0),1);
Xref = Xfull(P.roiRows,P.roiCols,:);
clear Xfull

assert(isequal(size(Xref),[256 256 103]), ...
    'Working reference must be exactly 256x256x103.');
assert(all(isfinite(Xref(:))),'Reference contains non-finite values.');

refSigma = std(Xref(:));
assert(refSigma>0,'Reference has zero variance.');

fprintf('Source variable : %s\n',sourceVar);
fprintf('Raw range       : [%.6g, %.6g]\n',rawMin,rawMax);
fprintf('Reference size  : %s\n',mat2str(size(Xref)));
fprintf('Reference sigma : %.8g\n',refSigma);

%% ------------------------------------------------------------------------
% 4. Reference identity guard + frozen rank-audit landmarks
% -------------------------------------------------------------------------
assert(abs(refSigma-P.referenceSigmaExpected) <= P.referenceSigmaAbsTol, ...
    ['Pavia reference/preprocessing does not match the audited reference. ' ...
     'Observed sigma=%.10g, expected %.10g. Stop rather than reusing the ' ...
     'previous rank landmarks on a different preprocessing.'], ...
    refSigma,P.referenceSigmaExpected);

Trank = fixed_rank_landmark_table(P);
writetable(Trank,fullfile(outRoot,'Pavia_RBD_rank_landmarks.csv'));

fprintf('\nFrozen reference-based rank landmarks:\n');
disp(Trank);

%% ------------------------------------------------------------------------
% 5. Active ranks, budgets, and seeds
% -------------------------------------------------------------------------
switch mode
    case 'SMOKE'
        activeSeeds   = P.seedSmoke;
        activeRanks   = P.primaryRank;
        activeBudgets = P.budgetGrid;

    case 'REPRO'
        activeSeeds   = P.seedsAll;
        activeRanks   = P.primaryRank;
        activeBudgets = P.primaryBudget;

    case 'PILOT'
        activeSeeds   = P.seedsPilot;
        activeRanks   = P.rankGridPilot;
        activeBudgets = P.budgetGrid;

    case 'DIAG'
        activeSeeds   = P.seedsAll;
        activeRanks   = P.rankGridDiag;
        activeBudgets = P.budgetGrid;

    otherwise
        error('Unexpected mode %s',mode);
end

maxRank = min(size(Xref,1),size(Xref,2));
activeRanks = activeRanks(activeRanks>=1 & activeRanks<=maxRank);
activeRanks = unique(round(activeRanks),'stable');
activeBudgets = sort(unique(activeBudgets));

fprintf('\nActive seeds   : %s\n',seed_text(activeSeeds));
fprintf('Active ranks   : %s\n',mat2str(activeRanks));
fprintf('Active budgets : %s\n',mat2str(activeBudgets));
fprintf('Reference landmarks: r90=%d, r95=%d, r98=%d, r99=%d\n', ...
    P.rank90,P.rank95,P.rank98,P.rank99);

%% ------------------------------------------------------------------------
% 7. Exact Mixed-protocol audit for active seeds
% -------------------------------------------------------------------------
Tprotocol = build_mixed_protocol_audit(Xref,activeSeeds,P,refSigma);
writetable(Tprotocol,fullfile(outRoot,'Pavia_RBD_protocol.csv'));

assert(all(Tprotocol.ProtocolPass==1), ...
    'Mixed corruption protocol failed for at least one seed.');
assert(all(Tprotocol.DisjointPass==1), ...
    'Mixed corruption supports are not disjoint.');

%% ------------------------------------------------------------------------
% 8. HCL dependency and final-config preflight
% -------------------------------------------------------------------------
assert(exist('hcl_trpca_final','file')==2, ...
    'hcl_trpca_final.m is required.');
interfaceReport = preflight_hcl_interface(C,P);

%% ------------------------------------------------------------------------
% 9. Resume previous raw rows if requested
% -------------------------------------------------------------------------
rawFile = fullfile(outRoot,'Pavia_RBD_raw.csv');
partialFile = fullfile(outRoot,'Pavia_RBD_raw_PARTIAL.csv');

Tall = table();
if exist(rawFile,'file')==2
    Tall = readtable(rawFile);
elseif exist(partialFile,'file')==2
    Tall = readtable(partialFile);
end

if ~isempty(Tall)
    fprintf('Resume: loaded %d completed rows.\n',height(Tall));
end

%% ------------------------------------------------------------------------
% 10. Rank x budget diagnosis
% -------------------------------------------------------------------------
for iseed = 1:numel(activeSeeds)
    seed = activeSeeds(iseed);
    generatorSeed = seed + 1000*P.mixedScenarioIndex;

    rng(generatorSeed,'twister');
    [Y,truth,proto] = generate_corruption_v3( ...
        Xref,P.mixedScenario,P,refSigma);

    fprintf('\n============================================================\n');
    fprintf('Mixed seed %d | generator %d | contamination %.4f%%\n', ...
        seed,generatorSeed,100*proto.totalRate);
    fprintf('============================================================\n');

    for ir = 1:numel(activeRanks)
        r = activeRanks(ir);

        doneThisRank = false(size(activeBudgets));
        for ib = 1:numel(activeBudgets)
            doneThisRank(ib) = has_completed_row( ...
                Tall,seed,r,activeBudgets(ib));
        end

        if all(doneThisRank)
            fprintf('  rank=%3d : all requested budgets already complete; skip.\n',r);
            continue;
        end

        % The reference frozen mask is the budget-80 mask. If budget 80 was
        % already completed in a previous interrupted session but later
        % budgets are missing, rerun budget 80 only to recover the exact mask
        % in memory. The repeated result is NOT appended.
        baseBudget = min(activeBudgets);
        baseMask = [];
        baseSelFrac = NaN;
        baseRESS = NaN;

        if has_completed_row(Tall,seed,r,baseBudget)
            fprintf('  rank=%3d : recovering budget-%d mask for resume check...\n', ...
                r,baseBudget);
            [~,baseMask,baseSelFrac,baseRESS] = run_one_setting( ...
                Xref,Y,truth,proto,seed,generatorSeed,r,baseBudget,C,P,false);
        end

        for ib = 1:numel(activeBudgets)
            budget = activeBudgets(ib);

            if has_completed_row(Tall,seed,r,budget)
                fprintf('  rank=%3d budget=%3d : already complete; skip row.\n', ...
                    r,budget);
                continue;
            end

            [row,M,selFrac,rESS] = run_one_setting( ...
                Xref,Y,truth,proto,seed,generatorSeed,r,budget,C,P,true);

            if isempty(baseMask)
                baseMask = M;
                baseSelFrac = selFrac;
                baseRESS = rESS;
                maskAgreement = 1;
                selDiff = 0;
                rESSDiff = 0;
            else
                maskAgreement = mean(M(:)==baseMask(:));
                selDiff = selFrac - baseSelFrac;
                rESSDiff = rESS - baseRESS;
            end

            row.MaskAgreementToBudget80 = maskAgreement;
            row.SelectedFractionDiffToBudget80 = selDiff;
            row.rESSDiffToBudget80 = rESSDiff;

            if isempty(Tall)
                Tall = row;
            else
                Tall = [Tall; row]; %#ok<AGROW>
            end

            Tall = dedupe_diag_table(Tall);
            writetable(Tall,partialFile);

            fprintf(['  rank=%3d budget=%3d | NRE %.6f | Det %.4f | ' ...
                     'Mask %.4f | Level %.4f | rel %.3e | tol=%d | ' ...
                     'maskAgree %.12f\n'], ...
                r,budget,row.NRE,row.DetectionF1,row.MaskF1, ...
                row.LevelMacroF1,row.FinalRelChange,row.ToleranceMet, ...
                row.MaskAgreementToBudget80);

            clear M row
        end
    end

    clear Y truth proto
end

%% ------------------------------------------------------------------------
% 11. Final raw and summary tables
% -------------------------------------------------------------------------
Tall = dedupe_diag_table(Tall);
writetable(Tall,rawFile);

Tsum = summarize_rank_budget(Tall);
writetable(Tsum,fullfile(outRoot,'Pavia_RBD_summary.csv'));

Tpub = compact_publication_summary(Tsum);
writetable(Tpub,fullfile(outRoot,'Pavia_RBD_publication_summary.csv'));

Trepro = evaluate_reproduction(Tall,P);
writetable(Trepro,fullfile(outRoot,'Pavia_RBD_reproduction.csv'));

Teffects = summarize_paired_effects(Tall,P);
writetable(Teffects,fullfile(outRoot,'Pavia_RBD_effects.csv'));

Tplateau = summarize_plateau(Tsum,P);
writetable(Tplateau,fullfile(outRoot,'Pavia_RBD_plateau.csv'));

%% ------------------------------------------------------------------------
% 12. Figures
% -------------------------------------------------------------------------
if ~isempty(Tsum)
    make_rank_budget_figures(Tsum,Trank,P,figRoot);
end
if ~isempty(Teffects)
    make_effect_figure(Teffects,P,figRoot);
end

%% ------------------------------------------------------------------------
% 13. Gate report
% -------------------------------------------------------------------------
reportFile = fullfile(outRoot,'Pavia_RBD_gate_report.txt');
write_gate_report(reportFile,Tall,Tsum,Trank,Trepro,Teffects,Tplateau,P, ...
    activeSeeds,activeRanks,activeBudgets);

%% ------------------------------------------------------------------------
% 14. Compact MAT
% -------------------------------------------------------------------------
R = struct();
R.mode = mode;
R.raw = Tall;
R.summary = Tsum;
R.publicationSummary = Tpub;
R.rankLandmarks = Trank;
R.protocol = Tprotocol;
R.reproduction = Trepro;
R.effects = Teffects;
R.plateau = Tplateau;
R.interface = interfaceReport;
R.config = P;
R.activeSeeds = activeSeeds;
R.activeRanks = activeRanks;
R.activeBudgets = activeBudgets;
R.outRoot = outRoot;
R.source = struct( ...
    'paviaFile',paviaFile, ...
    'sourceVariable',sourceVar, ...
    'rawMin',rawMin, ...
    'rawMax',rawMax, ...
    'referenceSigma',refSigma);

save(fullfile(outRoot,'Pavia_RBD_results_small.mat'), ...
    'R','P','Tall','Tsum','Tpub','Trank','Tprotocol','Trepro', ...
        'Teffects','Tplateau','-v7');

fprintf('\n============================================================\n');
fprintf(' PAVIAU RANK-BUDGET DIAG V2.1 COMPLETE | %s\n',mode);
fprintf('============================================================\n');
disp(Tsum);
fprintf('Gate report:\n  %s\n',reportFile);
fprintf('Results root:\n  %s\n',outRoot);

end


%% =========================================================================
% DIAGNOSTIC-SPECIFIC HELPERS
% =========================================================================
function txt = seed_text(seeds)
if numel(seeds)<=5
    txt = mat2str(seeds);
else
    txt = sprintf('%d ... %d (%d seeds)',seeds(1),seeds(end),numel(seeds));
end
end

function T = fixed_rank_landmark_table(P)
%FIXED_RANK_LANDMARK_TABLE
% Previously completed reference-based post-hoc rank audit.
%
% These are diagnostic landmarks, not deployable model-order estimates.
% They are frozen before all recovery runs in this V2 file.

EnergyTarget = P.rankAuditTargets(:);
EstimatedRank = [P.rank90; P.rank95; P.rank98; P.rank99];
RankStride = repmat(P.rankAuditStride,4,1);
NumSampledFourierSlices = repmat(P.rankAuditSlices,4,1);
ReferenceBasedPostHoc = ones(4,1);
Description = repmat({ ...
    'Median per-slice squared-singular-value energy rank over sampled mode-3 Fourier slices'},4,1);

T = table(EnergyTarget,EstimatedRank,RankStride, ...
    NumSampledFourierSlices,ReferenceBasedPostHoc,Description);
end

function r = rank_for_target(T,target)
[delta,idx] = min(abs(T.EnergyTarget-target));
assert(delta<1e-12,'Requested target %.4f not found in rank landmark table.',target);
r = T.EstimatedRank(idx);
end

function T = build_mixed_protocol_audit(Xref,seeds,P,refSigma)
rows = cell(numel(seeds),11);

for i=1:numel(seeds)
    seed = seeds(i);
    generatorSeed = seed + 1000*P.mixedScenarioIndex;
    rng(generatorSeed,'twister');

    [~,~,proto] = generate_corruption_v3( ...
        Xref,P.mixedScenario,P,refSigma);

    rows(i,:) = {seed,generatorSeed, ...
        proto.entryRate,proto.blockRate,proto.sliceRate,proto.totalRate, ...
        proto.numBlocks,proto.numSlices, ...
        proto.nTotal,double(proto.protocolPass),double(proto.disjointPass)};
end

T = cell2table(rows,'VariableNames',{ ...
    'Seed','GeneratorSeed', ...
    'ActualEntryRate','ActualBlockRate','ActualSliceRate','ActualTotalRate', ...
    'NumBlocks','NumSlices','NumTotal','ProtocolPass','DisjointPass'});
end

function tf = has_completed_row(T,seed,r,budget)
if isempty(T) || ~all(ismember({'Seed','Rank','Budget'},T.Properties.VariableNames))
    tf = false;
    return;
end
tf = any(T.Seed==seed & T.Rank==r & T.Budget==budget);
end

function [row,M,selFrac,rESS] = run_one_setting( ...
    Xref,Y,truth,proto,seed,generatorSeed,r,budget,C,P,verbose)

if nargin<12
    verbose = true;
end

t0 = tic;

H = run_hcl_diag_compat(Y,r,C,budget,P);
Xhat = require_xhat(H,'HCL-TRPCA');
verify_hcl_final_floor(H,P.epsRecExpected);

assert(isequal(size(Xhat),size(Xref)), ...
    'HCL returned the wrong tensor size.');
assert(all(isfinite(Xhat(:))), ...
    'HCL returned non-finite reconstruction values.');

met = recovery_metrics_raw(Xref,Xhat);
D = extract_hcl_diagnostics(H,size(Xref),P);

[detP,detR,detF] = binary_metrics(D.detectMask,truth.anyMask);
[maskP,maskR,maskF] = binary_metrics(D.recoveryMask,truth.anyMask);
[levelMacro,levelE,levelB,levelS] = ...
    level_f1(D.predLevel,truth.levelLabel);

runtime = toc(t0);

selFrac = D.selectedFraction;
rESS = D.rESS;
M = D.recoveryMask;

row = table( ...
    seed,generatorSeed,r,budget, ...
    double(r==P.primaryRank && budget==P.primaryBudget), ...
    proto.entryRate,proto.blockRate,proto.sliceRate,proto.totalRate, ...
    met.NRE,met.MPSNR,met.MSSIM,met.SAMdeg, ...
    detP,detR,detF, ...
    maskP,maskR,maskF, ...
    levelMacro,levelE,levelB,levelS, ...
    D.selectedFraction,D.rESS,D.iterations,D.finalRelChange, ...
    double(D.toleranceMet), ...
    NaN,NaN,NaN, ...
    met.OutOfRangeFraction,runtime, ...
    'VariableNames',{ ...
    'Seed','GeneratorSeed','Rank','Budget','IsPrimarySetting', ...
    'ActualEntryRate','ActualBlockRate','ActualSliceRate','ActualTotalRate', ...
    'NRE','MPSNR','MSSIM','SAMdeg', ...
    'DetectionPrecision','DetectionRecall','DetectionF1', ...
    'MaskPrecision','MaskRecall','MaskF1', ...
    'LevelMacroF1','LevelF1Entry','LevelF1Block','LevelF1Slice', ...
    'SelectedFraction','rESS','Iterations','FinalRelChange','ToleranceMet', ...
    'MaskAgreementToBudget80', ...
    'SelectedFractionDiffToBudget80','rESSDiffToBudget80', ...
    'OutOfRangeFraction','RuntimeSec'});

if verbose
    fprintf('    finished in %.1f s\n',runtime);
end

clear H Xhat D
end

function H = run_hcl_diag_compat(Y,r,C,budget,P)
% The only scientific changes relative to the frozen primary HCL config are:
%   cfg.rank / tubalRank = r
%   cfg.maxIterations / maxUpdates / maxIter = budget
% Everything else is inherited from the paper-final frozen configuration.

Ptmp = P;
Ptmp.maxUpdatesExpected = budget;

cfg = build_frozen_hcl_solver_cfg(C,r,Ptmp);

% Explicit invariants.
assert(cfg.warmIterations==10);
assert(abs(cfg.epsWarm-1e-2)<1e-14);
assert(abs(cfg.epsSel -1e-2)<1e-14);
assert(abs(cfg.epsRec -1e-3)<1e-14);
assert(abs(cfg.gamma-0.85)<1e-14);
assert(abs(cfg.tauDetect-0.50)<1e-14);
assert(abs(cfg.tauHigh-0.90)<1e-14);
assert(abs(cfg.fMax-0.32)<1e-14);
assert(abs(cfg.etaMin-0.68)<1e-14);
assert(cfg.maxIterations==budget);

H = hcl_trpca_final(Y,r,cfg);

Xhat = require_xhat(H,'HCL-TRPCA diagnostic');
assert(isequal(size(Xhat),size(Y)) && all(isfinite(Xhat(:))), ...
    'hcl_trpca_final returned an invalid reconstruction.');
end

function T = dedupe_diag_table(T)
if isempty(T)
    return;
end
keys = strcat(string(T.Seed),'_',string(T.Rank),'_',string(T.Budget));
[~,ia] = unique(keys,'last');
ia = sort(ia);
T = T(ia,:);
T = sortrows(T,{'Seed','Rank','Budget'});
end

function S = summarize_rank_budget(T)
if isempty(T)
    S = table();
    return;
end

ranks = unique(T.Rank(:))';
budgets = unique(T.Budget(:))';

metrics = { ...
    'NRE','MPSNR','MSSIM','SAMdeg', ...
    'DetectionF1','MaskF1','LevelMacroF1', ...
    'SelectedFraction','rESS','Iterations','FinalRelChange', ...
    'MaskAgreementToBudget80','OutOfRangeFraction','RuntimeSec'};

rows = {};

for ir=1:numel(ranks)
    for ib=1:numel(budgets)
        idx = T.Rank==ranks(ir) & T.Budget==budgets(ib);
        if ~any(idx)
            continue;
        end

        row = {ranks(ir),budgets(ib),nnz(idx)};
        for j=1:numel(metrics)
            x = T.(metrics{j})(idx);
            x = x(isfinite(x));
            if isempty(x)
                mu=NaN; sd=NaN;
            else
                mu=mean(x); sd=std(x);
            end
            row = [row,{mu,sd}]; %#ok<AGROW>
        end

        tolx = T.ToleranceMet(idx);
        tolx = tolx(isfinite(tolx));
        if isempty(tolx)
            tolRate = NaN;
        else
            tolRate = mean(tolx);
        end

        agr = T.MaskAgreementToBudget80(idx);
        agr = agr(isfinite(agr));
        if isempty(agr)
            minAgr = NaN;
        else
            minAgr = min(agr);
        end

        row = [row,{tolRate,minAgr}]; %#ok<AGROW>
        rows(end+1,:) = row; %#ok<AGROW>
    end
end

names = {'Rank','Budget','N'};
for j=1:numel(metrics)
    names = [names,{[metrics{j} '_Mean'],[metrics{j} '_Std']}]; %#ok<AGROW>
end
names = [names,{'ToleranceMetRate','MaskAgreement_Min'}];

S = cell2table(rows,'VariableNames',names);
S = sortrows(S,{'Rank','Budget'});
end

function T = evaluate_reproduction(raw,P)
idx = raw.Rank==P.primaryRank & raw.Budget==P.primaryBudget;

if ~any(idx)
    T = table(P.primaryRank,P.primaryBudget,0, ...
        NaN,P.paperMixedNREMean,NaN, ...
        NaN,P.paperMixedDetectionF1,NaN, ...
        NaN,P.paperMixedMaskF1,NaN, ...
        NaN,P.paperMixedLevelF1,NaN, ...
        {'NOT_AVAILABLE'}, ...
        'VariableNames',{ ...
        'Rank','Budget','N', ...
        'ObservedNREMean','PaperNREMean','NREAbsDiff', ...
        'ObservedDetectionF1','PaperDetectionF1','DetectionAbsDiff', ...
        'ObservedMaskF1','PaperMaskF1','MaskAbsDiff', ...
        'ObservedLevelF1','PaperLevelF1','LevelAbsDiff', ...
        'Status'});
    return;
end

sub = raw(idx,:);

nre = sub.NRE(isfinite(sub.NRE));
det = sub.DetectionF1(isfinite(sub.DetectionF1));
msk = sub.MaskF1(isfinite(sub.MaskF1));
lvl = sub.LevelMacroF1(isfinite(sub.LevelMacroF1));

nreMu = mean(nre);
detMu = mean(det);
mskMu = mean(msk);
lvlMu = mean(lvl);

dNRE = abs(nreMu-P.paperMixedNREMean);
dDet = abs(detMu-P.paperMixedDetectionF1);
dMsk = abs(mskMu-P.paperMixedMaskF1);
dLvl = abs(lvlMu-P.paperMixedLevelF1);

if height(sub)<numel(P.seedsAll)
    status = 'INCOMPLETE';
elseif dNRE <= P.reproGreenAbsTol && ...
       dDet <= P.reproStructuralTol && ...
       dMsk <= P.reproStructuralTol && ...
       dLvl <= P.reproStructuralTol
    status = 'PASS_GREEN';
elseif dNRE <= P.reproFailAbsTol
    status = 'PASS_AMBER';
else
    status = 'FAIL';
end

T = table(P.primaryRank,P.primaryBudget,height(sub), ...
    nreMu,P.paperMixedNREMean,dNRE, ...
    detMu,P.paperMixedDetectionF1,dDet, ...
    mskMu,P.paperMixedMaskF1,dMsk, ...
    lvlMu,P.paperMixedLevelF1,dLvl, ...
    {status}, ...
    'VariableNames',{ ...
    'Rank','Budget','N', ...
    'ObservedNREMean','PaperNREMean','NREAbsDiff', ...
    'ObservedDetectionF1','PaperDetectionF1','DetectionAbsDiff', ...
    'ObservedMaskF1','PaperMaskF1','MaskAbsDiff', ...
    'ObservedLevelF1','PaperLevelF1','LevelAbsDiff', ...
    'Status'});
end

function T = compact_publication_summary(S)
if isempty(S)
    T = table();
    return;
end

keep = { ...
    'Rank','Budget','N', ...
    'NRE_Mean','NRE_Std', ...
    'MPSNR_Mean','MPSNR_Std', ...
    'SAMdeg_Mean','SAMdeg_Std', ...
    'DetectionF1_Mean','DetectionF1_Std', ...
    'MaskF1_Mean','MaskF1_Std', ...
    'LevelMacroF1_Mean','LevelMacroF1_Std', ...
    'FinalRelChange_Mean','FinalRelChange_Std', ...
    'ToleranceMetRate','MaskAgreement_Min','RuntimeSec_Mean'};

keep = keep(ismember(keep,S.Properties.VariableNames));
T = S(:,keep);
end

function T = summarize_paired_effects(raw,P)
%SUMMARIZE_PAIRED_EFFECTS
% Quantify rank, budget, and rank x budget interaction using paired tensor
% instances. Negative MeanDeltaNRE means the SECOND setting has lower NRE.
%
% For PILOT (N<5), CIs are intentionally left NaN. For DIAG (N=10), a
% deterministic 10,000-resample paired bootstrap CI is reported.

if isempty(raw)
    T = table();
    return;
end

rows = {};
ranks = unique(raw.Rank(:))';
budgets = unique(raw.Budget(:))';

% 1) Budget effect at each rank: 80 -> 320.
for r = ranks
    [seed,x80,x320] = paired_values(raw,r,80,r,320);
    if isempty(seed), continue; end

    delta = x320-x80;
    relGain = 100*(x80-x320)./max(abs(x80),eps);

    [lo,hi] = maybe_bootstrap_ci(delta,P);
    rows(end+1,:) = { ... %#ok<AGROW>
        'BUDGET_80_TO_320',r,80,r,320,numel(seed), ...
        mean(delta),std(delta),lo,hi,mean(relGain),std(relGain)};
end

% 2) Rank effect at each available budget relative to primary r=25.
for b = budgets
    for r = ranks
        if r==P.primaryRank, continue; end

        [seed,x25,xr] = paired_values(raw,P.primaryRank,b,r,b);
        if isempty(seed), continue; end

        delta = xr-x25;
        relGain = 100*(x25-xr)./max(abs(x25),eps);

        [lo,hi] = maybe_bootstrap_ci(delta,P);
        rows(end+1,:) = { ... %#ok<AGROW>
            'RANK_VS_25',P.primaryRank,b,r,b,numel(seed), ...
            mean(delta),std(delta),lo,hi,mean(relGain),std(relGain)};
    end
end

% 3) Difference-in-differences interaction:
%    [(r,320)-(r,80)] - [(25,320)-(25,80)].
for r = ranks
    if r==P.primaryRank, continue; end

    [seeds,a80,a320,b80,b320] = fourway_values( ...
        raw,r,80,r,320,P.primaryRank,80,P.primaryRank,320);

    if isempty(seeds), continue; end

    interaction = (a320-a80) - (b320-b80);
    [lo,hi] = maybe_bootstrap_ci(interaction,P);

    rows(end+1,:) = { ... %#ok<AGROW>
        'RANK_X_BUDGET_DID',P.primaryRank,80,r,320,numel(seeds), ...
        mean(interaction),std(interaction),lo,hi,NaN,NaN};
end

names = { ...
    'EffectType','RankA','BudgetA','RankB','BudgetB','N', ...
    'MeanDeltaNRE','StdDeltaNRE','CI95Low','CI95High', ...
    'MeanRelativeImprovementPct','StdRelativeImprovementPct'};

if isempty(rows)
    % Typed 0-row table for modes (e.g. REPRO) where paired effects are
    % intentionally unavailable.
    T = table( ...
        cell(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',names);
else
    T = cell2table(rows,'VariableNames',names);
end
end

function T = summarize_plateau(S,P)
%SUMMARIZE_PLATEAU
% Compute the 160->320 NRE plateau diagnostic when both budgets exist.
%
% Modes such as REPRO intentionally contain only budget=80. In that case
% there is no plateau quantity to compute, so return a typed 0-row table.
% This is NOT a failed experiment and must not abort postprocessing.

names = { ...
    'Rank','NRE160','NRE320','RelativeChange160to320Pct', ...
    'ToleranceRate160','ToleranceRate320','FinalRelChange320','Plateau2Pct'};

% Typed empty numeric table: downstream expressions such as T.Rank==r
% remain valid even when height(T)==0.
emptyT = array2table(zeros(0,numel(names)),'VariableNames',names);

if isempty(S) || height(S)==0
    T = emptyT;
    return;
end

ranks = unique(S.Rank(:))';
rows = cell(0,numel(names));

for r = ranks
    i160 = S.Rank==r & S.Budget==160;
    i320 = S.Rank==r & S.Budget==320;

    % REPRO has only budget 80 by design, so this branch is expected there.
    if ~any(i160) || ~any(i320)
        continue;
    end

    n160 = S.NRE_Mean(i160);
    n320 = S.NRE_Mean(i320);

    % The summary should have one row per Rank x Budget. Guard explicitly
    % against accidental duplicates before converting to scalars.
    assert(numel(n160)==1 && numel(n320)==1, ...
        'Expected one summary row per Rank x Budget in summarize_plateau.');

    rel = abs(n160-n320)/max(abs(n320),eps);

    tol160 = S.ToleranceMetRate(i160);
    tol320 = S.ToleranceMetRate(i320);
    rel320 = S.FinalRelChange_Mean(i320);

    assert(numel(tol160)==1 && numel(tol320)==1 && numel(rel320)==1, ...
        'Duplicate Rank x Budget summary rows detected.');

    isPlateau = rel <= P.plateauRelativeTol;

    rows(end+1,:) = { ... %#ok<AGROW>
        double(r),double(n160),double(n320),double(100*rel), ...
        double(tol160),double(tol320),double(rel320),double(isPlateau)};
end

if isempty(rows)
    T = emptyT;
else
    T = cell2table(rows,'VariableNames',names);

    % Force all plateau columns to numeric doubles for stable downstream use.
    for j=1:numel(names)
        T.(names{j}) = double(T.(names{j}));
    end
end
end

function [seed,xA,xB] = paired_values(T,rA,bA,rB,bB)
A = T(T.Rank==rA & T.Budget==bA,{'Seed','NRE'});
B = T(T.Rank==rB & T.Budget==bB,{'Seed','NRE'});

if isempty(A) || isempty(B)
    seed=[]; xA=[]; xB=[];
    return;
end

A.Properties.VariableNames{'NRE'}='NREA';
B.Properties.VariableNames{'NRE'}='NREB';
J = innerjoin(A,B,'Keys','Seed');
J = sortrows(J,'Seed');

good = isfinite(J.NREA) & isfinite(J.NREB);
J = J(good,:);

seed = J.Seed;
xA = J.NREA;
xB = J.NREB;
end

function [seed,a1,a2,b1,b2] = fourway_values(T,rA1,bA1,rA2,bA2,rB1,bB1,rB2,bB2)
A1 = T(T.Rank==rA1 & T.Budget==bA1,{'Seed','NRE'});
A2 = T(T.Rank==rA2 & T.Budget==bA2,{'Seed','NRE'});
B1 = T(T.Rank==rB1 & T.Budget==bB1,{'Seed','NRE'});
B2 = T(T.Rank==rB2 & T.Budget==bB2,{'Seed','NRE'});

if isempty(A1) || isempty(A2) || isempty(B1) || isempty(B2)
    seed=[]; a1=[]; a2=[]; b1=[]; b2=[];
    return;
end

A1.Properties.VariableNames{'NRE'}='A1';
A2.Properties.VariableNames{'NRE'}='A2';
B1.Properties.VariableNames{'NRE'}='B1';
B2.Properties.VariableNames{'NRE'}='B2';

J = innerjoin(A1,A2,'Keys','Seed');
J = innerjoin(J,B1,'Keys','Seed');
J = innerjoin(J,B2,'Keys','Seed');
J = sortrows(J,'Seed');

good = isfinite(J.A1) & isfinite(J.A2) & isfinite(J.B1) & isfinite(J.B2);
J = J(good,:);

seed = J.Seed;
a1 = J.A1;
a2 = J.A2;
b1 = J.B1;
b2 = J.B2;
end

function [lo,hi] = maybe_bootstrap_ci(x,P)
x = double(x(:));
x = x(isfinite(x));

if numel(x)<5
    lo=NaN; hi=NaN;
    return;
end

oldRng = rng;
cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(P.bootstrapSeed + numel(x),'twister');

B = P.bootstrapReplicates;
n = numel(x);
means = zeros(B,1);

for ib=1:B
    idx = randi(n,n,1);
    means(ib) = mean(x(idx));
end

means = sort(means);
lo = empirical_percentile(means,0.025);
hi = empirical_percentile(means,0.975);
end

function q = empirical_percentile(sortedX,p)
n = numel(sortedX);
if n==0
    q=NaN;
    return;
end
pos = 1 + (n-1)*p;
lo = floor(pos);
hi = ceil(pos);
if lo==hi
    q = sortedX(lo);
else
    w = pos-lo;
    q = (1-w)*sortedX(lo) + w*sortedX(hi);
end
end

function make_rank_budget_figures(S,Trank,P,figRoot)
ranks = unique(S.Rank(:))';
budgets = unique(S.Budget(:))';
r98 = P.rank98;

% ---- Figure A: the direct rank-budget diagnosis --------------------------
f = figure('Color','w','Visible','off', ...
    'Position',[60 60 1500 470], ...
    'Name','Pavia rank-budget diagnostic');

tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');

nexttile(tl,1);
hold on;
for ib=1:numel(budgets)
    idx = S.Budget==budgets(ib);
    T = sortrows(S(idx,:), 'Rank');
    errorbar(T.Rank,T.NRE_Mean,T.NRE_Std,'-o', ...
        'LineWidth',1.2,'MarkerSize',5, ...
        'DisplayName',sprintf('%d updates',budgets(ib)));
end
xline(P.primaryRank,'--','Primary r=25','LabelVerticalAlignment','bottom');
if any(ranks==P.rank90)
    xline(P.rank90,':',sprintf('r_{90}=%d',P.rank90), ...
        'LabelVerticalAlignment','top');
end
if any(ranks==P.rank95)
    xline(P.rank95,':',sprintf('r_{95}=%d',P.rank95), ...
        'LabelVerticalAlignment','bottom');
end
if any(ranks==P.rank98)
    xline(P.rank98,':',sprintf('r_{98}=%d',P.rank98), ...
        'LabelVerticalAlignment','top');
end
xlabel('Projection rank');
ylabel('NRE');
title('(a) Reconstruction');
grid on;
legend('Location','best');

nexttile(tl,2);
hold on;
for ib=1:numel(budgets)
    idx = S.Budget==budgets(ib);
    T = sortrows(S(idx,:), 'Rank');
    y = max(T.FinalRelChange_Mean,realmin);
    semilogy(T.Rank,y,'-o', ...
        'LineWidth',1.2,'MarkerSize',5, ...
        'DisplayName',sprintf('%d updates',budgets(ib)));
end
set(gca,'YScale','log');
yline(P.relTolExpected,'--',sprintf('tol = %.0e',P.relTolExpected));
xlabel('Projection rank');
ylabel('Final relative change');
title('(b) Stopping behavior');
grid on;

nexttile(tl,3);
hold on;
for ib=1:numel(budgets)
    idx = S.Budget==budgets(ib);
    T = sortrows(S(idx,:), 'Rank');
    plot(T.Rank,100*T.ToleranceMetRate,'-o', ...
        'LineWidth',1.2,'MarkerSize',5, ...
        'DisplayName',sprintf('%d updates',budgets(ib)));
end
xlabel('Projection rank');
ylabel('Tolerance-met rate (%)');
ylim([-2 102]);
title('(c) Convergence frequency');
grid on;

base = fullfile(figRoot,'Pavia_rank_budget_diagnostic');
savefig(f,[base '.fig']);
try
    exportgraphics(f,[base '.png'],'Resolution',300);
    exportgraphics(f,[base '.pdf'],'ContentType','vector');
catch
    print(f,[base '.png'],'-dpng','-r300');
end
close(f);

% ---- Figure B: structural metrics at each budget -------------------------
f2 = figure('Color','w','Visible','off', ...
    'Position',[60 60 820 560], ...
    'Name','Pavia structural stability');

hold on;
bShow = max(budgets);
idx = S.Budget==bShow;
T = sortrows(S(idx,:),'Rank');

errorbar(T.Rank,T.DetectionF1_Mean,T.DetectionF1_Std,'-o', ...
    'LineWidth',1.2,'MarkerSize',5,'DisplayName','Detection F1');
errorbar(T.Rank,T.MaskF1_Mean,T.MaskF1_Std,'-s', ...
    'LineWidth',1.2,'MarkerSize',5,'DisplayName','Frozen-mask F1');
errorbar(T.Rank,T.LevelMacroF1_Mean,T.LevelMacroF1_Std,'-d', ...
    'LineWidth',1.2,'MarkerSize',5,'DisplayName','Level Macro-F1');

xline(P.primaryRank,'--','Primary r=25');
xlabel('Projection rank');
ylabel('Score');
ylim([0 1.02]);
title(sprintf('Structural diagnostics at %d-update budget',bShow));
grid on;
legend('Location','best');

base2 = fullfile(figRoot,'Pavia_rank_structural_stability');
savefig(f2,[base2 '.fig']);
try
    exportgraphics(f2,[base2 '.png'],'Resolution',300);
    exportgraphics(f2,[base2 '.pdf'],'ContentType','vector');
catch
    print(f2,[base2 '.png'],'-dpng','-r300');
end
close(f2);
end

function make_effect_figure(T,P,figRoot)
%MAKE_EFFECT_FIGURE
% Compact diagnostic visualization of the paired NRE effects.

if isempty(T)
    return;
end

f = figure('Color','w','Visible','off', ...
    'Position',[80 80 1350 430], ...
    'Name','Pavia rank-budget paired effects');

tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');

% (a) Budget gain 80 -> 320 for each rank.
nexttile(tl,1);
A = T(strcmp(T.EffectType,'BUDGET_80_TO_320'),:);
if ~isempty(A)
    A = sortrows(A,'RankA');
    plot(A.RankA,A.MeanRelativeImprovementPct,'-o','LineWidth',1.2);
    yline(0,'--');
    xlabel('Projection rank');
    ylabel('Relative NRE improvement (%)');
    title('(a) Budget effect: 80 \rightarrow 320');
    grid on;
end

% (b) Rank effect at budget 320 relative to r=25.
nexttile(tl,2);
B = T(strcmp(T.EffectType,'RANK_VS_25') & T.BudgetA==320 & T.BudgetB==320,:);
if ~isempty(B)
    B = sortrows(B,'RankB');
    x = [P.primaryRank; B.RankB];
    y = [0; B.MeanRelativeImprovementPct];
    plot(x,y,'-o','LineWidth',1.2);
    yline(0,'--');
    xlabel('Projection rank');
    ylabel('Relative NRE improvement vs r=25 (%)');
    title('(b) Rank effect at 320 updates');
    grid on;
end

% (c) Rank x budget difference-in-differences.
nexttile(tl,3);
C = T(strcmp(T.EffectType,'RANK_X_BUDGET_DID'),:);
if ~isempty(C)
    C = sortrows(C,'RankB');
    plot(C.RankB,C.MeanDeltaNRE,'-o','LineWidth',1.2);
    yline(0,'--');
    xlabel('Projection rank');
    ylabel('Difference-in-differences in NRE');
    title('(c) Rank \times budget interaction');
    grid on;
end

base = fullfile(figRoot,'Pavia_rank_budget_effects');
savefig(f,[base '.fig']);
try
    exportgraphics(f,[base '.png'],'Resolution',300);
    exportgraphics(f,[base '.pdf'],'ContentType','vector');
catch
    print(f,[base '.png'],'-dpng','-r300');
end
close(f);
end

function write_gate_report(filename,raw,S,Trank,Trepro,Teffects,Tplateau,P,seeds,ranks,budgets)
fid = fopen(filename,'w');
assert(fid>0,'Could not open gate report for writing.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'PAVIAU RANK-BUDGET DIAGNOSTIC V2.1 REPORT\n');
fprintf(fid,'Generated: %s\n\n',datestr(now));

fprintf(fid,'SCIENTIFIC PURPOSE\n');
fprintf(fid,['Failure-mechanism diagnosis only. The frozen V3.2 primary Pavia ' ...
    'comparison is not replaced by this diagnostic.\n\n']);

fprintf(fid,'ACTIVE DESIGN\n');
fprintf(fid,'Seeds   : %s\n',seed_text(seeds));
fprintf(fid,'Ranks   : %s\n',mat2str(ranks));
fprintf(fid,'Budgets : %s\n\n',mat2str(budgets));

fprintf(fid,'FROZEN REFERENCE-BASED RANK LANDMARKS\n');
for i=1:height(Trank)
    fprintf(fid,'Target %.2f -> median sampled-slice energy rank %d\n', ...
        Trank.EnergyTarget(i),Trank.EstimatedRank(i));
end
fprintf(fid,['These are post-hoc reference-based diagnostic landmarks, not true ' ...
    'tubal ranks and not deployable rank estimates.\n\n']);

fprintf(fid,'GATE D0: PRIMARY-SETTING REPRODUCTION\n');
if ~isempty(Trepro)
    fprintf(fid,['r=%d, budget=%d, N=%d\n' ...
        '  NRE: observed %.6f | paper %.6f | abs diff %.6f\n' ...
        '  DetectionF1: observed %.4f | paper %.4f | abs diff %.4f\n' ...
        '  MaskF1: observed %.4f | paper %.4f | abs diff %.4f\n' ...
        '  LevelMacroF1: observed %.4f | paper %.4f | abs diff %.4f\n' ...
        '  Status: %s\n\n'], ...
        Trepro.Rank(1),Trepro.Budget(1),Trepro.N(1), ...
        Trepro.ObservedNREMean(1),Trepro.PaperNREMean(1),Trepro.NREAbsDiff(1), ...
        Trepro.ObservedDetectionF1(1),Trepro.PaperDetectionF1(1),Trepro.DetectionAbsDiff(1), ...
        Trepro.ObservedMaskF1(1),Trepro.PaperMaskF1(1),Trepro.MaskAbsDiff(1), ...
        Trepro.ObservedLevelF1(1),Trepro.PaperLevelF1(1),Trepro.LevelAbsDiff(1), ...
        Trepro.Status{1});
end

fprintf(fid,'GATE D1: FROZEN-MASK INVARIANCE ACROSS BUDGETS\n');
agr = raw.MaskAgreementToBudget80;
agr = agr(isfinite(agr));
if isempty(agr)
    fprintf(fid,'No agreement values available.\n\n');
else
    fprintf(fid,'Minimum observed mask agreement to budget 80: %.12f\n',min(agr));
    if min(agr) >= P.maskAgreementWarn
        fprintf(fid,'Status: PASS. Budget does not change the frozen mask.\n\n');
    else
        fprintf(fid,['Status: WARNING. The selected mask changed with budget; ' ...
            'the budget comparison is not clean.\n\n']);
    end
end

fprintf(fid,'GATE D2/D3: CONVERGENCE AND NRE PLATEAU\n');
for ir=1:numel(ranks)
    r=ranks(ir);
    fprintf(fid,'\nRank %d\n',r);

    for ib=1:numel(budgets)
        b=budgets(ib);
        idx=S.Rank==r & S.Budget==b;
        if any(idx)
            q=S(idx,:);
            fprintf(fid,['  budget %d: NRE %.6f +/- %.6f | rel %.3e | ' ...
                'tol rate %.1f%% | DetectionF1 %.4f | MaskF1 %.4f | LevelF1 %.4f\n'], ...
                b,q.NRE_Mean,q.NRE_Std,q.FinalRelChange_Mean, ...
                100*q.ToleranceMetRate,q.DetectionF1_Mean, ...
                q.MaskF1_Mean,q.LevelMacroF1_Mean);
        end
    end

    ip = Tplateau.Rank==r;
    if any(ip)
        q=Tplateau(ip,:);
        fprintf(fid,'  relative NRE change 160->320: %.3f%%', ...
            q.RelativeChange160to320Pct);
        if q.Plateau2Pct==1
            fprintf(fid,'  [PLATEAU <= %.1f%%]\n',100*P.plateauRelativeTol);
        else
            fprintf(fid,'  [NOT PLATEAU]\n');
        end
    end
end

fprintf(fid,'\nPAIRED EFFECT DECOMPOSITION\n');
if isempty(Teffects)
    fprintf(fid,'No complete paired effect table available yet.\n');
else
    B = Teffects(strcmp(Teffects.EffectType,'BUDGET_80_TO_320'),:);
    for i=1:height(B)
        fprintf(fid,['Budget effect r=%d, 80->320: mean delta NRE %.6f, ' ...
            'relative improvement %.2f%%'], ...
            B.RankA(i),B.MeanDeltaNRE(i),B.MeanRelativeImprovementPct(i));
        if isfinite(B.CI95Low(i))
            fprintf(fid,', 95%% bootstrap CI [%.6f, %.6f]', ...
                B.CI95Low(i),B.CI95High(i));
        end
        fprintf(fid,'\n');
    end

    R = Teffects(strcmp(Teffects.EffectType,'RANK_VS_25') & ...
        Teffects.BudgetA==320 & Teffects.BudgetB==320,:);
    for i=1:height(R)
        fprintf(fid,['Rank effect at 320, r=25 -> r=%d: mean delta NRE %.6f, ' ...
            'relative improvement %.2f%%'], ...
            R.RankB(i),R.MeanDeltaNRE(i),R.MeanRelativeImprovementPct(i));
        if isfinite(R.CI95Low(i))
            fprintf(fid,', 95%% bootstrap CI [%.6f, %.6f]', ...
                R.CI95Low(i),R.CI95High(i));
        end
        fprintf(fid,'\n');
    end
end

fprintf(fid,'\nINTERPRETATION RULES\n');
fprintf(fid,['1) A materially lower NRE at larger budgets with an invariant frozen ' ...
    'mask identifies a finite-budget component.\n']);
fprintf(fid,['2) A materially lower NRE at larger ranks after 160->320 NRE has ' ...
    'plateaued identifies a model-order component.\n']);
fprintf(fid,['3) If convergence/plateau is reached and reconstruction remains weak ' ...
    'across the pre-specified rank grid, retain natural-data model/recovery ' ...
    'mismatch as the remaining explanation.\n']);
fprintf(fid,['4) Detection-F1, mask-F1, and level Macro-F1 must be read separately ' ...
    'from NRE. Stable structural metrics with changing NRE directly support ' ...
    'the attribution/recovery separation claimed in the manuscript.\n']);
fprintf(fid,['5) Do NOT report the best diagnostic setting as independent performance. ' ...
    'A changed primary protocol requires fresh Pavia confirmation seeds.\n']);
end




%% =========================================================================
% FROZEN V3.2 HELPERS COPIED VERBATIM FROM THE CURRENT PRIMARY RUNNER
% =========================================================================
function s = mk_scenario(name,e,b,sl)
s = struct('name',name,'entryRate',e,'blockRate',b,'sliceRate',sl);
end

function paviaFile = resolve_pavia_file(C,projectRoot)
candidates = {};

if isstruct(C) && isfield(C,'real')
    if isfield(C.real,'paviaMat') && ~isempty(C.real.paviaMat)
        candidates{end+1} = C.real.paviaMat; %#ok<AGROW>
    end
    if isfield(C.real,'dataRoot') && ~isempty(C.real.dataRoot)
        candidates{end+1} = fullfile(C.real.dataRoot,'PaviaU.mat'); %#ok<AGROW>
    end
end

% Public path fallbacks only; local-machine absolute paths are intentionally
% excluded from the release. projectRoot is <repository>/code/pavia here.
repoRoot = fileparts(fileparts(projectRoot));
candidates = [candidates, { ...
    fullfile(repoRoot,'data','PaviaU.mat')}];

paviaFile = '';
for i=1:numel(candidates)
    f = char(candidates{i});
    if exist(f,'file')==2
        paviaFile = f;
        break;
    end
end
assert(~isempty(paviaFile), ...
    'Could not locate PaviaU.mat. Edit resolve_pavia_file() or C.real.paviaMat.');
end

function [X,varName] = load_largest_3d_numeric(matFile)
info = whos('-file',matFile);
bestBytes = -inf;
varName = '';

numericClasses = {'double','single','uint8','uint16','uint32','uint64', ...
                  'int8','int16','int32','int64'};

for i=1:numel(info)
    if numel(info(i).size)==3 && ismember(info(i).class,numericClasses)
        if info(i).bytes > bestBytes
            bestBytes = info(i).bytes;
            varName = info(i).name;
        end
    end
end

assert(~isempty(varName), ...
    'No numeric 3-D tensor was found in %s.',matFile);

S = load(matFile,varName);
X = S.(varName);
end


%% =========================================================================
% CONTROLLED CORRUPTION V3
% =========================================================================
function [Y,truth,proto] = generate_corruption_v3(X,sc,P,refSigma)

sz = size(X);
n1 = sz(1); n2 = sz(2); n3 = sz(3);
N = numel(X);

entryMask = false(sz);
blockMask = false(sz);
sliceMask = false(sz);

% ---- 1) Whole-band / slice corruption first -------------------------------
ns = round(sc.sliceRate*n3);
if sc.sliceRate>0 && ns<1
    ns = 1;
end
ns = min(ns,n3);

if ns>0
    sliceIds = sort(randperm(n3,ns));
    sliceMask(:,:,sliceIds) = true;
else
    sliceIds = [];
end

% ---- 2) Non-overlapping local spatial blocks on non-slice bands ----------
targetBlock = round(sc.blockRate*N);
numBlocks = 0;
attempts = 0;

if targetBlock>0
    validBands = setdiff(1:n3,sliceIds);
    assert(~isempty(validBands),'No non-slice bands available for block corruption.');

    % 10%% occupancy is moderate; this cap is intentionally generous.
    maxAttempts = 200000;

    while nnz(blockMask) < targetBlock
        attempts = attempts + 1;
        if attempts > maxAttempts
            error(['Block placement failed before reaching the requested rate. ' ...
                   'Placed %.4f%%, target %.4f%%.'], ...
                   100*nnz(blockMask)/N,100*sc.blockRate);
        end

        b = P.blockSizes(randi(numel(P.blockSizes)));
        b = min([b,n1,n2]);

        k = validBands(randi(numel(validBands)));
        i0 = randi(n1-b+1);
        j0 = randi(n2-b+1);

        region = blockMask(i0:i0+b-1,j0:j0+b-1,k);
        if any(region(:))
            continue;
        end

        % Entire square is accepted; the final block may create a tiny
        % overshoot. We report the realized rate rather than trimming a block.
        blockMask(i0:i0+b-1,j0:j0+b-1,k) = true;
        numBlocks = numBlocks + 1;
    end
end

% ---- 3) Isolated entry corruption on still-clean support ------------------
targetEntry = round(sc.entryRate*N);
if targetEntry>0
    avail = find(~sliceMask & ~blockMask);
    assert(numel(avail)>=targetEntry, ...
        'Not enough clean support remains for the requested entry corruption.');
    pick = avail(randperm(numel(avail),targetEntry));
    entryMask(pick) = true;
end

% ---- 4) Protocol validity -------------------------------------------------
disjointPass = ...
    ~any(entryMask(:) & blockMask(:)) && ...
    ~any(entryMask(:) & sliceMask(:)) && ...
    ~any(blockMask(:) & sliceMask(:));
assert(disjointPass,'Entry/block/slice supports must be disjoint.');

anyMask = entryMask | blockMask | sliceMask;

levelLabel = zeros(sz,'uint8');
levelLabel(entryMask) = 1;
levelLabel(blockMask) = 2;
levelLabel(sliceMask) = 3;

% ---- 5) Additive gross corruption ----------------------------------------
E = zeros(sz,'double');
E = add_signed_sigma_corruption(E,entryMask,refSigma,P);
E = add_signed_sigma_corruption(E,blockMask,refSigma,P);
E = add_signed_sigma_corruption(E,sliceMask,refSigma,P);

Y = X + E;
assert(all(isfinite(Y(:))),'Generated observation contains non-finite values.');

% ---- 6) Exact realized rates ----------------------------------------------
proto = struct();
proto.entryRate = nnz(entryMask)/N;
proto.blockRate = nnz(blockMask)/N;
proto.sliceRate = nnz(sliceMask)/N;
proto.totalRate = nnz(anyMask)/N;

proto.nEntry = nnz(entryMask);
proto.nBlock = nnz(blockMask);
proto.nSlice = nnz(sliceMask);
proto.nTotal = nnz(anyMask);
proto.numBlocks = numBlocks;
proto.numSlices = numel(sliceIds);
proto.disjointPass = disjointPass;

% Entry is exact up to integer rounding.
tolEntry = 1/N + 10*eps;

% Slice rate is constrained by integer bands.
tolSlice = 0.5/n3 + 10*eps;

% Block is allowed to overshoot by at most one largest square.
maxBlockOvershoot = max(P.blockSizes)^2/N + 10*eps;

passE = abs(proto.entryRate-sc.entryRate) <= tolEntry;
passS = abs(proto.sliceRate-sc.sliceRate) <= tolSlice;

if sc.blockRate==0
    passB = proto.blockRate==0;
else
    passB = proto.blockRate >= sc.blockRate - 10*eps && ...
            proto.blockRate <= sc.blockRate + maxBlockOvershoot;
end

proto.protocolPass = passE && passB && passS && disjointPass;

truth = struct();
truth.entryMask = entryMask;
truth.blockMask = blockMask;
truth.sliceMask = sliceMask;
truth.anyMask = anyMask;
truth.levelLabel = levelLabel;
truth.sliceIds = sliceIds;
truth.actualEntryRate = proto.entryRate;
truth.actualBlockRate = proto.blockRate;
truth.actualSliceRate = proto.sliceRate;
truth.actualTotalRate = proto.totalRate;
end

function E = add_signed_sigma_corruption(E,M,refSigma,P)
n = nnz(M);
if n==0
    return;
end
sgn = ones(n,1);
sgn(rand(n,1)<0.5) = -1;

mult = P.ampSigmaMin + ...
       (P.ampSigmaMax-P.ampSigmaMin).*rand(n,1);

E(M) = sgn .* refSigma .* mult;
end


%% =========================================================================
% FINAL METHOD ADAPTERS
% =========================================================================
function report = preflight_hcl_interface(C,P)
%PREFLIGHT_HCL_INTERFACE
% Fail fast on an interface/config mismatch BEFORE any expensive Pavia run.

hfile = which('hcl_trpca_final');
assert(~isempty(hfile),'hcl_trpca_final.m is not on the MATLAB path.');

fprintf('\n============================================================\n');
fprintf(' HCL FINAL INTERFACE PREFLIGHT\n');
fprintf('============================================================\n');
fprintf('Solver file : %s\n',hfile);

cfg = build_frozen_hcl_solver_cfg(C,3,P);

[requiredFields,arg3] = third_argument_fields(hfile);
fprintf('Confirmed call: hcl_trpca_final(Y,r,%s)\n',arg3);
fprintf('Direct cfg fields detected in source: %d\n',numel(requiredFields));

% Tiny deterministic tensor. This is an interface/numerical sanity check,
% NOT an experimental result and is never pooled with Pavia statistics.
oldRng = rng;
cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(2026091520,'twister');

Ytoy = randn(20,18,9);
Ytoy = local_rank_project_for_preflight(Ytoy,3);
Ytoy(4:8,5:10,3) = Ytoy(4:8,5:10,3) + 2.5;
Ytoy(3,4,7) = Ytoy(3,4,7) - 3.0;

t0 = tic;
Htoy = hcl_trpca_final(Ytoy,3,cfg);
elapsed = toc(t0);
Xtoy = require_xhat(Htoy,'HCL-TRPCA preflight');

assert(isequal(size(Xtoy),size(Ytoy)), ...
    'Preflight output size mismatch.');
assert(all(isfinite(Xtoy(:))), ...
    'Preflight output contains non-finite values.');

fprintf('Preflight tiny call: PASS (%.3f s)\n',elapsed);
fprintf('Frozen solver cfg:\n');
fprintf('  warmIterations = %d\n',cfg.warmIterations);
fprintf('  maxIterations  = %d\n',cfg.maxIterations);
fprintf('  epsWarm        = %.4g\n',cfg.epsWarm);
fprintf('  epsSel         = %.4g\n',cfg.epsSel);
fprintf('  epsRec         = %.4g\n',cfg.epsRec);
fprintf('  gamma          = %.4g\n',cfg.gamma);
fprintf('  relTol/tol     = %.4g / %.4g\n',cfg.relTol,cfg.tol);
fprintf('============================================================\n');

report = struct();
report.pass = true;
report.solverFile = hfile;
report.thirdArgumentName = arg3;
report.requiredFields = requiredFields;
report.cfg = cfg;
report.toyRuntimeSec = elapsed;
end

function cfg = build_frozen_hcl_solver_cfg(C,r,P)
%BUILD_FROZEN_HCL_SOLVER_CFG
% Construct the solver cfg explicitly rather than guessing a nested project
% path. Values below are the already-frozen FINAL method specification.
%
% Important: this function repairs ONLY the runner/config interface.
% It does not alter the corruption generator, Pavia seeds, thresholds, or
% algorithm parameters in response to Pavia outcomes.

cfg = struct();

% ---- Core iteration/stopping ---------------------------------------------
cfg.warmIterations = 10;
cfg.warmStartIter   = 10;       % compatibility alias
cfg.maxIterations  = P.maxUpdatesExpected;
cfg.maxUpdates     = P.maxUpdatesExpected;
cfg.maxIter        = P.maxUpdatesExpected;
cfg.tol            = P.relTolExpected;
cfg.relTol         = P.relTolExpected;
cfg.tolerance      = P.relTolExpected;
cfg.gamma          = 0.85;
cfg.damping        = 0.85;

% ---- Rank aliases --------------------------------------------------------
cfg.rank            = r;
cfg.tubalRank       = r;

% ---- Warm residual / hierarchy geometry ---------------------------------
cfg.b               = 8;
cfg.blockSize       = 8;
cfg.warmBlockSize   = 8;
cfg.h               = 5;
cfg.attrWindow      = 5;
cfg.attributionWindow = 5;

% ---- Selection -----------------------------------------------------------
cfg.tauDetect       = 0.50;
cfg.tau_d           = 0.50;
cfg.detectProbThreshold = 0.50;
cfg.detectionThreshold  = 0.50;

cfg.tauHigh         = 0.90;
cfg.tau_h           = 0.90;
cfg.highThreshold   = 0.90;
cfg.selectionThreshold = 0.90;

cfg.fMax            = 0.32;
cfg.maxSelectedFraction = 0.32;
cfg.trimCap         = 0.32;

cfg.etaMin          = 0.68;
cfg.essMin          = 0.68;
cfg.minESS          = 0.68;
cfg.essLowerBound   = 0.68;

% ---- Three distinct floors: DO NOT collapse to one weightFloor -----------
cfg.epsWarm         = 1e-2;
cfg.warmFloor       = 1e-2;
cfg.epsSel          = 1e-2;
cfg.selectionFloor  = 1e-2;
cfg.epsRec          = P.epsRecExpected;
cfg.recoveryFloor   = P.epsRecExpected;

% ---- Evidence enhancement ------------------------------------------------
cfg.betaBlock       = 0.65;
cfg.betaSlice       = 0.75;
cfg.beta_b          = 0.65;
cfg.beta_s          = 0.75;
cfg.nu              = 0.50;

% ---- Evidence calibration ------------------------------------------------
cfg.entryEvidenceCenter = 2.5;
cfg.evidenceCenterEntry = 2.5;
cfg.cEvidenceEntry      = 2.5;
cfg.entryEvidenceScale  = 0.45;
cfg.evidenceScaleEntry  = 0.45;
cfg.entryTemperature    = 0.45;

cfg.blockEvidenceCenter = 0.45;
cfg.evidenceCenterBlock = 0.45;
cfg.blockEvidenceScale  = 0.06;
cfg.evidenceScaleBlock  = 0.06;

cfg.sliceEvidenceCenter = 0.50;
cfg.evidenceCenterSlice = 0.50;
cfg.sliceEvidenceScale  = 0.06;
cfg.evidenceScaleSlice  = 0.06;

% ---- Warm Tukey constants ------------------------------------------------
cfg.tukeyEntry      = 4.685;
cfg.tukeyBlock      = 4.0;
cfg.tukeySlice      = 3.5;
cfg.cEntry          = 4.685;
cfg.cBlock          = 4.0;
cfg.cSlice          = 3.5;

% ---- Scale / numerical conventions --------------------------------------
cfg.epsScale        = 1e-8;
cfg.scaleEps        = 1e-8;
cfg.epsilonScale    = 1e-8;

% ---- Harmless diagnostics/control flags ---------------------------------
cfg.verbose         = false;
cfg.showProgress    = false;
cfg.saveHistory     = false;
cfg.storeHistory    = false;
cfg.returnDiagnostics = true;

% Copy same-named extra fields from the project configuration ONLY when the
% frozen map above does not already define them. This is useful for purely
% implementation-specific fields (e.g. FFT options), but it never overwrites
% a frozen scientific parameter.
hfile = which('hcl_trpca_final');
[requiredFields,~] = third_argument_fields(hfile);

for i=1:numel(requiredFields)
    f = requiredFields{i};
    if ~isfield(cfg,f)
        [found,val,pathText] = find_field_recursive(C,f,'C'); %#ok<ASGLU>
        if found
            cfg.(f) = val;
        end
    end
end

% A single old-style weightFloor cannot represent the FINAL method because
% warm/selection floors are 1e-2 while recovery floor is 1e-3.
% If the local "final" solver still requires only weightFloor and does not
% expose the separate final fields, stop rather than silently changing the
% method.
if any(strcmp(requiredFields,'weightFloor')) && ...
        ~all(ismember({'epsWarm','epsSel','epsRec'},requiredFields))
    error(['The local hcl_trpca_final.m still appears to use the legacy ' ...
           'single field cfg.weightFloor. That interface cannot faithfully ' ...
           'represent epsWarm=epsSel=1e-2 and epsRec=1e-3. Use the same ' ...
           'final hcl_trpca_final.m that produced the Urban/SBI final runs.']);
end

missing = requiredFields(~isfield(cfg,requiredFields));
if ~isempty(missing)
    fprintf('\nThe local final solver directly references cfg fields that the\n');
    fprintf('runner cannot safely infer:\n');
    for i=1:numel(missing)
        fprintf('  - %s\n',missing{i});
    end
    error(['HCL solver-config preflight failed BEFORE the Pavia audit. ' ...
           'Upload hcl_trpca_final.m and HCL_DSP_FINAL_config.m if these ' ...
           'implementation-specific fields need to be mapped.']);
end

% Scientific invariant check.
assert(abs(cfg.epsWarm-1e-2)<1e-14);
assert(abs(cfg.epsSel -1e-2)<1e-14);
assert(abs(cfg.epsRec -1e-3)<1e-14);
assert(cfg.warmIterations==10);
assert(cfg.maxIterations==P.maxUpdatesExpected);
assert(abs(cfg.gamma-0.85)<1e-14);
assert(abs(cfg.tauDetect-0.50)<1e-14);
assert(abs(cfg.tauHigh-0.90)<1e-14);
assert(abs(cfg.fMax-0.32)<1e-14);
assert(abs(cfg.etaMin-0.68)<1e-14);
end

function [found,val,pathText] = find_field_recursive(S,target,path0)
found = false;
val = [];
pathText = '';

if ~isstruct(S) || numel(S)~=1
    return;
end

if isfield(S,target)
    found = true;
    val = S.(target);
    pathText = [path0 '.' target];
    return;
end

fn = fieldnames(S);
for i=1:numel(fn)
    v = S.(fn{i});
    if isstruct(v) && numel(v)==1
        [found,val,pathText] = find_field_recursive(v,target,[path0 '.' fn{i}]);
        if found
            return;
        end
    end
end
end

function X = local_rank_project_for_preflight(Y,r)
% Small self-contained truncated t-SVD projection for interface testing.
Yf = fft(double(Y),[],3);
Xf = zeros(size(Yf),'like',Yf);
for k=1:size(Yf,3)
    [U,S,V] = svd(Yf(:,:,k),'econ');
    rr = min([r,size(S,1),size(S,2)]);
    Xf(:,:,k) = U(:,1:rr)*S(1:rr,1:rr)*V(:,1:rr)';
end
X = real(ifft(Xf,[],3));
end

function [fields,arg3] = third_argument_fields(hfile)
% Parse cfg.<field> accesses for the third argument of hcl_trpca_final.
fields = {};
arg3 = '';
try
    txt = fileread(hfile);
catch ME
    error('Could not read local solver source %s: %s',hfile,ME.message);
end

tok = regexp(txt, ...
    'function[^\n\r]*hcl_trpca_final\s*\(([^\)]*)\)', ...
    'tokens','once');
if isempty(tok)
    error('Could not parse hcl_trpca_final function declaration.');
end

args = strtrim(strsplit(tok{1},','));
assert(numel(args)>=3,'hcl_trpca_final must have at least three inputs.');
arg3 = strtrim(args{3});

pat = [regexptranslate('escape',arg3) '\.([A-Za-z]\w*)'];
t = regexp(txt,pat,'tokens');
if ~isempty(t)
    fields = unique(cellfun(@(x)x{1},t,'UniformOutput',false),'stable');
end
end

function Xhat = require_xhat(O,method)
if isnumeric(O)
    Xhat = O;
    return;
end

assert(isstruct(O), ...
    '%s output must be a tensor or struct.',method);

candidates = {'Xhat','X','L','reconstruction','Xrec'};
Xhat = [];
for i=1:numel(candidates)
    if isfield(O,candidates{i}) && isnumeric(O.(candidates{i})) ...
            && ndims(O.(candidates{i}))==3
        Xhat = O.(candidates{i});
        break;
    end
end

assert(~isempty(Xhat), ...
    '%s output does not contain Xhat/X/L/reconstruction/Xrec.',method);
Xhat = double(Xhat);
end

function verify_hcl_final_floor(H,expected)
% Final Pavia must use epsRec=1e-3 while warm/selection remain at 1e-2.

verified = false;

if isfield(H,'epsRec') && isscalar(H.epsRec)
    assert(abs(double(H.epsRec)-expected)<=1e-12, ...
        'HCL returned epsRec=%g, expected %g.',double(H.epsRec),expected);
    verified = true;
end

if isfield(H,'W0') && isnumeric(H.W0) && ~isempty(H.W0)
    W = double(H.W0);
    low = W(W<0.5);
    if ~isempty(low)
        floorVal = median(low(:));
        assert(abs(floorVal-expected)<=1e-10, ...
            'Frozen HCL recovery weights use floor %g, expected %g.',floorVal,expected);
        verified = true;
    end
end

if ~verified
    warning(['Could not independently verify epsRec from HCL output fields. ' ...
             'Confirm hcl_trpca_final is the same final core used by Urban/SBI.']);
end
end

function D = extract_hcl_diagnostics(H,sz,P)

[piE,piB,piS] = final_attribution_maps(H,sz);

po = min(1,max(0,piE+piB+piS));
detectMask = po >= P.detectThreshold;

stack = cat(4,piE,piB,piS);
[~,predLevel] = max(stack,[],4);
predLevel = uint8(predLevel);
predLevel(~detectMask) = 0;

recoveryMask = find_recovery_mask(H,sz);

D = struct();
D.piEntry = piE;
D.piBlock = piB;
D.piSlice = piS;
D.detectMask = detectMask;
D.predLevel = predLevel;
D.recoveryMask = recoveryMask;
D.selectedFraction = mean(recoveryMask(:));

D.rESS = get_scalar_field(H,{'recoveryESS','rESS','ESSRatio'},NaN);
if ~isfinite(D.rESS)
    W = find_weight_tensor(H,sz);
    if ~isempty(W)
        D.rESS = ess_ratio(W);
    else
        % Exact for the final two-level frozen mask.
        f = D.selectedFraction;
        e = P.epsRecExpected;
        D.rESS = (1-f+e*f)^2 / max(1-f+e^2*f,eps);
    end
end

D.iterations = get_scalar_field(H,{'iterations','iter','numIterations'},NaN);
D.finalRelChange = get_scalar_field(H, ...
    {'finalRelChange','relChange','finalRelativeChange'},NaN);

tolField = get_scalar_field(H,{'toleranceMet','converged'},NaN);
if isfinite(tolField)
    D.toleranceMet = logical(tolField);
elseif isfinite(D.finalRelChange)
    D.toleranceMet = D.finalRelChange < P.relTolExpected;
else
    D.toleranceMet = false;
end
end

function [piE,piB,piS] = final_attribution_maps(H,sz)

piE=[]; piB=[]; piS=[];

if isfield(H,'final') && isstruct(H.final)
    F = H.final;
    if isfield(F,'piEntry'), piE=F.piEntry; end
    if isfield(F,'piBlock'), piB=F.piBlock; end
    if isfield(F,'piSlice'), piS=F.piSlice; end
end

if isempty(piE) && isfield(H,'pEntry'), piE=H.pEntry; end
if isempty(piB) && isfield(H,'pBlock'), piB=H.pBlock; end
if isempty(piS) && isfield(H,'pSlice'), piS=H.pSlice; end

assert(~isempty(piE) && ~isempty(piB) && ~isempty(piS), ...
    ['Final HCL entry/block/slice attribution maps were not found. ' ...
     'V3 requires them for the paper-final Pavia attribution test.']);

assert(isequal(size(piE),sz) && isequal(size(piB),sz) && isequal(size(piS),sz), ...
    'Final HCL attribution maps have incorrect sizes.');

piE = double(piE);
piB = double(piB);
piS = double(piS);
end

function M = find_recovery_mask(H,sz)
M=[];

direct = {'recoveryMask','frozenMask','mask','O0'};
for i=1:numel(direct)
    if isfield(H,direct{i}) && isequal(size(H.(direct{i})),sz)
        M = logical(H.(direct{i}));
        return;
    end
end

if isfield(H,'selection') && isstruct(H.selection)
    S=H.selection;
    cand={'recoveryMask','frozenMask','mask','O0'};
    for i=1:numel(cand)
        if isfield(S,cand{i}) && isequal(size(S.(cand{i})),sz)
            M=logical(S.(cand{i}));
            return;
        end
    end
end

W = find_weight_tensor(H,sz);
if ~isempty(W)
    M = W < 0.5;
    return;
end

error(['Could not locate the ACTUAL frozen recovery mask in the HCL output. ' ...
       'Do not reconstruct it from a score threshold. Add the local output ' ...
       'field name to find_recovery_mask().']);
end

function W = find_weight_tensor(H,sz)
W=[];
candidates={'W0','W','weights','recoveryWeights'};
for i=1:numel(candidates)
    if isfield(H,candidates{i}) && isnumeric(H.(candidates{i})) ...
            && isequal(size(H.(candidates{i})),sz)
        W=double(H.(candidates{i}));
        return;
    end
end
if isfield(H,'selection') && isstruct(H.selection)
    S=H.selection;
    candidates={'W0','W','weights','recoveryWeights'};
    for i=1:numel(candidates)
        if isfield(S,candidates{i}) && isnumeric(S.(candidates{i})) ...
                && isequal(size(S.(candidates{i})),sz)
            W=double(S.(candidates{i}));
            return;
        end
    end
end
end

function v = get_scalar_field(S,names,defaultValue)
v=defaultValue;
for i=1:numel(names)
    if isfield(S,names{i}) && isnumeric(S.(names{i})) ...
            && isscalar(S.(names{i}))
        v=double(S.(names{i}));
        return;
    end
end
end

function e = ess_ratio(W)
w=double(W(:));
e=(sum(w)^2)/max(numel(w)*sum(w.^2),eps);
end


%% =========================================================================
% METRICS -- RAW RECONSTRUCTION, NO SILENT BOX PROJECTION
% =========================================================================
function met = recovery_metrics_raw(Xref,Xhat)
A=double(Xref);
B=double(Xhat);

d=B-A;
met.NRE=norm(d(:))/max(norm(A(:)),eps);
met.MPSNR=mean_psnr_raw(A,B);
met.MSSIM=mean_ssim_raw(A,B);
met.SAMdeg=mean_sam_deg(A,B);
met.OutOfRangeFraction=mean(B(:)<0 | B(:)>1);
end

function p = mean_psnr_raw(A,B)
K=size(A,3);
v=nan(K,1);
for k=1:K
    D=A(:,:,k)-B(:,:,k);
    mse=mean(D(:).^2);
    if mse<=eps
        v(k)=Inf;
    else
        v(k)=10*log10(1/mse); % operational reference range is [0,1]
    end
end
if any(isinf(v))
    p=Inf;
else
    p=mean(v,'omitnan');
end
end

function s = mean_ssim_raw(A,B)
if exist('ssim','file')~=2
    warning('ssim() unavailable; MSSIM will be NaN.');
    s=NaN;
    return;
end

K=size(A,3);
v=nan(K,1);
for k=1:K
    v(k)=ssim(B(:,:,k),A(:,:,k),'DynamicRange',1);
end
s=mean(v,'omitnan');
end

function a = mean_sam_deg(A,B)
P=reshape(A,[],size(A,3));
Q=reshape(B,[],size(B,3));

dotv=sum(P.*Q,2);
den=sqrt(sum(P.^2,2).*sum(Q.^2,2));
good=den>1e-12;

if ~any(good)
    a=NaN;
    return;
end

c=dotv(good)./den(good);
c=max(-1,min(1,c));
a=mean(acos(c))*180/pi;
end

function [p,r,f] = binary_metrics(pred,truth)
pred=logical(pred(:));
truth=logical(truth(:));

tp=nnz(pred & truth);
fp=nnz(pred & ~truth);
fn=nnz(~pred & truth);

if ~any(truth)
    p=NaN; r=NaN; f=NaN;
    return;
end

p=tp/max(tp+fp,1);
r=tp/max(tp+fn,1);
f=2*p*r/max(p+r,eps);
end

function [macro,e,b,s] = level_f1(pred,trueL)
pred=uint8(pred(:));
trueL=uint8(trueL(:));

vals=nan(3,1);
for c=1:3
    t=(trueL==c);
    if ~any(t)
        continue;
    end
    q=(pred==c);
    tp=nnz(t&q);
    fp=nnz(~t&q);
    fn=nnz(t&~q);
    pp=tp/max(tp+fp,1);
    rr=tp/max(tp+fn,1);
    vals(c)=2*pp*rr/max(pp+rr,eps);
end

present=vals(isfinite(vals));
if isempty(present)
    macro=NaN;
else
    macro=mean(present);
end

e=vals(1); b=vals(2); s=vals(3);
end