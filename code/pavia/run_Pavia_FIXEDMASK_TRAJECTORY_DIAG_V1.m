function R = run_Pavia_FIXEDMASK_TRAJECTORY_DIAG_V1(mode)
% RUN_PAVIA_FIXEDMASK_TRAJECTORY_DIAG_V1
%
% Pavia Mixed fixed-mask recovery trajectory diagnosis.
%
% PURPOSE
% -------
% Diagnose why longer frozen-weight recovery can increase clean-reference
% NRE even when the recovery iteration continues to reduce its own fitting
% criterion.
%
% This script does NOT retune HCL.
% It does NOT modify the primary Pavia protocol.
% It does NOT replace any result in the primary paper table.
%
% Logged trajectory after HCL warm start + mask selection:
%
%   1) Frozen weighted objective
%      f_W0(X) = 0.5 * || sqrt(W0).*(Y-X) ||_F^2
%
%   2) Clean-reference NRE
%
%   3) Relative update
%
%   4) Eq.(20)-based normalized rank-feasibility upper bound
%
% MODES
% -----
% SMOKE:
%   seed 20280001, rank 25, total budget 80
%
% PILOT:
%   seed 20280001, ranks [25 30], total budget 320
%
% FINAL:
%   seeds 20280001:20280010, ranks [25 30], total budget 320
%
% RECOMMENDED:
%
%   R0 = run_Pavia_FIXEDMASK_TRAJECTORY_DIAG_V1('SMOKE');
%   R1 = run_Pavia_FIXEDMASK_TRAJECTORY_DIAG_V1('PILOT');
%   R2 = run_Pavia_FIXEDMASK_TRAJECTORY_DIAG_V1('FINAL');
%
% IMPORTANT
% ---------
% Total update count follows the paper convention:
% 10 warm-start iterations are INCLUDED in Tmax.
%
% The diagnostic reproduces ONLY the frozen-weight recovery stage after
% obtaining the official warm state / frozen mask from hcl_trpca_final.
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(mode)
    mode = 'SMOKE';
end

mode = upper(char(string(mode)));

assert(ismember(mode,{'SMOKE','PILOT','FINAL'}), ...
    'mode must be SMOKE, PILOT, or FINAL.');

%% ========================================================================
% 0. Project / frozen settings
% =========================================================================

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));

assert(exist('HCL_DSP_FINAL_config','file')==2, ...
    'HCL_DSP_FINAL_config.m not found.');

assert(exist('hcl_trpca_final','file')==2, ...
    'hcl_trpca_final.m not found.');

C = HCL_DSP_FINAL_config();

P = struct();

P.version = 'PAVIA_FIXEDMASK_TRAJECTORY_DIAG_V1';

P.roiRows = 178:433;
P.roiCols = 43:298;
P.normalizationMax = 8000;

P.blockSizes = [8 12 16 24];

P.ampSigmaMin = 6;
P.ampSigmaMax = 10;

P.entryRate = 0.05;
P.blockRate = 0.10;
P.sliceRate = 0.10;

P.warmIterations = 10;

P.gamma = 0.85;
P.epsWarm = 1e-2;
P.epsSel  = 1e-2;
P.epsRec  = 1e-3;

P.relTol = 1e-5;

switch mode

    case 'SMOKE'

        P.seeds = 20280001;
        P.ranks = 25;
        P.maxTotalUpdates = 80;

    case 'PILOT'

        P.seeds = 20280001;
        P.ranks = [25 30];
        P.maxTotalUpdates = 320;

    case 'FINAL'

        P.seeds = 20280001:20280010;
        P.ranks = [25 30];
        P.maxTotalUpdates = 320;
end

P.checkpoints = intersect( ...
    [10 20 40 80 160 320], ...
    10:P.maxTotalUpdates);

%% ========================================================================
% 1. Output
% =========================================================================

stamp = datestr(now,'yyyymmddTHHMMSS');

outBase = fullfile(projectRoot,'results_HCL_DSP_FINAL');

if exist(outBase,'dir')~=7
    mkdir(outBase);
end

outRoot = fullfile(outBase, ...
    sprintf('PAVIA_FIXEDMASK_TRAJECTORY_V1_%s_%s',mode,stamp));

mkdir(outRoot);

figRoot = fullfile(outRoot,'figures');
mkdir(figRoot);

fprintf('\n============================================================\n');
fprintf(' PAVIA FIXED-MASK TRAJECTORY DIAG V1 | %s\n',mode);
fprintf('============================================================\n');
fprintf('Seeds       : %s\n',mat2str(P.seeds));
fprintf('Ranks       : %s\n',mat2str(P.ranks));
fprintf('Warm updates: %d\n',P.warmIterations);
fprintf('Max budget  : %d total updates\n',P.maxTotalUpdates);
fprintf('gamma       : %.4f\n',P.gamma);
fprintf('epsRec      : %.4g\n',P.epsRec);
fprintf('Output      : %s\n',outRoot);
fprintf('============================================================\n\n');

%% ========================================================================
% 2. Load Pavia reference
% =========================================================================

paviaFile = resolve_pavia_file_diag(C,projectRoot);

[Xfull,varName] = load_largest_3d_numeric_diag(paviaFile);

Xfull = double(Xfull) ./ P.normalizationMax;
Xfull = min(max(Xfull,0),1);

Xref = Xfull(P.roiRows,P.roiCols,:);

clear Xfull

assert(isequal(size(Xref),[256 256 103]), ...
    'Expected working Pavia tensor 256x256x103.');

refSigma = std(Xref(:));

fprintf('Pavia file      : %s\n',paviaFile);
fprintf('Source variable : %s\n',varName);
fprintf('Reference size  : %s\n',mat2str(size(Xref)));
fprintf('Reference sigma : %.8f\n\n',refSigma);

%% ========================================================================
% 3. Run trajectories
% =========================================================================

allRows = {};
runSummary = {};

runCounter = 0;

for iseed = 1:numel(P.seeds)

    seed = P.seeds(iseed);

    % Mixed is scenario #4 in final V3.2.
    % Preserve exact final generator mapping:
    generatorSeed = seed + 4000;

    rng(generatorSeed,'twister');

    [Y,truth,proto] = generate_mixed_v3_diag( ...
        Xref,P,refSigma);

    fprintf('\n============================================================\n');
    fprintf('Mixed seed %d | generator %d | contamination %.4f%%\n', ...
        seed,generatorSeed,100*proto.totalRate);
    fprintf('============================================================\n');

    for ir = 1:numel(P.ranks)

        r = P.ranks(ir);

        runCounter = runCounter + 1;

        fprintf('\nRank %d\n',r);

        %% -----------------------------------------------------------------
        % A. Obtain official warm-start state and frozen recovery mask.
        %
        % We call the same frozen HCL implementation with total budget
        % equal to the warm-start budget. No scientific parameter is changed.
        % ------------------------------------------------------------------

        cfgWarm = build_frozen_cfg_diag( ...
            C,r,P,P.warmIterations);

        t0 = tic;

        H0 = hcl_trpca_final(Y,r,cfgWarm);

        warmRuntime = toc(t0);

        X = require_xhat_diag(H0);

        W0 = find_weight_tensor_diag(H0,size(Y));

        if isempty(W0)

            M0 = find_recovery_mask_diag(H0,size(Y));

            W0 = ones(size(Y));
            W0(M0) = P.epsRec;
        end

        W0 = double(W0);

        assert(isequal(size(W0),size(Y)), ...
            'W0 size mismatch.');

        lowWeights = W0(W0 < 0.5);

        assert(~isempty(lowWeights), ...
            'No selected recovery locations were found.');

        observedFloor = median(lowWeights(:));

        assert(abs(observedFloor-P.epsRec)<1e-10, ...
            'Recovery floor is %.6g, expected %.6g.', ...
            observedFloor,P.epsRec);

        selectedFraction = mean(W0(:)<0.5);

        fprintf('  warm state obtained in %.1f s\n',warmRuntime);
        fprintf('  selected fraction = %.6f\n',selectedFraction);

        %% -----------------------------------------------------------------
        % B. Allocate trajectory.
        % ------------------------------------------------------------------

        Tmax = P.maxTotalUpdates;

        obj = nan(Tmax,1);
        nre = nan(Tmax,1);
        rel = nan(Tmax,1);
        rankBound = nan(Tmax,1);

        % warm-end state = total update index Tw
        tWarm = P.warmIterations;

        obj(tWarm) = weighted_objective_diag(Y,X,W0);
        nre(tWarm) = nre_diag(Xref,X);

        obj0 = obj(tWarm);

        fprintf('  t=%3d | obj %.6e | NRE %.6f\n', ...
            tWarm,obj(tWarm),nre(tWarm));

        %% -----------------------------------------------------------------
        % C. Frozen-weight recovery trajectory.
        % ------------------------------------------------------------------

        recoveryStart = tic;

        stoppedEarly = false;
        actualFinalT = tWarm;

        for t = tWarm:(Tmax-1)

            V = W0 .* Y + (1-W0) .* X;

            Z = tubal_rank_project_diag(V,r);

            Xnew = ...
                P.gamma .* Z + ...
                (1-P.gamma) .* X;

            stepNorm = norm(Xnew(:)-X(:));

            relChange = ...
                stepNorm / ...
                max(norm(X(:)),eps);

            % Eq.(20)-based NORMALIZED upper bound
            rkBound = ...
                ((1-P.gamma)/P.gamma) * ...
                stepNorm / ...
                max(norm(Xnew(:)),eps);

            X = Xnew;

            totalT = t+1;

            obj(totalT) = weighted_objective_diag(Y,X,W0);
            nre(totalT) = nre_diag(Xref,X);
            rel(totalT) = relChange;
            rankBound(totalT) = rkBound;

            actualFinalT = totalT;

            if ismember(totalT,P.checkpoints)

                fprintf(['  t=%3d | obj/obj10 %.6f | NRE %.6f | ' ...
                         'rel %.3e | rankBound %.3e\n'], ...
                    totalT, ...
                    obj(totalT)/obj0, ...
                    nre(totalT), ...
                    rel(totalT), ...
                    rankBound(totalT));
            end

            if relChange < P.relTol

                fprintf('  tolerance reached at total t=%d\n',totalT);

                stoppedEarly = true;
                break;
            end
        end

        recoveryRuntime = toc(recoveryStart);

        %% -----------------------------------------------------------------
        % D. Convert to table.
        % ------------------------------------------------------------------

        validT = (tWarm:actualFinalT)';

        objNorm = obj(validT) ./ obj0;

        for jj=1:numel(validT)

            tt = validT(jj);

            allRows(end+1,:) = { ... %#ok<AGROW>
                seed, ...
                generatorSeed, ...
                r, ...
                tt, ...
                obj(tt), ...
                objNorm(jj), ...
                nre(tt), ...
                rel(tt), ...
                rankBound(tt), ...
                selectedFraction, ...
                proto.totalRate};
        end

        %% -----------------------------------------------------------------
        % E. Run-level diagnosis.
        % ------------------------------------------------------------------

        use = validT(validT>tWarm);

        if isempty(use)

            objMonotoneFraction = NaN;
            objEndRatio = NaN;

        else

            dObj = diff(obj(tWarm:actualFinalT));

            objMonotoneFraction = ...
                mean(dObj <= max(1e-12,1e-10*abs(obj0)));

            objEndRatio = obj(actualFinalT)/obj0;
        end

        [bestNRE,bestIdx] = ...
            min(nre(tWarm:actualFinalT));

        bestT = tWarm + bestIdx - 1;

        finalNRE = nre(actualFinalT);

        runSummary(end+1,:) = { ... %#ok<AGROW>
            seed, ...
            generatorSeed, ...
            r, ...
            proto.totalRate, ...
            selectedFraction, ...
            actualFinalT, ...
            stoppedEarly, ...
            obj0, ...
            obj(actualFinalT), ...
            objEndRatio, ...
            objMonotoneFraction, ...
            bestNRE, ...
            bestT, ...
            finalNRE, ...
            rel(actualFinalT), ...
            rankBound(actualFinalT), ...
            warmRuntime, ...
            recoveryRuntime};

        fprintf(['  DONE rank=%d | final t=%d | best NRE %.6f @%d | ' ...
                 'final NRE %.6f | obj ratio %.6f | mono %.4f\n'], ...
            r,actualFinalT,bestNRE,bestT, ...
            finalNRE,objEndRatio,objMonotoneFraction);

        %% -----------------------------------------------------------------
        % F. Save one run immediately.
        % ------------------------------------------------------------------

        runFile = fullfile(outRoot, ...
            sprintf('trajectory_seed%d_rank%d.mat',seed,r));

        Trace = struct();

        Trace.seed = seed;
        Trace.generatorSeed = generatorSeed;
        Trace.rank = r;
        Trace.totalUpdate = validT;
        Trace.weightedObjective = obj(validT);
        Trace.objectiveNormalized = objNorm;
        Trace.NRE = nre(validT);
        Trace.relativeUpdate = rel(validT);
        Trace.rankFeasibilityBound = rankBound(validT);
        Trace.selectedFraction = selectedFraction;
        Trace.actualContamination = proto.totalRate;
        Trace.W0Floor = observedFloor;
        Trace.actualFinalT = actualFinalT;

        save(runFile,'Trace','-v7');

        clear H0 W0 X Xnew V Z Trace
    end

    clear Y truth proto
end

%% ========================================================================
% 4. Build CSV tables
% =========================================================================

trajNames = { ...
    'Seed', ...
    'GeneratorSeed', ...
    'Rank', ...
    'TotalUpdate', ...
    'WeightedObjective', ...
    'ObjectiveNormalized', ...
    'NRE', ...
    'RelativeUpdate', ...
    'RankFeasibilityBound', ...
    'SelectedFraction', ...
    'ActualContamination'};

T = cell2table(allRows,'VariableNames',trajNames);

summaryNames = { ...
    'Seed', ...
    'GeneratorSeed', ...
    'Rank', ...
    'ActualContamination', ...
    'SelectedFraction', ...
    'FinalTotalUpdate', ...
    'ToleranceMet', ...
    'ObjectiveAtWarmEnd', ...
    'ObjectiveFinal', ...
    'ObjectiveFinalRatio', ...
    'ObjectiveMonotoneFraction', ...
    'BestNRE', ...
    'BestNREUpdate', ...
    'FinalNRE', ...
    'FinalRelativeUpdate', ...
    'FinalRankFeasibilityBound', ...
    'WarmRuntimeSec', ...
    'RecoveryRuntimeSec'};

S = cell2table(runSummary,'VariableNames',summaryNames);

writetable(T, ...
    fullfile(outRoot,'Pavia_fixedmask_trajectory_raw.csv'));

writetable(S, ...
    fullfile(outRoot,'Pavia_fixedmask_trajectory_run_summary.csv'));

%% ========================================================================
% 5. Aggregate by rank and iteration
% =========================================================================

A = aggregate_trajectory_diag(T,P);

writetable(A, ...
    fullfile(outRoot,'Pavia_fixedmask_trajectory_aggregate.csv'));

%% ========================================================================
% 6. Generate paper diagnostic figure
% =========================================================================

make_trajectory_figure_diag(A,P,figRoot);

%% ========================================================================
% 7. Automatic scientific interpretation table
% =========================================================================

D = diagnose_behavior_diag(S);

writetable(D, ...
    fullfile(outRoot,'Pavia_fixedmask_diagnosis.csv'));

%% ========================================================================
% 8. Save
% =========================================================================

R = struct();

R.mode = mode;
R.config = P;
R.raw = T;
R.runSummary = S;
R.aggregate = A;
R.diagnosis = D;
R.outRoot = outRoot;

save(fullfile(outRoot, ...
    'Pavia_fixedmask_trajectory_results.mat'), ...
    'R','P','T','S','A','D','-v7.3');

fprintf('\n============================================================\n');
fprintf(' FIXED-MASK TRAJECTORY DIAGNOSTIC COMPLETE\n');
fprintf('============================================================\n');
disp(S);
disp(D);
fprintf('\nResults:\n%s\n',outRoot);

end


%% ========================================================================
% HELPERS
% =========================================================================

function cfg = build_frozen_cfg_diag(C,r,P,maxIterations)

cfg = struct();

cfg.rank = r;
cfg.tubalRank = r;

cfg.warmIterations = P.warmIterations;
cfg.warmStartIter = P.warmIterations;

cfg.maxIterations = maxIterations;
cfg.maxUpdates = maxIterations;
cfg.maxIter = maxIterations;

cfg.tol = P.relTol;
cfg.relTol = P.relTol;
cfg.tolerance = P.relTol;

cfg.gamma = P.gamma;
cfg.damping = P.gamma;

cfg.b = 8;
cfg.blockSize = 8;
cfg.warmBlockSize = 8;

cfg.h = 5;
cfg.attrWindow = 5;
cfg.attributionWindow = 5;

cfg.tauDetect = 0.50;
cfg.tau_d = 0.50;
cfg.detectProbThreshold = 0.50;
cfg.detectionThreshold = 0.50;

cfg.tauHigh = 0.90;
cfg.tau_h = 0.90;
cfg.highThreshold = 0.90;
cfg.selectionThreshold = 0.90;

cfg.fMax = 0.32;
cfg.maxSelectedFraction = 0.32;
cfg.trimCap = 0.32;

cfg.etaMin = 0.68;
cfg.essMin = 0.68;
cfg.minESS = 0.68;
cfg.essLowerBound = 0.68;

cfg.epsWarm = P.epsWarm;
cfg.warmFloor = P.epsWarm;

cfg.epsSel = P.epsSel;
cfg.selectionFloor = P.epsSel;

cfg.epsRec = P.epsRec;
cfg.recoveryFloor = P.epsRec;

cfg.betaBlock = 0.65;
cfg.betaSlice = 0.75;
cfg.beta_b = 0.65;
cfg.beta_s = 0.75;
cfg.nu = 0.50;

cfg.entryEvidenceCenter = 2.5;
cfg.evidenceCenterEntry = 2.5;
cfg.cEvidenceEntry = 2.5;

cfg.entryEvidenceScale = 0.45;
cfg.evidenceScaleEntry = 0.45;
cfg.entryTemperature = 0.45;

cfg.blockEvidenceCenter = 0.45;
cfg.evidenceCenterBlock = 0.45;

cfg.blockEvidenceScale = 0.06;
cfg.evidenceScaleBlock = 0.06;

cfg.sliceEvidenceCenter = 0.50;
cfg.evidenceCenterSlice = 0.50;

cfg.sliceEvidenceScale = 0.06;
cfg.evidenceScaleSlice = 0.06;

cfg.tukeyEntry = 4.685;
cfg.tukeyBlock = 4.0;
cfg.tukeySlice = 3.5;

cfg.cEntry = 4.685;
cfg.cBlock = 4.0;
cfg.cSlice = 3.5;

cfg.epsScale = 1e-8;
cfg.scaleEps = 1e-8;
cfg.epsilonScale = 1e-8;

cfg.verbose = false;
cfg.showProgress = false;
cfg.returnDiagnostics = true;

% Copy implementation-specific fields if required.
solverFile = which('hcl_trpca_final');

if ~isempty(solverFile)

    [requiredFields,~] = third_argument_fields_diag(solverFile);

    for i=1:numel(requiredFields)

        f = requiredFields{i};

        if ~isfield(cfg,f)

            [found,val] = find_field_recursive_diag(C,f);

            if found
                cfg.(f) = val;
            end
        end
    end
end

end


function [Y,truth,proto] = generate_mixed_v3_diag(X,P,refSigma)

sz = size(X);

n1 = sz(1);
n2 = sz(2);
n3 = sz(3);

N = numel(X);

entryMask = false(sz);
blockMask = false(sz);
sliceMask = false(sz);

%% slice first

ns = round(P.sliceRate*n3);

sliceIds = sort(randperm(n3,ns));

sliceMask(:,:,sliceIds) = true;

%% block next

targetBlock = round(P.blockRate*N);

validBands = setdiff(1:n3,sliceIds);

numBlocks = 0;
attempts = 0;

while nnz(blockMask) < targetBlock

    attempts = attempts + 1;

    assert(attempts < 200000, ...
        'Block generation exceeded maximum attempts.');

    b = P.blockSizes(randi(numel(P.blockSizes)));

    k = validBands(randi(numel(validBands)));

    i0 = randi(n1-b+1);
    j0 = randi(n2-b+1);

    region = blockMask( ...
        i0:i0+b-1, ...
        j0:j0+b-1, ...
        k);

    if any(region(:))
        continue;
    end

    blockMask( ...
        i0:i0+b-1, ...
        j0:j0+b-1, ...
        k) = true;

    numBlocks = numBlocks + 1;
end

%% entry last

targetEntry = round(P.entryRate*N);

avail = find(~sliceMask & ~blockMask);

pick = avail(randperm(numel(avail),targetEntry));

entryMask(pick) = true;

%% disjoint

assert(~any(entryMask(:)&blockMask(:)));
assert(~any(entryMask(:)&sliceMask(:)));
assert(~any(blockMask(:)&sliceMask(:)));

anyMask = entryMask | blockMask | sliceMask;

%% additive corruption

E = zeros(sz);

E = add_corruption_diag(E,entryMask,refSigma,P);
E = add_corruption_diag(E,blockMask,refSigma,P);
E = add_corruption_diag(E,sliceMask,refSigma,P);

Y = X + E;

%% outputs

truth = struct();

truth.entryMask = entryMask;
truth.blockMask = blockMask;
truth.sliceMask = sliceMask;
truth.anyMask = anyMask;

proto = struct();

proto.entryRate = nnz(entryMask)/N;
proto.blockRate = nnz(blockMask)/N;
proto.sliceRate = nnz(sliceMask)/N;
proto.totalRate = nnz(anyMask)/N;

proto.numBlocks = numBlocks;
proto.numSlices = numel(sliceIds);

end


function E = add_corruption_diag(E,M,refSigma,P)

n = nnz(M);

if n==0
    return;
end

sgn = ones(n,1);

sgn(rand(n,1)<0.5) = -1;

mag = ...
    P.ampSigmaMin + ...
    (P.ampSigmaMax-P.ampSigmaMin).*rand(n,1);

E(M) = sgn .* refSigma .* mag;

end


function X = tubal_rank_project_diag(Y,r)

if exist('hcl_project_tubal_rank','file')==2

    X = hcl_project_tubal_rank(Y,r);
    return;
end

Yf = fft(double(Y),[],3);

Xf = zeros(size(Yf),'like',Yf);

for k=1:size(Yf,3)

    [U,S,V] = svd(Yf(:,:,k),'econ');

    rr = min([r,size(S,1),size(S,2)]);

    Xf(:,:,k) = ...
        U(:,1:rr) * ...
        S(1:rr,1:rr) * ...
        V(:,1:rr)';
end

X = real(ifft(Xf,[],3));

end


function f = weighted_objective_diag(Y,X,W)

R = Y-X;

f = 0.5 * sum(W(:).*R(:).^2);

end


function v = nre_diag(Xref,X)

v = norm(X(:)-Xref(:)) / max(norm(Xref(:)),eps);

end


function X = require_xhat_diag(H)

if isnumeric(H)

    X = double(H);
    return;
end

cand = {'Xhat','X','L','reconstruction','Xrec'};

X = [];

for i=1:numel(cand)

    if isfield(H,cand{i}) && ...
            isnumeric(H.(cand{i})) && ...
            ndims(H.(cand{i}))==3

        X = double(H.(cand{i}));
        return;
    end
end

error('Could not find reconstructed tensor in HCL output.');

end


function W = find_weight_tensor_diag(H,sz)

W = [];

cand = {'W0','W','weights','recoveryWeights'};

for i=1:numel(cand)

    if isfield(H,cand{i}) && ...
            isnumeric(H.(cand{i})) && ...
            isequal(size(H.(cand{i})),sz)

        W = double(H.(cand{i}));
        return;
    end
end

if isfield(H,'selection') && isstruct(H.selection)

    S = H.selection;

    for i=1:numel(cand)

        if isfield(S,cand{i}) && ...
                isnumeric(S.(cand{i})) && ...
                isequal(size(S.(cand{i})),sz)

            W = double(S.(cand{i}));
            return;
        end
    end
end

end


function M = find_recovery_mask_diag(H,sz)

M = [];

cand = {'recoveryMask','frozenMask','mask','O0'};

for i=1:numel(cand)

    if isfield(H,cand{i}) && ...
            isequal(size(H.(cand{i})),sz)

        M = logical(H.(cand{i}));
        return;
    end
end

if isfield(H,'selection') && isstruct(H.selection)

    S = H.selection;

    for i=1:numel(cand)

        if isfield(S,cand{i}) && ...
                isequal(size(S.(cand{i})),sz)

            M = logical(S.(cand{i}));
            return;
        end
    end
end

error('Could not locate frozen HCL recovery mask.');

end


function A = aggregate_trajectory_diag(T,P)

rows = {};

for ir = 1:numel(P.ranks)

    r = P.ranks(ir);

    Tr = T(T.Rank==r,:);

    tVals = unique(Tr.TotalUpdate);

    for j=1:numel(tVals)

        tt = tVals(j);

        idx = Tr.TotalUpdate==tt;

        n = nnz(idx);

        rows(end+1,:) = { ... %#ok<AGROW>
            r, ...
            tt, ...
            n, ...
            mean(Tr.ObjectiveNormalized(idx),'omitnan'), ...
            std(Tr.ObjectiveNormalized(idx),'omitnan'), ...
            mean(Tr.NRE(idx),'omitnan'), ...
            std(Tr.NRE(idx),'omitnan'), ...
            mean(Tr.RelativeUpdate(idx),'omitnan'), ...
            std(Tr.RelativeUpdate(idx),'omitnan'), ...
            mean(Tr.RankFeasibilityBound(idx),'omitnan'), ...
            std(Tr.RankFeasibilityBound(idx),'omitnan')};
    end
end

A = cell2table(rows,'VariableNames',{ ...
    'Rank', ...
    'TotalUpdate', ...
    'N', ...
    'ObjectiveNormMean', ...
    'ObjectiveNormSD', ...
    'NREMean', ...
    'NRESD', ...
    'RelativeUpdateMean', ...
    'RelativeUpdateSD', ...
    'RankBoundMean', ...
    'RankBoundSD'});

end


function make_trajectory_figure_diag(A,P,figRoot)

f = figure( ...
    'Color','w', ...
    'Position',[100 100 1200 850], ...
    'Visible','off');

tl = tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% objective

nexttile;
hold on;

for r=P.ranks

    B=A(A.Rank==r,:);

    plot(B.TotalUpdate, ...
         B.ObjectiveNormMean, ...
         '-o', ...
         'LineWidth',1.5, ...
         'MarkerSize',3);
end

xlabel('Total update index');
ylabel('Normalized weighted objective');
title('(a) Frozen weighted objective');
grid on;

%% NRE

nexttile;
hold on;

for r=P.ranks

    B=A(A.Rank==r,:);

    plot(B.TotalUpdate, ...
         B.NREMean, ...
         '-o', ...
         'LineWidth',1.5, ...
         'MarkerSize',3);
end

xlabel('Total update index');
ylabel('NRE');
title('(b) Clean-reference recovery error');
grid on;

%% relative update

nexttile;
hold on;

for r=P.ranks

    B=A(A.Rank==r,:);

    semilogy(B.TotalUpdate, ...
         B.RelativeUpdateMean, ...
         '-o', ...
         'LineWidth',1.5, ...
         'MarkerSize',3);
end

yline(P.relTol,'--');

xlabel('Total update index');
ylabel('Relative update');
title('(c) Stopping diagnostic');
grid on;

%% rank bound

nexttile;
hold on;

for r=P.ranks

    B=A(A.Rank==r,:);

    semilogy(B.TotalUpdate, ...
         B.RankBoundMean, ...
         '-o', ...
         'LineWidth',1.5, ...
         'MarkerSize',3);
end

xlabel('Total update index');
ylabel('Normalized feasibility bound');
title('(d) Rank-feasibility upper bound');
grid on;

legend( ...
    arrayfun(@(x)sprintf('r = %d',x), ...
    P.ranks,'UniformOutput',false), ...
    'Location','best');

title(tl, ...
    'Pavia Mixed: fixed-mask recovery trajectory diagnosis');

base = fullfile(figRoot, ...
    'Pavia_fixedmask_trajectory');

savefig(f,[base '.fig']);

exportgraphics( ...
    f,[base '.png'], ...
    'Resolution',400);

exportgraphics( ...
    f,[base '.pdf'], ...
    'ContentType','vector');

close(f);

end


function D = diagnose_behavior_diag(S)

rows = {};

for i=1:height(S)

    objectiveMostlyDescending = ...
        S.ObjectiveMonotoneFraction(i) >= 0.95;

    nreWorsenedFromBest = ...
        S.FinalNRE(i) > ...
        S.BestNRE(i) * 1.01;

    if objectiveMostlyDescending && nreWorsenedFromBest

        diagnosis = ...
            "OBJECTIVE_REFERENCE_MISMATCH_CANDIDATE";

    elseif ~objectiveMostlyDescending

        diagnosis = ...
            "SOLVER_NONMONOTONE_DIAGNOSTIC";

    else

        diagnosis = ...
            "NO_CLEAR_PATHOLOGY";
    end

    rows(end+1,:) = { ... %#ok<AGROW>
        S.Seed(i), ...
        S.Rank(i), ...
        objectiveMostlyDescending, ...
        nreWorsenedFromBest, ...
        diagnosis};
end

D = cell2table(rows,'VariableNames',{ ...
    'Seed', ...
    'Rank', ...
    'ObjectiveMostlyDescending', ...
    'NREWorsenedAfterBest', ...
    'Diagnosis'});

end


function paviaFile = resolve_pavia_file_diag(C,projectRoot)

candidates = {};

if isfield(C,'real')

    if isfield(C.real,'paviaMat')
        candidates{end+1}=C.real.paviaMat; %#ok<AGROW>
    end

    if isfield(C.real,'dataRoot')
        candidates{end+1}= ...
            fullfile(C.real.dataRoot,'PaviaU.mat'); %#ok<AGROW>
    end
end

% Public path fallback only; local-machine absolute paths are intentionally
% excluded from the release. projectRoot is <repository>/code/pavia here.
repoRoot = fileparts(fileparts(projectRoot));
candidates = [candidates,{ ...
    fullfile(repoRoot,'data','PaviaU.mat')}];

paviaFile='';

for i=1:numel(candidates)

    f=char(candidates{i});

    if exist(f,'file')==2

        paviaFile=f;
        return;
    end
end

error('Could not locate PaviaU.mat.');

end


function [X,varName] = load_largest_3d_numeric_diag(matFile)

info=whos('-file',matFile);

varName='';
best=-Inf;

numericClasses={ ...
    'double','single', ...
    'uint8','uint16','uint32','uint64', ...
    'int8','int16','int32','int64'};

for i=1:numel(info)

    if numel(info(i).size)==3 && ...
            ismember(info(i).class,numericClasses)

        if info(i).bytes>best

            best=info(i).bytes;
            varName=info(i).name;
        end
    end
end

assert(~isempty(varName), ...
    'No numeric 3-D tensor found.');

S=load(matFile,varName);

X=S.(varName);

end


function [fields,arg3] = third_argument_fields_diag(hfile)

txt=fileread(hfile);

tok=regexp(txt, ...
    'function[^\n\r]*hcl_trpca_final\s*\(([^\)]*)\)', ...
    'tokens','once');

if isempty(tok)

    fields={};
    arg3='';
    return;
end

args=strtrim(strsplit(tok{1},','));

arg3=strtrim(args{3});

pat=[regexptranslate('escape',arg3) ...
    '\.([A-Za-z]\w*)'];

t=regexp(txt,pat,'tokens');

if isempty(t)

    fields={};

else

    fields=unique( ...
        cellfun(@(x)x{1},t,'UniformOutput',false), ...
        'stable');
end

end


function [found,val] = find_field_recursive_diag(S,target)

found=false;
val=[];

if ~isstruct(S) || numel(S)~=1
    return;
end

if isfield(S,target)

    found=true;
    val=S.(target);
    return;
end

fn=fieldnames(S);

for i=1:numel(fn)

    v=S.(fn{i});

    if isstruct(v) && numel(v)==1

        [found,val]= ...
            find_field_recursive_diag(v,target);

        if found
            return;
        end
    end
end

end