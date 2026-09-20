function R = run_PaviaU_CONTROLLED_FINAL_V3_2(mode)
%RUN_PAVIAU_CONTROLLED_FINAL_V3_2 Controlled Pavia University benchmark.
%
% Modes:
%   INTERFACE : load Pavia data and check the HCL solver/config interface.
%   PROTOCOL  : generate and validate the 10-by-4 corruption supports.
%   AUDIT     : evaluate HCL on all 40 inputs.
%   FULL      : evaluate tSVD-LRA, TNN-TRPCA, p-TRPCA, and HCL.
%
% Reference: rows 178:433, columns 43:298 of PaviaU; divide by 8000
% and clip to [0,1]. The working tensor is 256-by-256-by-103.
% Rank is prespecified as 25; the optional rank helper checks the stored
% 98% Fourier-energy rule with cap 25 when it is available.
%
% Supports are generated in slice, block, entry order and are disjoint.
% Scenarios use 5% entry, 10% block, 10% slice, or their mixture.
% Block sizes are [8 12 16 24], with tensor-wide support fractions.
% Corruption has independent random signs and magnitude Uniform[6,10]
% times the standard deviation of the reference.
%
% Metrics use raw, unclipped reconstructions: NRE, mean bandwise PSNR
% and SSIM with data range 1, and mean spectral angle in degrees.
% Detection uses final po >= 0.5; mask metrics use the frozen mask.
% LevelMacroF1 denotes granularity Macro-F1 over present true classes.
%
% Scientific settings remain frozen, including epsWarm=epsSel=1e-2,
% epsRec=1e-3, recovery budget 80, and relative-step tolerance 1e-5.
% FULL requires both external baseline implementations.
% Outputs are placed in a new timestamped directory.

if nargin < 1 || isempty(mode)
    mode = 'PROTOCOL';
end
mode = upper(char(string(mode)));
assert(ismember(mode, {'INTERFACE','PROTOCOL','AUDIT','FULL'}), ...
    'Mode must be INTERFACE, PROTOCOL, AUDIT, or FULL.');

%% ------------------------------------------------------------------------
% 0. Locate project/config/data
% -------------------------------------------------------------------------
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));

assert(exist('HCL_DSP_FINAL_config','file')==2, ...
    'HCL_DSP_FINAL_config.m is not on the MATLAB path.');
C = HCL_DSP_FINAL_config();

paviaFile = resolve_pavia_file(C, projectRoot);
fprintf('\n============================================================\n');
fprintf(' PAVIAU CONTROLLED FINAL V3.2 | %s\n', mode);
fprintf('============================================================\n');
fprintf('Project root : %s\n', projectRoot);
fprintf('Pavia file   : %s\n', paviaFile);

%% ------------------------------------------------------------------------
% 1. Frozen V3 protocol -- do not alter after looking at V3 outcomes
% -------------------------------------------------------------------------
P = struct();
P.version             = 'PAVIA_CONTROLLED_FINAL_V3_2';
P.seeds               = 20280001:20280010;
P.rank                = 25;
P.rankEnergy          = 0.98;
P.rankCap             = 25;
P.rankStride          = 4;
P.roiRows             = 178:433;
P.roiCols             = 43:298;
P.normalizationMax    = 8000;
P.blockSizes          = [8 12 16 24];
P.ampSigmaMin         = 6;
P.ampSigmaMax         = 10;
P.detectThreshold     = 0.50;
P.epsRecExpected      = 1e-3;
P.maxUpdatesExpected  = 80;
P.relTolExpected      = 1e-5;
P.figureSeed          = 20280001;
P.figureScenario      = 'Mixed';

P.scenarios = [ ...
    mk_scenario('Entry', 0.05, 0.00, 0.00), ...
    mk_scenario('Block', 0.00, 0.10, 0.00), ...
    mk_scenario('Slice', 0.00, 0.00, 0.10), ...
    mk_scenario('Mixed', 0.05, 0.10, 0.10) ...
    ];

% Use a NEW timestamped output directory so V2 and V3 are never pooled.
stamp = datestr(now,'yyyymmddTHHMMSS');
outBase = fullfile(projectRoot,'results_HCL_DSP_FINAL');
if exist(outBase,'dir')~=7, mkdir(outBase); end
outRoot = fullfile(outBase, sprintf('PAVIA_CONTROLLED_V3_2_%s_%s',mode,stamp));
mkdir(outRoot);
figRoot = fullfile(outRoot,'figures');
mkdir(figRoot);

fprintf('Output root  : %s\n', outRoot);
fprintf('Seeds        : %d ... %d (%d)\n',P.seeds(1),P.seeds(end),numel(P.seeds));
fprintf('Scenarios    : Entry5 | Block10 | Slice10 | Mixed E5+B10+S10\n');
fprintf('Block sizes  : %s\n', mat2str(P.blockSizes));
fprintf('Amplitude    : Uniform[%g,%g] x std(reference), random signs\n', ...
    P.ampSigmaMin,P.ampSigmaMax);
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% 2. Load exact PaviaU reference and freeze rank
% -------------------------------------------------------------------------
[Xfull, sourceVar] = load_largest_3d_numeric(paviaFile);

assert(ndims(Xfull)==3, 'PaviaU source must be a numeric 3-D tensor.');
assert(size(Xfull,1)>=max(P.roiRows) && ...
       size(Xfull,2)>=max(P.roiCols) && ...
       size(Xfull,3)==103, ...
       'Expected a PaviaU-like cube with at least 433x298x103 entries.');

rawMin = min(Xfull(:));
rawMax = max(Xfull(:));

Xfull = double(Xfull) ./ P.normalizationMax;
Xfull = min(max(Xfull,0),1);   % reference normalization only
Xref = Xfull(P.roiRows,P.roiCols,:);
clear Xfull

assert(isequal(size(Xref),[256 256 103]), ...
    'Working PaviaU tensor must be exactly 256x256x103.');
assert(all(isfinite(Xref(:))), 'Reference contains non-finite values.');

refSigma = std(Xref(:));
assert(refSigma > 0, 'PaviaU reference has zero variance.');

fprintf('Source variable : %s\n',sourceVar);
fprintf('Raw range       : [%.6g, %.6g]\n',rawMin,rawMax);
fprintf('Working tensor  : %s\n',mat2str(size(Xref)));
fprintf('Reference sigma : %.8g\n',refSigma);
fprintf('Amplitude range : [%.8g, %.8g]\n', ...
    P.ampSigmaMin*refSigma,P.ampSigmaMax*refSigma);

% Verify the previously frozen rank rule if the helper is available.
if exist('hcl_estimate_tubal_rank','file')==2
    rCheck = hcl_estimate_tubal_rank( ...
        Xref,P.rankEnergy,P.rankCap,P.rankStride);
    fprintf('98%% energy rank check (cap %d) = %d\n',P.rankCap,rCheck);
    assert(rCheck==P.rank, ...
        ['Frozen Pavia rank mismatch: expected r=25 but helper returned %d. ' ...
         'Do not silently continue; reconcile preprocessing/rank helper first.'],rCheck);
else
    warning('hcl_estimate_tubal_rank not found; using prespecified frozen r=25.');
end

%% ------------------------------------------------------------------------
% 3. Verify frozen method dependencies before any expensive run
% -------------------------------------------------------------------------
if ismember(mode,{'INTERFACE','AUDIT','FULL'})
    assert(exist('hcl_trpca_final','file')==2, ...
        ['hcl_trpca_final.m is required for the paper-final HCL run. ' ...
         'Do NOT silently fall back to the old single-floor V0.4 core, ' ...
         'because it cannot enforce epsWarm=epsSel=1e-2 and epsRec=1e-3 separately.']);
end

if strcmp(mode,'FULL')
    assert(exist('my_TNN_adapter','file')==2, ...
        'FULL mode requires my_TNN_adapter.m.');
    assert(exist('my_pTRPCA_adapter','file')==2, ...
        'FULL mode requires my_pTRPCA_adapter.m.');
    assert_pTRPCA_frozen(C);
end

%% ------------------------------------------------------------------------
% 3.5 HCL interface/config preflight BEFORE the 40-task audit
% -------------------------------------------------------------------------
if ismember(mode,{'INTERFACE','AUDIT','FULL'})
    interfaceReport = preflight_hcl_interface(C,P);
else
    interfaceReport = struct();
end

if strcmp(mode,'INTERFACE')
    R = struct();
    R.mode = mode;
    R.interface = interfaceReport;
    R.config = P;
    R.outRoot = outRoot;
    save(fullfile(outRoot,'Pavia_V3_2_interface_report.mat'), ...
        'R','P','interfaceReport','-v7');
    fprintf('\nINTERFACE SELF-TEST PASSED.\n');
    fprintf('Now run:\n');
    fprintf('  R1 = run_PaviaU_CONTROLLED_FINAL_V3_2(''AUDIT'');\n');
    return;
end

%% ------------------------------------------------------------------------
% 4. Generate all 40 controlled inputs and audit the protocol first
% -------------------------------------------------------------------------
protocolRows = cell(0,18);
protocolCache = cell(numel(P.seeds),numel(P.scenarios));

for iseed = 1:numel(P.seeds)
    seed = P.seeds(iseed);

    for isc = 1:numel(P.scenarios)
        sc = P.scenarios(isc);

        % Same seed/scenario mapping in PROTOCOL, AUDIT, and FULL.
        generatorSeed = seed + 1000*isc;
        rng(generatorSeed,'twister');

        [Y,truth,proto] = generate_corruption_v3(Xref,sc,P,refSigma);

        protocolCache{iseed,isc} = struct( ...
            'generatorSeed',generatorSeed, ...
            'truth',truth, ...
            'proto',proto);

        protocolRows(end+1,:) = { ... %#ok<AGROW>
            sc.name, seed, generatorSeed, ...
            sc.entryRate,sc.blockRate,sc.sliceRate, ...
            proto.entryRate,proto.blockRate,proto.sliceRate,proto.totalRate, ...
            proto.nEntry,proto.nBlock,proto.nSlice,proto.nTotal, ...
            proto.numBlocks,proto.numSlices, ...
            proto.protocolPass,proto.disjointPass};

        fprintf(['Protocol | seed=%d | %-5s | E %.3f%% | B %.3f%% | ' ...
                 'S %.3f%% | total %.3f%% | blocks=%d | PASS=%d\n'], ...
            seed,sc.name, ...
            100*proto.entryRate,100*proto.blockRate,100*proto.sliceRate, ...
            100*proto.totalRate,proto.numBlocks,proto.protocolPass);

        clear Y
    end
end

protocolNames = { ...
    'Scenario','Seed','GeneratorSeed', ...
    'TargetEntryRate','TargetBlockRate','TargetSliceRate', ...
    'ActualEntryRate','ActualBlockRate','ActualSliceRate','ActualTotalRate', ...
    'NumEntry','NumBlock','NumSlice','NumTotal', ...
    'NumBlocks','NumSlices','ProtocolPass','DisjointPass'};

Tprotocol = cell2table(protocolRows,'VariableNames',protocolNames);
writetable(Tprotocol,fullfile(outRoot,'Pavia_V3_protocol_raw.csv'));

TprotocolSummary = summarize_protocol(Tprotocol,P.scenarios);
writetable(TprotocolSummary,fullfile(outRoot,'Pavia_V3_protocol_summary.csv'));

assert(all(Tprotocol.ProtocolPass), ...
    'V3 corruption protocol audit failed. Inspect Pavia_V3_protocol_raw.csv.');
assert(all(Tprotocol.DisjointPass), ...
    'V3 supports are not mutually exclusive.');

fprintf('\n================ PROTOCOL AUDIT PASSED ================\n');
disp(TprotocolSummary);

if strcmp(mode,'PROTOCOL')
    R = struct();
    R.mode = mode;
    R.protocol = Tprotocol;
    R.protocolSummary = TprotocolSummary;
    R.config = P;
    R.outRoot = outRoot;

    save(fullfile(outRoot,'Pavia_V3_results_small.mat'), ...
        'R','P','Tprotocol','TprotocolSummary','-v7');

    fprintf('\nPROTOCOL mode finished. No recovery was run.\n');
    fprintf('Next recommended command:\n');
    fprintf('  R1 = run_PaviaU_CONTROLLED_FINAL_V3_2(''AUDIT'');\n');
    return;
end

%% ------------------------------------------------------------------------
% 5. Recovery methods
% -------------------------------------------------------------------------
if strcmp(mode,'AUDIT')
    methods = {'HCL-TRPCA'};
else
    methods = {'tSVD-LRA','TNN-TRPCA','p-TRPCA','HCL-TRPCA'};
end

resultRows = cell(0,33);
representativeSaved = false;

for iseed = 1:numel(P.seeds)
    seed = P.seeds(iseed);

    for isc = 1:numel(P.scenarios)
        sc = P.scenarios(isc);

        % Regenerate exactly; do not store 40 full cubes in memory.
        generatorSeed = seed + 1000*isc;
        rng(generatorSeed,'twister');
        [Y,truth,proto] = generate_corruption_v3(Xref,sc,P,refSigma);

        fprintf('\n============================================================\n');
        fprintf('Pavia V3 | seed=%d | %s | total contamination %.3f%%\n', ...
            seed,sc.name,100*proto.totalRate);
        fprintf('============================================================\n');

        for im = 1:numel(methods)
            method = methods{im};
            fprintf('  %-12s : ',method);

            t0 = tic;
            try
                switch method
                    case 'tSVD-LRA'
                        Xhat = run_tsvd_lra(Y,P.rank);
                        H = [];

                    case 'TNN-TRPCA'
                        O = my_TNN_adapter(Y,Xref,truth,C);
                        Xhat = require_xhat(O,method);
                        H = [];
                        clear O

                    case 'p-TRPCA'
                        O = my_pTRPCA_adapter(Y,Xref,truth,C);
                        Xhat = require_xhat(O,method);
                        H = [];
                        clear O

                    case 'HCL-TRPCA'
                        H = run_hcl_final_compat(Y,P.rank,C);
                        Xhat = require_xhat(H,method);
                        verify_hcl_final_floor(H,P.epsRecExpected);

                    otherwise
                        error('Unknown method %s',method);
                end

                runtime = toc(t0);
                assert(isequal(size(Xhat),size(Xref)), ...
                    '%s returned the wrong tensor size.',method);
                assert(all(isfinite(Xhat(:))), ...
                    '%s returned non-finite reconstruction values.',method);

                met = recovery_metrics_raw(Xref,Xhat);

                % HCL-only structural diagnostics.
                detP=NaN; detR=NaN; detF=NaN;
                maskP=NaN; maskR=NaN; maskF=NaN;
                levelMacro=NaN; levelE=NaN; levelB=NaN; levelS=NaN;
                selectedFraction=NaN; rESS=NaN;
                iterations=NaN; finalRelChange=NaN; toleranceMet=NaN;

                if strcmp(method,'HCL-TRPCA')
                    D = extract_hcl_diagnostics(H,size(Xref),P);
                    [detP,detR,detF] = binary_metrics(D.detectMask,truth.anyMask);
                    [maskP,maskR,maskF] = binary_metrics(D.recoveryMask,truth.anyMask);
                    [levelMacro,levelE,levelB,levelS] = ...
                        level_f1(D.predLevel,truth.levelLabel);

                    selectedFraction = D.selectedFraction;
                    rESS = D.rESS;
                    iterations = D.iterations;
                    finalRelChange = D.finalRelChange;
                    toleranceMet = D.toleranceMet;

                    % Fixed, prespecified representative visualization.
                    if ~representativeSaved && ...
                            seed==P.figureSeed && strcmp(sc.name,P.figureScenario)
                        save_representative_figures( ...
                            Xref,Y,Xhat,truth,D,seed,sc.name,figRoot);
                        representativeSaved = true;
                    end
                end

                resultRows(end+1,:) = { ... %#ok<AGROW>
                    sc.name,seed,generatorSeed,method,P.rank, ...
                    sc.entryRate,sc.blockRate,sc.sliceRate, ...
                    proto.entryRate,proto.blockRate,proto.sliceRate,proto.totalRate, ...
                    met.NRE,met.MPSNR,met.MSSIM,met.SAMdeg, ...
                    detP,detR,detF, ...
                    maskP,maskR,maskF, ...
                    levelMacro,levelE,levelB,levelS, ...
                    selectedFraction,rESS,iterations,finalRelChange,toleranceMet, ...
                    met.OutOfRangeFraction,runtime};

                fprintf(['NRE %.6f | MPSNR %.3f | SAM %.3f deg | ' ...
                         'DetF1 %.4f | MaskF1 %.4f | LevelF1 %.4f | %.1fs\n'], ...
                    met.NRE,met.MPSNR,met.SAMdeg,detF,maskF,levelMacro,runtime);

                clear Xhat H D

            catch ME
                runtime = toc(t0);
                fprintf('\nFAILED after %.1fs: %s\n',runtime,ME.message);

                % Save all completed rows before rethrowing.
                Tpartial = build_result_table(resultRows);
                if ~isempty(Tpartial)
                    writetable(Tpartial,fullfile(outRoot,'Pavia_V3_raw_PARTIAL.csv'));
                    save(fullfile(outRoot,'Pavia_V3_partial.mat'), ...
                        'Tpartial','P','Tprotocol','TprotocolSummary','-v7');
                end
                rethrow(ME);
            end
        end

        % Keep memory bounded between large real-HSI runs.
        clear Y truth proto
    end
end

%% ------------------------------------------------------------------------
% 6. Tables and compact output
% -------------------------------------------------------------------------
T = build_result_table(resultRows);
writetable(T,fullfile(outRoot,'Pavia_V3_raw.csv'));

Tsum = summarize_results(T,methods,P.scenarios);
writetable(Tsum,fullfile(outRoot,'Pavia_V3_summary.csv'));

R = struct();
R.mode = mode;
R.raw = T;
R.summary = Tsum;
R.protocol = Tprotocol;
R.protocolSummary = TprotocolSummary;
R.config = P;
R.outRoot = outRoot;
R.source = struct( ...
    'paviaFile',paviaFile, ...
    'sourceVariable',sourceVar, ...
    'rawMin',rawMin, ...
    'rawMax',rawMax, ...
    'referenceSigma',refSigma);

save(fullfile(outRoot,'Pavia_V3_results_small.mat'), ...
    'R','P','T','Tsum','Tprotocol','TprotocolSummary','-v7');

fprintf('\n============================================================\n');
fprintf(' PAVIAU CONTROLLED V3.2 COMPLETE | %s\n',mode);
fprintf('============================================================\n');
disp(Tsum);
fprintf('Results saved to:\n  %s\n',outRoot);

if strcmp(mode,'AUDIT')
    fprintf('\nAUDIT completed with the corrected protocol.\n');
    fprintf(['Inspect Block/Mixed Detection-F1, Mask-F1 and LevelMacro-F1. ' ...
             'Do NOT retune HCL from these results.\n']);
    fprintf('Full baseline comparison command:\n');
    fprintf('  R2 = run_PaviaU_CONTROLLED_FINAL_V3_2(''FULL'');\n');
end
end


%% =========================================================================
% SCENARIO / DATA HELPERS
% =========================================================================
function s = mk_scenario(name,e,b,sl)
s = struct('name',name,'entryRate',e,'blockRate',b,'sliceRate',sl);
end

function paviaFile = resolve_pavia_file(C,projectRoot)
% Public-release resolver: no machine-specific absolute paths.
% Search priority:
%   1) C.real.paviaMat
%   2) C.real.dataRoot/PaviaU.mat
%   3) configs/local_paths.m (local-only, gitignored)
%   4) HCL_TRPCA_DATA_ROOT environment variable
%   5) repository-local data/PaviaU.mat

candidates = {};

if isstruct(C) && isfield(C,'real')
    if isfield(C.real,'paviaMat') && ~isempty(C.real.paviaMat)
        candidates{end+1} = C.real.paviaMat; %#ok<AGROW>
    end
    if isfield(C.real,'dataRoot') && ~isempty(C.real.dataRoot)
        candidates{end+1} = fullfile(C.real.dataRoot,'PaviaU.mat'); %#ok<AGROW>
    end
end

% Optional local, gitignored machine-specific path configuration.
if exist('local_paths','file') == 2
    try
        P = local_paths();
        if isfield(P,'pavia') && ~isempty(P.pavia)
            candidates{end+1} = P.pavia; %#ok<AGROW>
        end
    catch ME
        warning('local_paths.m could not be read: %s',ME.message);
    end
end

% Optional environment-variable convention.
dataRoot = getenv('HCL_TRPCA_DATA_ROOT');
if ~isempty(dataRoot)
    candidates{end+1} = fullfile(dataRoot,'PaviaU.mat'); %#ok<AGROW>
end

% Optional repository-local data directory.
% In the public tree projectRoot normally points to the repository root
% or to code/pavia; both cases are handled conservatively below.
candidateRoots = {projectRoot, fileparts(projectRoot), fileparts(fileparts(projectRoot))};
for ir = 1:numel(candidateRoots)
    if ~isempty(candidateRoots{ir})
        candidates{end+1} = fullfile(candidateRoots{ir},'data','PaviaU.mat'); %#ok<AGROW>
    end
end

paviaFile = '';
for i = 1:numel(candidates)
    f = char(candidates{i});
    if exist(f,'file') == 2
        paviaFile = f;
        break;
    end
end

assert(~isempty(paviaFile), [ ...
    'Could not locate PaviaU.mat. Configure C.real.paviaMat, ' ...
    'configs/local_paths.m, HCL_TRPCA_DATA_ROOT, or repo/data/.']);
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
function H = run_hcl_final_compat(Y,r,C)
%RUN_HCL_FINAL_COMPAT  Exact three-input adapter for the local final solver.
%
% The local error log established the interface:
%       H = hcl_trpca_final(Y,r,cfg)
% The failure was not a solver failure; the project-level struct C was passed
% where the solver-level cfg was required.
%
% V3.2 therefore NEVER passes C directly. It builds an explicit solver cfg,
% verifies every cfg.<field> directly referenced by the local solver source,
% and then calls the confirmed 3-input interface.

Ptmp = struct('maxUpdatesExpected',80,'relTolExpected',1e-5, ...
              'epsRecExpected',1e-3);
cfg = build_frozen_hcl_solver_cfg(C,r,Ptmp);

H = hcl_trpca_final(Y,r,cfg);

Xhat = require_xhat(H,'HCL-TRPCA');
assert(isequal(size(Xhat),size(Y)) && all(isfinite(Xhat(:))), ...
    'hcl_trpca_final returned an invalid reconstruction.');
end

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
           'Check hcl_trpca_final.m and HCL_DSP_FINAL_config.m for ' ...
           'missing or incompatible solver configuration fields.']);
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

function assert_pTRPCA_frozen(C)
assert(isfield(C,'pTRPCA') && isfield(C.pTRPCA,'finalized') ...
    && logical(C.pTRPCA.finalized), ...
    'p-TRPCA parameters must be frozen before FULL Pavia evaluation.');

expectedW = [1 1.1 1.5];
assert(isfield(C.pTRPCA,'wValues') && ...
    max(abs(double(C.pTRPCA.wValues(:)')-expectedW))<1e-12, ...
    'Unexpected p-TRPCA wValues; do not retune on Pavia.');
assert(abs(double(C.pTRPCA.p1)-1)<1e-12 && ...
       abs(double(C.pTRPCA.p2)-0.8)<1e-12, ...
    'Unexpected p-TRPCA p1/p2; do not retune on Pavia.');
end

function Xhat = run_tsvd_lra(Y,r)
if exist('hcl_project_tubal_rank','file')==2
    Xhat = hcl_project_tubal_rank(Y,r);
    return;
end

% Self-contained exact dense truncated t-SVD fallback.
Yf = fft(double(Y),[],3);
Xf = zeros(size(Yf),'like',Yf);

for k=1:size(Yf,3)
    [U,S,V] = svd(Yf(:,:,k),'econ');
    rr = min([r,size(S,1),size(S,2)]);
    Xf(:,:,k) = U(:,1:rr)*S(1:rr,1:rr)*V(:,1:rr)';
end

Xhat = real(ifft(Xf,[],3));
end


%% =========================================================================
% HCL DIAGNOSTIC EXTRACTION
% =========================================================================
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


%% =========================================================================
% TABLES / SUMMARIES
% =========================================================================
function T = build_result_table(rows)
names={ ...
    'Scenario','Seed','GeneratorSeed','Method','Rank', ...
    'TargetEntryRate','TargetBlockRate','TargetSliceRate', ...
    'ActualEntryRate','ActualBlockRate','ActualSliceRate','ActualTotalRate', ...
    'NRE','MPSNR','MSSIM','SAMdeg', ...
    'DetectionPrecision','DetectionRecall','DetectionF1', ...
    'MaskPrecision','MaskRecall','MaskF1', ...
    'LevelMacroF1','LevelF1Entry','LevelF1Block','LevelF1Slice', ...
    'SelectedFraction','rESS','Iterations','FinalRelChange','ToleranceMet', ...
    'OutOfRangeFraction','RuntimeSec'};

if isempty(rows)
    T=cell2table(cell(0,numel(names)),'VariableNames',names);
else
    T=cell2table(rows,'VariableNames',names);
end
end

function S = summarize_protocol(T,scenarios)
rows={};
for i=1:numel(scenarios)
    sc=scenarios(i);
    idx=strcmp(T.Scenario,sc.name);

    rows(end+1,:)={ ... %#ok<AGROW>
        sc.name,nnz(idx), ...
        mean(T.ActualEntryRate(idx)),std(T.ActualEntryRate(idx)), ...
        mean(T.ActualBlockRate(idx)),std(T.ActualBlockRate(idx)), ...
        mean(T.ActualSliceRate(idx)),std(T.ActualSliceRate(idx)), ...
        mean(T.ActualTotalRate(idx)),std(T.ActualTotalRate(idx)), ...
        min(T.ProtocolPass(idx)),min(T.DisjointPass(idx))};
end

S=cell2table(rows,'VariableNames',{ ...
    'Scenario','N', ...
    'EntryRate_Mean','EntryRate_Std', ...
    'BlockRate_Mean','BlockRate_Std', ...
    'SliceRate_Mean','SliceRate_Std', ...
    'TotalRate_Mean','TotalRate_Std', ...
    'AllProtocolPass','AllDisjointPass'});
end

function S = summarize_results(T,methods,scenarios)

metrics={ ...
    'ActualTotalRate', ...
    'NRE','MPSNR','MSSIM','SAMdeg', ...
    'DetectionF1','MaskF1','LevelMacroF1', ...
    'SelectedFraction','rESS','Iterations','FinalRelChange', ...
    'OutOfRangeFraction','RuntimeSec'};

rows={};
for is=1:numel(scenarios)
    sc=scenarios(is);

    for im=1:numel(methods)
        method=methods{im};
        idx=strcmp(T.Scenario,sc.name) & strcmp(T.Method,method);

        row={sc.name,method,nnz(idx)};
        for j=1:numel(metrics)
            x=T.(metrics{j})(idx);
            x=x(isfinite(x));
            if isempty(x)
                mu=NaN; sd=NaN;
            else
                mu=mean(x); sd=std(x);
            end
            row=[row,{mu,sd}]; %#ok<AGROW>
        end
        rows(end+1,:)=row; %#ok<AGROW>
    end
end

names={'Scenario','Method','N'};
for j=1:numel(metrics)
    names=[names,{[metrics{j} '_Mean'],[metrics{j} '_Std']}]; %#ok<AGROW>
end

S=cell2table(rows,'VariableNames',names);
end


%% =========================================================================
% REPRESENTATIVE FIGURES -- FIXED BY RULE, NOT VISUAL CHERRY-PICKING
% =========================================================================
function save_representative_figures(Xref,Y,Xhat,truth,D,seed,scName,figRoot)

% Local non-slice band: deterministic argmax of true entry+block support.
localScore=squeeze(sum(sum(truth.entryMask | truth.blockMask,1),2));
sliceBandFlag=squeeze(any(any(truth.sliceMask,1),2));
localScore(sliceBandFlag)=-Inf;
[bestLocal,kLocal]=max(localScore);
assert(isfinite(bestLocal) && bestLocal>0, ...
    'No local-corruption band found for representative figure.');

% Slice band: deterministic first true slice index.
sliceIds=find(sliceBandFlag);
assert(~isempty(sliceIds),'No slice band found for Mixed representative.');
kSlice=sliceIds(1);

make_one_mechanism_figure( ...
    Xref,Y,Xhat,truth,D,kLocal,seed,scName,'local',figRoot);

make_one_mechanism_figure( ...
    Xref,Y,Xhat,truth,D,kSlice,seed,scName,'slice',figRoot);

% Save only compact 2-D representative evidence, never the whole 3-D run.
Rep=struct();
Rep.seed=seed;
Rep.scenario=scName;
Rep.localBand=kLocal;
Rep.sliceBand=kSlice;
Rep.local=compact_band(Xref,Y,Xhat,truth,D,kLocal);
Rep.slice=compact_band(Xref,Y,Xhat,truth,D,kSlice);
save(fullfile(figRoot,sprintf('%s_seed%d_representative_small.mat',scName,seed)), ...
    'Rep','-v7');
end

function S=compact_band(Xref,Y,Xhat,truth,D,k)
S=struct();
S.band=k;
S.reference=Xref(:,:,k);
S.observation=Y(:,:,k);
S.recovery=Xhat(:,:,k);
S.residual=abs(Y(:,:,k)-Xhat(:,:,k));
S.trueLevel=truth.levelLabel(:,:,k);
S.predLevel=D.predLevel(:,:,k);
S.recoveryMask=D.recoveryMask(:,:,k);
S.piEntry=D.piEntry(:,:,k);
S.piBlock=D.piBlock(:,:,k);
S.piSlice=D.piSlice(:,:,k);
end

function make_one_mechanism_figure(Xref,Y,Xhat,truth,D,k,seed,scName,tag,figRoot)

f=figure('Color','w','Visible','off', ...
    'Position',[50 50 1900 720], ...
    'Name',sprintf('Pavia V3 %s %s band %d',scName,tag,k));

subplot(2,5,1);
imagesc(Xref(:,:,k),[0 1]); axis image off; colorbar; title('(a) Reference');

subplot(2,5,2);
imagesc(min(max(Y(:,:,k),0),1),[0 1]); axis image off; colorbar; title('(b) Corrupted');

subplot(2,5,3);
imagesc(min(max(Xhat(:,:,k),0),1),[0 1]); axis image off; colorbar; title('(c) HCL recovery');

subplot(2,5,4);
imagesc(double(truth.levelLabel(:,:,k)),[0 3]); axis image off; colorbar; title('(d) True level');

subplot(2,5,5);
imagesc(double(D.predLevel(:,:,k)),[0 3]); axis image off; colorbar; title('(e) Predicted level');

subplot(2,5,6);
imagesc(double(D.recoveryMask(:,:,k)),[0 1]); axis image off; colorbar; title('(f) Frozen mask');

subplot(2,5,7);
imagesc(D.piEntry(:,:,k),[0 1]); axis image off; colorbar; title('(g) \pi_e');

subplot(2,5,8);
imagesc(D.piBlock(:,:,k),[0 1]); axis image off; colorbar; title('(h) \pi_b');

subplot(2,5,9);
imagesc(D.piSlice(:,:,k),[0 1]); axis image off; colorbar; title('(i) \pi_s');

subplot(2,5,10);
imagesc(abs(Y(:,:,k)-Xhat(:,:,k))); axis image off; colorbar; title('(j) |Residual|');

sgtitle(sprintf('PaviaU V3 | %s | seed %d | band %d (%s)', ...
    scName,seed,k,tag),'Interpreter','none');

base=fullfile(figRoot,sprintf('%s_seed%d_%s_band%03d',scName,seed,tag,k));
savefig(f,[base '.fig']);

try
    exportgraphics(f,[base '.png'],'Resolution',300);
    exportgraphics(f,[base '.pdf'],'ContentType','vector');
catch
    print(f,[base '.png'],'-dpng','-r300');
end

close(f);
end
