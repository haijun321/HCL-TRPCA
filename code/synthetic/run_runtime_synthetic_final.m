function R = run_runtime_synthetic_final(mode)
% RUN_RUNTIME_SYNTHETIC_FINAL
% Paper-facing runtime benchmark for the FROZEN HCL protocol.
%
% Usage:
%   R = run_runtime_synthetic_final('PRECHECK');
%   R = run_runtime_synthetic_final('FINAL');
%
% PRECHECK:
%   May run with only the currently registered methods.
%   Useful for checking code, but NOT the final paper runtime table.
%
% FINAL:
%   Requires exact registered TNN-TRPCA and p-TRPCA adapters.
%   Refuses to run if more than one MATLAB process is detected on Windows
%   when that check is available.
%
% FINAL RUNTIME RULES:
%   - Do not run PaviaU, SBI, Urban, Rank, Weak/Overlap, or Stability
%     simultaneously.
%   - Use the same machine and MATLAB session conditions for every method.
%   - Do not modify the frozen HCL configuration.
%   - Use exact cited external baseline implementations, not substitutes.

if nargin<1, mode='PRECHECK'; end
mode = upper(string(mode));

C = HCL_DSP_FINAL_config();

if mode=="FINAL"
    assert(~isempty(C.baselines.TNN), ...
        ['FINAL runtime requires C.baselines.TNN to be registered ', ...
         'to the exact verified TNN-TRPCA implementation.']);
    assert(~isempty(C.baselines.pTRPCA), ...
        ['FINAL runtime requires C.baselines.pTRPCA to be registered ', ...
         'to the exact verified p-TRPCA implementation.']);

    nMatlab = count_matlab_processes();
    if isfinite(nMatlab) && nMatlab > 1
        error(['FINAL runtime aborted: %d MATLAB processes detected. ', ...
               'Close the other MATLAB windows and rerun runtime alone.'], ...
               nMatlab);
    end
end

outDir = fullfile(C.outputRoot,'RUNTIME_FINAL');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Same paired EBS instances for all methods.
seeds = C.confirmationSeeds(1:min(5,numel(C.confirmationSeeds)));

% Three repeated timings per seed; the per-seed MEDIAN is used for the
% paper-facing summary to reduce incidental scheduler jitter.
nRepeats = 3;

methods = build_methods(C);

% Warm-up on the first paired instance, untimed.
D0 = hcl_synthetic_case(seeds(1),'EBS',C,'DISJOINT');
fprintf('\n=== RUNTIME WARM-UP (not timed) ===\n');
for im = 1:numel(methods)
    fprintf('Warm-up: %s\n',methods(im).name);
    rng(12345+im,'twister');
    Xhat = methods(im).runner(D0); %#ok<NASGU>
end

rows = {};
rid = 0;

fprintf('\n=== RUNTIME BENCHMARK: %s ===\n',mode);

for iseed = 1:numel(seeds)
    seed = seeds(iseed);
    D = hcl_synthetic_case(seed,'EBS',C,'DISJOINT');

    % Rotate method order across seeds to reduce systematic order bias.
    nM = numel(methods);
    order = circshift(1:nM,[0,mod(iseed-1,nM)]);

    for io = 1:numel(order)
        im = order(io);
        name = methods(im).name;

        times = zeros(nRepeats,1);

        for rep = 1:nRepeats
            % Reproducible method-specific RNG state for methods that use RNG.
            rr = mod(double(seed) + 10000*im + 100*rep,2^32-1);
            rng(rr,'twister');

            drawnow;
            tic;
            Xhat = methods(im).runner(D); %#ok<NASGU>
            times(rep) = toc;
        end

        for rep = 1:nRepeats
            rid = rid+1;
            rows(rid,:) = { ...
                seed, char(name), rep, times(rep), ...
                median(times), mean(times), std(times)};
        end

        fprintf('seed=%d | %-16s | median %.4f s\n', ...
            seed,name,median(times));
    end
end

Trep = cell2table(rows,'VariableNames', { ...
    'Seed','Method','Repeat','RuntimeSec', ...
    'SeedMedianRuntimeSec','SeedMeanRuntimeSec','SeedSDRuntimeSec'});

writetable(Trep,fullfile(outDir,'runtime_final_repeats.csv'));

% One row per method x seed
[Tseed,ia] = unique(Trep(:,{'Seed','Method','SeedMedianRuntimeSec', ...
    'SeedMeanRuntimeSec','SeedSDRuntimeSec'}),'rows','stable'); %#ok<ASGLU>
writetable(Tseed,fullfile(outDir,'runtime_final_seed_summary.csv'));

% Paper-facing summary: mean +/- SD of per-seed medians.
methodsU = unique(Tseed.Method,'stable');
sumRows = {};
sid = 0;

for i=1:numel(methodsU)
    idx = strcmp(Tseed.Method,methodsU{i});
    x = Tseed.SeedMedianRuntimeSec(idx);

    sid=sid+1;
    sumRows(sid,:) = { ...
        methodsU{i},numel(x),mean(x),std(x),median(x),min(x),max(x)};
end

S = cell2table(sumRows,'VariableNames', { ...
    'Method','NumSeeds','MeanOfSeedMediansSec','SDOfSeedMediansSec', ...
    'MedianOfSeedMediansSec','MinSeedMedianSec','MaxSeedMedianSec'});

writetable(S,fullfile(outDir,'runtime_final_summary.csv'));

Meta = struct();
Meta.mode = char(mode);
Meta.matlabVersion = version;
Meta.computer = computer;
Meta.numDetectedMatlabProcesses = count_matlab_processes();
Meta.numSeeds = numel(seeds);
Meta.repeatsPerSeed = nRepeats;
Meta.epsRec = C.hcl.epsRec;
Meta.gamma = C.hcl.gamma;
Meta.relTol = C.hcl.relTol;
Meta.formalReady = ...
    mode=="FINAL" && ...
    ~isempty(C.baselines.TNN) && ...
    ~isempty(C.baselines.pTRPCA);

save(fullfile(outDir,'runtime_final_results.mat'),'Trep','Tseed','S','Meta');

fprintf('\nRuntime summary:\n');
disp(S);

if Meta.formalReady
    fprintf('STATUS: FORMAL RUNTIME TABLE READY (subject to isolated machine conditions).\n');
else
    fprintf(['STATUS: PRECHECK ONLY. Register exact TNN and p-TRPCA adapters ', ...
             'and rerun in FINAL mode with no concurrent heavy jobs.\n']);
end

R = struct('repeats',Trep,'seedSummary',Tseed,'summary',S,'meta',Meta);

end

% ========================================================================
function methods = build_methods(C)
methods = struct('name',{},'runner',{});

methods(end+1).name = "tSVD-LRA";
methods(end).runner = @(D) hcl_project_tubal_rank(D.Y,D.rank);

methods(end+1).name = "HCL-TRPCA";
methods(end).runner = @(D) run_hcl_xhat(D,C);

if ~isempty(C.baselines.TNN)
    methods(end+1).name = "TNN-TRPCA";
    methods(end).runner = @(D) run_adapter_xhat(C.baselines.TNN,D,C);
end

if ~isempty(C.baselines.pTRPCA)
    methods(end+1).name = "p-TRPCA";
    methods(end).runner = @(D) run_adapter_xhat(C.baselines.pTRPCA,D,C);
end

if ~isempty(C.baselines.EntryRobust)
    methods(end+1).name = "Entry-Robust";
    methods(end).runner = @(D) run_adapter_xhat(C.baselines.EntryRobust,D,C);
end
end

function Xhat = run_hcl_xhat(D,C)
H = hcl_trpca_final(D.Y,D.rank,C.hcl);
Xhat = H.Xhat;
end

function Xhat = run_adapter_xhat(f,D,C)
O = f(D.Y,D.Xstar,D.truth,C);
assert(isstruct(O) && isfield(O,'Xhat'), ...
    'Baseline adapter must return a struct containing out.Xhat.');
Xhat = O.Xhat;
end

function n = count_matlab_processes()
% Best-effort check. NaN means unavailable.
n = NaN;

try
    if ispc
        [status,out] = system('tasklist /FI "IMAGENAME eq MATLAB.exe" /NH');
        if status==0
            hits = regexpi(out,'MATLAB\.exe','match');
            n = numel(hits);
        end
    elseif isunix
        [status,out] = system('pgrep -fc matlab');
        if status==0
            n = str2double(strtrim(out));
        end
    end
catch
    n = NaN;
end
end
