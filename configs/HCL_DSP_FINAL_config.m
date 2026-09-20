function C = HCL_DSP_FINAL_config()
%HCL_DSP_FINAL_CONFIG Public-release final configuration.
%
% Scientific parameters are frozen to the final manuscript settings.
% Machine-specific paths are loaded only from configs/local_paths.m,
% which must remain untracked.

C = struct();
C.version = 'HCL_DSP_FINAL_1_0_PUBLIC';

cfgDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(cfgDir);

C.projectRoot = repoRoot;
C.outputRoot  = fullfile(repoRoot,'generated_results');

%% ------------------------ Frozen HCL parameters -------------------------
C.hcl = struct();

C.hcl.epsWarm = 1e-2;
C.hcl.epsSel  = 1e-2;
C.hcl.epsRec  = 1e-3;

C.hcl.warmIterations = 10;
C.hcl.maxIterations  = 80;

C.hcl.b = 8;
C.hcl.entryCenter = 2.5;
C.hcl.entryTemp   = 0.45;
C.hcl.blockCenter = 0.45;
C.hcl.blockTemp   = 0.06;
C.hcl.sliceCenter = 0.50;
C.hcl.sliceTemp   = 0.06;

C.hcl.betaBlock = 0.65;
C.hcl.betaSlice = 0.75;
C.hcl.nu        = 0.50;

C.hcl.tauDetect = 0.50;
C.hcl.tauHigh   = 0.90;
C.hcl.fMax      = 0.32;
C.hcl.etaMin    = 0.68;

C.hcl.gamma  = 0.85;
C.hcl.relTol = 1e-5;

C.hcl.tukeyEntry = 4.685;
C.hcl.tukeyBlock = 4.0;
C.hcl.tukeySlice = 3.5;
C.hcl.scaleEps   = 1e-8;

%% ---------------------- Synthetic confirmation --------------------------
C.synthetic = struct();
C.synthetic.I = 64;
C.synthetic.J = 64;
C.synthetic.K = 30;
C.synthetic.trueRank = 5;
C.synthetic.snrDb = 35;
C.synthetic.ampRange = [6 10];
C.synthetic.blockSizes = [6 10 14];

C.synthetic.scenarios.EB  = struct('entry',0.05,'block',0.10,'slice',0.00);
C.synthetic.scenarios.ES  = struct('entry',0.05,'block',0.00,'slice',0.10);
C.synthetic.scenarios.BS  = struct('entry',0.00,'block',0.10,'slice',0.10);
C.synthetic.scenarios.EBS = struct('entry',0.05,'block',0.10,'slice',0.10);

C.confirmationSeeds = 20270001:20270030;
C.developmentSeeds  = 20261001:20261005;
C.robustnessSeeds   = 20271001:20271020;
C.failureSeeds      = 20272001:20272010;

C.weakAmpRange = [3 5];
C.rankGrid = [3 4 5 6 8];

%% --------------------------- Experiment flags ---------------------------
C.flags = struct();
C.flags.saveFullSyntheticState = true;
C.flags.savePerRunMat = true;
C.flags.overwrite = false;

%% -------------------------- p-TRPCA configuration -----------------------
C.pTRPCA = struct();
C.pTRPCA.finalized = true;
C.pTRPCA.wValues = [1 1.1 1.5];
C.pTRPCA.p1 = 1.0;
C.pTRPCA.p2 = 0.8;
C.pTRPCA.mu = 1e-4;
C.pTRPCA.tol = 1e-8;
C.pTRPCA.rho = 1.1;
C.pTRPCA.maxIter = 500;

%% -------------------------- External baselines --------------------------
C.baselines = struct();
C.baselines.TNN = @my_TNN_adapter;
C.baselines.pTRPCA = @my_pTRPCA_adapter;
C.baselines.EntryRobust = [];

%% -------------------------- Real-data paths ------------------------------
C.real = struct();
C.real.dataRoot = '';
C.real.paviaMat = '';
C.real.urbanMat = '';
C.real.sbiRoot  = '';

localFile = fullfile(cfgDir,'local_paths.m');

if exist(localFile,'file') == 2
    addpath(cfgDir);
    P = local_paths();

    if isfield(P,'dataRoot') && ~isempty(P.dataRoot)
        C.real.dataRoot = P.dataRoot;
    end
    if isfield(P,'pavia') && ~isempty(P.pavia)
        C.real.paviaMat = P.pavia;
    end
    if isfield(P,'urban') && ~isempty(P.urban)
        C.real.urbanMat = P.urban;
    end
    if isfield(P,'sbiRoot') && ~isempty(P.sbiRoot)
        C.real.sbiRoot = P.sbiRoot;
    end
end

% Optional environment-variable fallback.
if isempty(C.real.dataRoot)
    envRoot = getenv('HCL_TRPCA_DATA_ROOT');
    if ~isempty(envRoot)
        C.real.dataRoot = envRoot;
    end
end

if isempty(C.real.paviaMat) && ~isempty(C.real.dataRoot)
    C.real.paviaMat = fullfile(C.real.dataRoot,'PaviaU.mat');
end

if isempty(C.real.urbanMat) && ~isempty(C.real.dataRoot)
    C.real.urbanMat = fullfile(C.real.dataRoot,'Urban_F210.mat');
end

if isempty(C.real.sbiRoot) && ~isempty(C.real.dataRoot)
    candidateMat = fullfile(C.real.dataRoot,'MAT');
    if exist(candidateMat,'dir') == 7
        C.real.sbiRoot = candidateMat;
    else
        C.real.sbiRoot = C.real.dataRoot;
    end
end

C.real.sbiNames = { ...
    'SBI_Board', ...
    'SBI_CAVIAR1', ...
    'SBI_Hall_Monitor', ...
    'SBI_People_Foliage'};

end
