function D = hcl_synthetic_case(seed, scenarioName, C, mode)
% HCL_SYNTHETIC_CASE
% Synthetic tensor generator matching the manuscript moderate protocol.
%
% mode:
%   'DISJOINT' : slice -> non-overlapping blocks -> isolated entries
%   'WEAK'     : same geometry, weaker amplitude
%   'OVERLAP'  : allows entry/block/slice supports to overlap
%
% D.truth.level labels for disjoint/weak:
%   0 clean, 1 entry, 2 block, 3 slice
% For overlap, use D.truth.activeLevels instead of a unique label.

if nargin < 4, mode = 'DISJOINT'; end
rng(seed,'twister');

S = C.synthetic;
spec = S.scenarios.(scenarioName);

I=S.I; J=S.J; K=S.K; r=S.trueRank; M=I*J*K;

%% Clean tubal-rank-r tensor
A = randn(I,r,K);
B = randn(r,J,K);
Xstar = tprod3(A,B);
sx = std(Xstar(:),0);
if sx <= 0, error('Degenerate synthetic draw.'); end
Xstar = Xstar/sx;

%% Dense Gaussian noise at requested SNR
sigma2 = norm(Xstar(:))^2 / (M * 10^(S.snrDb/10));
N = sqrt(sigma2)*randn(I,J,K);

%% Corruption supports
entryMask = false(I,J,K);
blockMask = false(I,J,K);
sliceMask = false(I,J,K);

switch upper(mode)
    case {'DISJOINT','WEAK'}
        % Whole slices first
        nSlices = round(spec.slice*K);
        if nSlices > 0
            ids = randperm(K,nSlices);
            sliceMask(:,:,ids) = true;
        end

        % Blocks in currently clean support
        targetBlock = round(spec.block*M);
        blockMask = place_blocks(blockMask | sliceMask, targetBlock, ...
                                 S.blockSizes, false);

        % Isolated entries in remaining clean positions
        targetEntry = round(spec.entry*M);
        available = find(~sliceMask & ~blockMask);
        nTake = min(targetEntry,numel(available));
        if nTake > 0
            pick = available(randperm(numel(available),nTake));
            entryMask(pick)=true;
        end

    case 'OVERLAP'
        % Slice support
        nSlices = round(spec.slice*K);
        if nSlices > 0
            ids = randperm(K,nSlices);
            sliceMask(:,:,ids)=true;
        end

        % Blocks may be placed anywhere, including slice-supported regions.
        targetBlock = round(spec.block*M);
        blockMask = place_blocks(false(I,J,K), targetBlock, ...
                                 S.blockSizes, true);

        % Entries may occur anywhere.
        targetEntry = round(spec.entry*M);
        nTake = min(targetEntry,M);
        pick = randperm(M,nTake);
        entryMask(pick)=true;

    otherwise
        error('Unknown mode: %s',mode);
end

%% Additive corruption
if strcmpi(mode,'WEAK')
    amp = C.weakAmpRange;
else
    amp = S.ampRange;
end

E = zeros(I,J,K);

E = E + signed_uniform(entryMask,amp);
E = E + signed_uniform(blockMask,amp);
E = E + signed_uniform(sliceMask,amp);

Y = Xstar + N + E;

%% Truth
unionMask = entryMask | blockMask | sliceMask;
truth = struct();
truth.entryMask = entryMask;
truth.blockMask = blockMask;
truth.sliceMask = sliceMask;
truth.unionMask = unionMask;
truth.actualContamination = nnz(unionMask)/M;
truth.activeLevels = cat(4,entryMask,blockMask,sliceMask);

if strcmpi(mode,'OVERLAP')
    truth.level = [];
else
    level = zeros(I,J,K,'uint8');
    level(entryMask)=1;
    level(blockMask)=2;
    level(sliceMask)=3;
    truth.level=level;
end

D = struct();
D.seed = seed;
D.scenario = scenarioName;
D.mode = upper(mode);
D.Xstar = Xstar;
D.N = N;
D.E = E;
D.Y = Y;
D.truth = truth;
D.rank = r;
end

% ========================================================================
function Bmask = place_blocks(forbidden, targetCount, blockSizes, allowOverlap)
[I,J,K] = size(forbidden);
Bmask = false(I,J,K);
maxAttempts = max(10000,20*ceil(targetCount/max(blockSizes)^2));
attempt = 0;

while nnz(Bmask) < targetCount && attempt < maxAttempts
    attempt = attempt + 1;
    side = blockSizes(randi(numel(blockSizes)));
    if side > I || side > J, continue; end
    k = randi(K);
    i0 = randi(I-side+1);
    j0 = randi(J-side+1);
    rr = i0:(i0+side-1);
    cc = j0:(j0+side-1);

    regionForbidden = forbidden(rr,cc,k);
    regionExisting = Bmask(rr,cc,k);

    if allowOverlap
        % overlap with other mechanisms allowed; avoid duplicating same block mask too much
        if mean(regionExisting(:)) > 0.25, continue; end
    else
        if any(regionForbidden(:)) || any(regionExisting(:)), continue; end
    end

    Bmask(rr,cc,k)=true;
end
end

function E = signed_uniform(mask,amp)
E = zeros(size(mask));
idx = find(mask);
if isempty(idx), return; end
sgn = 2*(rand(numel(idx),1)>0.5)-1;
mag = amp(1) + (amp(2)-amp(1))*rand(numel(idx),1);
E(idx) = sgn.*mag;
end

function C = tprod3(A,B)
[na,ra,n3] = size(A);
[rb,nb,n3b] = size(B);
if ra~=rb || n3~=n3b, error('t-product dimension mismatch.'); end
Af=fft(A,[],3); Bf=fft(B,[],3);
Cf=zeros(na,nb,n3);
for k=1:n3
    Cf(:,:,k)=Af(:,:,k)*Bf(:,:,k);
end
C=real(ifft(Cf,[],3));
end
