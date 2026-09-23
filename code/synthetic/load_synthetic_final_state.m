function X = load_synthetic_final_state(file)
%LOAD_SYNTHETIC_FINAL_STATE
% Strict loader for frozen synthetic HCL states.
% All attribution quantities are read explicitly from H.final.

S = load(file,'D','H','M');

assert(isfield(S,'D'), 'Missing D in %s',file);
assert(isfield(S,'H'), 'Missing H in %s',file);
assert(isfield(S,'M'), 'Missing M in %s',file);
assert(isfield(S.H,'final'), 'Missing H.final in %s',file);

F = S.H.final;

requiredFinal = { ...
    'ae','ab','as', ...
    'po', ...
    'qe','qb','qs', ...
    'piEntry','piBlock','piSlice'};

for k = 1:numel(requiredFinal)
    assert(isfield(F,requiredFinal{k}), ...
        'Missing H.final.%s in %s',requiredFinal{k},file);
end

% Final-stage evidence and attribution
X.ae = double(F.ae);
X.ab = double(F.ab);
X.as = double(F.as);
X.po = double(F.po);
X.qe = double(F.qe);
X.qb = double(F.qb);
X.qs = double(F.qs);
X.piEntry = double(F.piEntry);
X.piBlock = double(F.piBlock);
X.piSlice = double(F.piSlice);

% Ground truth
assert(isfield(S.D,'truth'), 'Missing D.truth in %s',file);
assert(isfield(S.D.truth,'level'), 'Missing D.truth.level in %s',file);
assert(isfield(S.D.truth,'unionMask'), 'Missing D.truth.unionMask in %s',file);

X.gtLabel  = uint8(S.D.truth.level);
X.trueMask = logical(S.D.truth.unionMask);
X.gtEntry = X.gtLabel == 1;
X.gtBlock = X.gtLabel == 2;
X.gtSlice = X.gtLabel == 3;

% Exact final HCL detection support
X.detMask = X.po >= 0.50;

% Exact official HCL granularity labels
A = cat(4,X.piEntry,X.piBlock,X.piSlice);
[~,lab] = max(A,[],4);

X.predHCL = zeros(size(X.detMask),'uint8');
X.predHCL(X.detMask) = uint8(lab(X.detMask));

% Stored official metrics
X.officialDetF1  = double(S.M.DetF1);
X.officialGranF1 = double(S.M.LevelMacroF1);

X.file = string(file);
end
