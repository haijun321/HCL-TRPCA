function T = run_sbi_standardized()
%RUN_SBI_STANDARDIZED Run the four SBI structural diagnostics.
%
% The SBI experiments in this release are structural and operational
% diagnostics. The distributed SBI ground-truth images audited for Board
% and CAVIAR1 are frame-wise foreground segmentation masks, not clean
% background reconstruction references. Their presence must therefore not
% trigger PSNR, SSIM, or AGE evaluation.
%
% Observation fields:
%   Yuint8 : H x W x T uint8
%   Y      : H x W x T numeric tensor, normally scaled to [0,1]
%
% Optional semantic metadata:
%   GroundTruthType = 'foreground_segmentation_mask'
%   HasSegmentationGT = true
%   ValidReconstructionReference = false
%
% A user-supplied clean background is evaluated only when BOTH
% ValidReconstructionReference is true and GroundTruthType is exactly
% 'clean_background_reference'. The image may then be supplied as
% CleanBackgroundUint8, CleanBackground, BgtUint8, or GT. The last two
% names are supported only for compatibility with legacy converted files.

C = HCL_DSP_FINAL_config();
assert(~isempty(C.real.sbiRoot), 'SBI root is empty.');

outDir = fullfile(C.outputRoot, 'SBI_FINAL');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

rows = {};
rid = 0;

for i = 1:numel(C.real.sbiNames)
    name = C.real.sbiNames{i};
    f = fullfile(C.real.sbiRoot, [name '.mat']);

    fprintf('\n');
    fprintf('================================================\n');
    fprintf(' Loading SBI sequence: %s\n', name);
    fprintf('================================================\n');
    assert(isfile(f), 'Missing SBI file: %s', f);

    S = load(f);
    Y = load_video_tensor(S, f);
    fprintf('Video size = %d x %d x %d\n', ...
        size(Y,1), size(Y,2), size(Y,3));

    semantic = inspect_reference_semantics(S, Y, name);
    fprintf('Ground-truth type: %s\n', semantic.groundTruthType);
    fprintf('Segmentation annotation available: %s\n', ...
        yes_no(semantic.hasSegmentationGT));
    fprintf('Valid clean reconstruction reference: %s\n', ...
        yes_no(semantic.validReconstructionReference));

    r = energy_rank_fft(Y, 0.98, 15);
    fprintf('Rank = %d\n', r);

    H = hcl_trpca_final(Y, r, C.hcl);

    PSNR = NaN;
    SSIMv = NaN;
    AGE = NaN;

    if semantic.validReconstructionReference
        Bg = median(H.Xhat, 3);
        Bg = min(max(Bg, 0), 1);
        GT = semantic.cleanReference;

        mse = mean((Bg(:) - GT(:)).^2);
        PSNR = 10 * log10(1 / max(mse, eps));
        AGE = mean(abs(Bg(:) - GT(:)));

        if exist('ssim', 'file') == 2
            SSIMv = ssim(Bg, GT, 'DynamicRange', 1);
        end
    end

    R = struct();
    R.name = name;
    R.Y = Y;
    R.H = H;
    R.rank = r;
    R.hasAnnotation = semantic.hasAnnotation;
    R.hasSegmentationGT = semantic.hasSegmentationGT;
    R.spatialAlignmentVerified = semantic.spatialAlignmentVerified;
    R.frameIndexAlignmentVerified = ...
        semantic.frameIndexAlignmentVerified;
    R.groundTruthType = semantic.groundTruthType;
    R.hasCleanReference = semantic.hasCleanReference;
    R.validReconstructionReference = ...
        semantic.validReconstructionReference;
    R.reconstructionMetricsStatus = ...
        semantic.reconstructionMetricsStatus;
    R.PSNR = PSNR;
    R.SSIM = SSIMv;
    R.AGE = AGE;

    if semantic.hasAnnotation
        R.annotation = semantic.annotation;
    end
    if semantic.validReconstructionReference
        R.cleanReference = semantic.cleanReference;
    end

    save(fullfile(outDir, [name '_HCL_FINAL.mat']), 'R', '-v7.3');

    rid = rid + 1;
    rows(rid,:) = { ...
        name, ...
        size(Y,1), ...
        size(Y,2), ...
        size(Y,3), ...
        r, ...
        H.selectedFraction, ...
        H.recoveryESS, ...
        H.iterations, ...
        H.finalRelChange, ...
        H.toleranceMet, ...
        semantic.hasAnnotation, ...
        semantic.hasSegmentationGT, ...
        semantic.spatialAlignmentVerified, ...
        semantic.frameIndexAlignmentVerified, ...
        semantic.hasCleanReference, ...
        semantic.validReconstructionReference, ...
        semantic.groundTruthType, ...
        PSNR, ...
        SSIMv, ...
        AGE, ...
        semantic.reconstructionMetricsStatus};
end

T = cell2table(rows, 'VariableNames', { ...
    'Dataset', ...
    'Height', ...
    'Width', ...
    'Frames', ...
    'Rank', ...
    'SelectedFraction', ...
    'rESS', ...
    'Iterations', ...
    'FinalRelChange', ...
    'ToleranceMet', ...
    'HasAnnotation', ...
    'HasSegmentationGT', ...
    'SpatialAlignmentVerified', ...
    'FrameIndexAlignmentVerified', ...
    'HasCleanReference', ...
    'ValidReconstructionReference', ...
    'GroundTruthType', ...
    'PSNR', ...
    'SSIM', ...
    'AGE', ...
    'ReconstructionMetricsStatus'});

writetable(T, fullfile(outDir, 'sbi_final_summary.csv'));
disp(T);
end


function Y = load_video_tensor(S, sourceFile)
if isfield(S, 'Yuint8')
    Y = double(S.Yuint8) / 255;
elseif isfield(S, 'Y')
    Y = double(S.Y);
    if max(Y(:)) > 1
        Y = Y / 255;
    end
else
    error('SBI file %s does not contain Yuint8 or Y.', sourceFile);
end

assert(ndims(Y) == 3, ...
    'SBI observation must be an H-by-W-by-T grayscale tensor.');
assert(all(isfinite(Y(:))), ...
    'SBI observation contains non-finite values.');
end


function info = inspect_reference_semantics(S, Y, sequenceName)
info = struct();
info.hasAnnotation = false;
info.hasSegmentationGT = false;
info.spatialAlignmentVerified = false;
info.frameIndexAlignmentVerified = false;
info.groundTruthType = read_text_field( ...
    S, 'GroundTruthType', 'none');
info.hasCleanReference = false;
info.validReconstructionReference = false;
info.reconstructionMetricsStatus = ...
    'not_applicable_no_valid_clean_background_reference';
info.annotation = [];
info.cleanReference = [];

[annotation, annotationFound] = first_numeric_field(S, { ...
    'SegmentationGTUint8', 'SegmentationGT'});

[legacyCandidate, legacyFound] = first_numeric_field(S, { ...
    'BgtUint8', 'GT'});

auditedMaskSequence = any(strcmp(sequenceName, { ...
    'SBI_Board', 'SBI_CAVIAR1'}));

if annotationFound
    info.hasSegmentationGT = true;
    info.groundTruthType = 'foreground_segmentation_mask';
elseif strcmpi(info.groundTruthType, 'foreground_segmentation_mask')
    info.hasSegmentationGT = true;
elseif read_logical_scalar(S, 'HasSegmentationGT', false)
    info.hasSegmentationGT = true;
    info.groundTruthType = 'foreground_segmentation_mask';
elseif auditedMaskSequence && legacyFound
    % v1.0.1 semantic clarification for the two audited sequences.
    info.hasSegmentationGT = true;
    info.groundTruthType = 'foreground_segmentation_mask';
end

if annotationFound
    info.annotation = normalize_image(annotation, ...
        'segmentation annotation');
    info.hasAnnotation = true;
elseif info.hasSegmentationGT && legacyFound
    info.annotation = normalize_image(legacyCandidate, ...
        'legacy segmentation annotation');
    info.hasAnnotation = true;
end

if info.hasAnnotation
    annotationSize = size(info.annotation);
    info.spatialAlignmentVerified = ...
        numel(annotationSize) >= 2 && ...
        annotationSize(1) == size(Y,1) && ...
        annotationSize(2) == size(Y,2);
    info.frameIndexAlignmentVerified = ...
        ndims(info.annotation) == 3 && ...
        size(info.annotation,3) == size(Y,3);
end

validityRequested = read_logical_scalar( ...
    S, 'ValidReconstructionReference', false);
cleanTypeDeclared = strcmpi( ...
    info.groundTruthType, 'clean_background_reference');

[cleanReference, cleanFound] = first_numeric_field(S, { ...
    'CleanBackgroundUint8', 'CleanBackground'});
if ~cleanFound && validityRequested && cleanTypeDeclared && legacyFound
    cleanReference = legacyCandidate;
    cleanFound = true;
end

info.hasCleanReference = cleanFound && cleanTypeDeclared;

if validityRequested && cleanTypeDeclared && cleanFound
    cleanReference = normalize_image(cleanReference, ...
        'clean background reference');
    assert(ismatrix(cleanReference), ...
        'Clean background reference must be a two-dimensional image.');
    assert(isequal(size(cleanReference), [size(Y,1), size(Y,2)]), ...
        'Clean background reference does not match the video dimensions.');

    info.cleanReference = cleanReference;
    info.validReconstructionReference = true;
    info.spatialAlignmentVerified = true;
    info.reconstructionMetricsStatus = ...
        'computed_from_explicitly_validated_clean_background_reference';
elseif validityRequested
    warning('HCL:SBIReferenceDisabled', [ ...
        'Reference-based metrics were disabled for %s because both ' ...
        'GroundTruthType=''clean_background_reference'' and a numeric ' ...
        'clean-background field are required.'], sequenceName);
end
end


function [value, found] = first_numeric_field(S, fieldNames)
value = [];
found = false;

for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if isfield(S, fieldName) && ...
            (isnumeric(S.(fieldName)) || islogical(S.(fieldName)))
        value = S.(fieldName);
        found = true;
        return;
    end
end
end


function X = normalize_image(X, label)
X = double(X);
assert(all(isfinite(X(:))), '%s contains non-finite values.', label);
if max(X(:)) > 1
    X = X / 255;
end
end


function value = read_logical_scalar(S, fieldName, defaultValue)
value = defaultValue;
if ~isfield(S, fieldName)
    return;
end

candidate = S.(fieldName);
if islogical(candidate) && isscalar(candidate)
    value = logical(candidate);
elseif isnumeric(candidate) && isscalar(candidate) && ...
        isfinite(candidate) && any(candidate == [0, 1])
    value = logical(candidate);
else
    warning('HCL:InvalidSemanticMetadata', ...
        '%s must be a logical or numeric scalar; using false.', fieldName);
    value = false;
end
end


function value = read_text_field(S, fieldName, defaultValue)
value = defaultValue;
if ~isfield(S, fieldName)
    return;
end

candidate = S.(fieldName);
if ischar(candidate)
    value = strtrim(candidate);
elseif isstring(candidate) && isscalar(candidate)
    value = strtrim(char(candidate));
else
    warning('HCL:InvalidSemanticMetadata', ...
        '%s must be a character vector or scalar string.', fieldName);
end
end


function value = yes_no(tf)
if tf
    value = 'yes';
else
    value = 'no';
end
end


function r = energy_rank_fft(X, target, cap)
Xf = fft(X, [], 3);
m = min(size(X,1), size(X,2));
E = zeros(m,1);

for k = 1:size(Xf,3)
    s = svd(Xf(:,:,k), 'econ');
    E(1:numel(s)) = E(1:numel(s)) + abs(s(:)).^2;
end

totalEnergy = sum(E);
assert(totalEnergy > 0 && isfinite(totalEnergy), ...
    'Cannot estimate the rank of a zero-energy or invalid tensor.');
c = cumsum(E) / totalEnergy;
r = find(c >= target, 1, 'first');

if isempty(r)
    r = m;
end
r = min(r, cap);
end
