function R = run_realdata_standardized(kind)
% RUN_REALDATA_STANDARDIZED
%
% Paper-final standardized real-data runner for HCL-TRPCA.
%
% Supported modes:
%
%   R = run_realdata_standardized('PAVIA');
%   R = run_realdata_standardized('URBAN');
%
% -------------------------------------------------------------------------
% Pavia University
% -------------------------------------------------------------------------
% Expected source:
%   PaviaU.mat containing a numeric 3-D cube, typically 610 x 340 x 103.
%
% Paper protocol:
%   crop rows 178:433
%   crop cols  43:298
%   -> 256 x 256 x 103
%   divide by 8000 and clip to [0,1]
%   rank: 98% Fourier energy, cap 25
%
% -------------------------------------------------------------------------
% HYDICE Urban
% -------------------------------------------------------------------------
% Verified source format:
%
%   Y      : 210 x 94249
%   nRow   : 307
%   nCol   : 307
%   nBand  : 210
%
% where each ROW Y(k,:) is one vectorized spectral band.
%
% Verified reconstruction:
%
%   Y
%   -> 307 x 307 x 210
%   -> global min-max normalization
%   -> crop rows/cols 26:281
%   -> 256 x 256 x 210
%
% All 210 bands are retained.
% SlectBands is stored only as source metadata and is NOT used.
%
% Urban paper rank:
%   r = 25
%
% Urban recovery budget:
%   80 total HCL updates
%
% -------------------------------------------------------------------------

C = HCL_DSP_FINAL_config();

%% Make sure output directory exists
if ~exist(C.outputRoot,'dir')
    mkdir(C.outputRoot);
end

switch upper(kind)

    %% ====================================================================
    %  PAVIA UNIVERSITY
    % =====================================================================
    case 'PAVIA'

        assert(~isempty(C.real.paviaMat), ...
            'Set C.real.paviaMat first.');

        assert(isfile(C.real.paviaMat), ...
            'Pavia MAT file does not exist: %s', ...
            C.real.paviaMat);

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf(' Loading Pavia University FINAL protocol\n');
        fprintf('============================================================\n');

        %% ---------------------------------------------------------------
        % 1. Load largest numeric 3-D cube
        %% ---------------------------------------------------------------

        Xraw = largest_3d_numeric(C.real.paviaMat);

        fprintf('Raw Pavia cube = %d x %d x %d\n', ...
            size(Xraw,1), ...
            size(Xraw,2), ...
            size(Xraw,3));

        assert(size(Xraw,1) >= 433 && ...
               size(Xraw,2) >= 298, ...
            'Pavia raw cube too small.');

        %% ---------------------------------------------------------------
        % 2. Paper-fixed crop
        %% ---------------------------------------------------------------

        X = Xraw(178:433,43:298,:);

        clear Xraw

        fprintf('Cropped Pavia cube = %d x %d x %d\n', ...
            size(X,1), ...
            size(X,2), ...
            size(X,3));

        assert(size(X,1) == 256 && ...
               size(X,2) == 256, ...
            'Unexpected Pavia crop size.');

        %% ---------------------------------------------------------------
        % 3. Paper normalization
        %% ---------------------------------------------------------------

        X = double(X) ./ 8000;

        X = min(max(X,0),1);

        fprintf('Pavia normalized range = [%g, %g]\n', ...
            min(X(:)),max(X(:)));

        %% ---------------------------------------------------------------
        % 4. Fourier-energy rank
        %% ---------------------------------------------------------------

        r = energy_rank_fft(X,0.98,25);

        fprintf('Pavia fitted rank = %d\n',r);

        %% ---------------------------------------------------------------
        % 5. Final frozen HCL-TRPCA
        %% ---------------------------------------------------------------

        Ppavia = C.hcl;

        fprintf('Pavia epsRec = %.4g\n', ...
            Ppavia.epsRec);

        H = hcl_trpca_final( ...
            X, ...
            r, ...
            Ppavia);

        %% ---------------------------------------------------------------
        % 6. Save result
        %% ---------------------------------------------------------------

        R = struct();

        R.X = X;
        R.rank = r;
        R.H = H;

        R.protocol = Ppavia;

        R.sourceFile = C.real.paviaMat;

        R.rawDescription = ...
            'Pavia University source cube';

        R.workingSize = size(X);

        R.cropRows = 178:433;
        R.cropCols = 43:298;

        R.normalization = ...
            'divide by 8000 and clip to [0,1]';

        outFile = fullfile( ...
            C.outputRoot, ...
            'PAVIA_FINAL.mat');

        save( ...
            outFile, ...
            'R', ...
            '-v7.3');

        fprintf('\n');
        fprintf('Pavia final result saved:\n');
        fprintf('%s\n',outFile);

        fprintf('\n');
        fprintf('Pavia diagnostics:\n');
        fprintf('selected fraction = %.6f\n', ...
            H.selectedFraction);
        fprintf('recovery ESS      = %.6f\n', ...
            H.recoveryESS);
        fprintf('iterations        = %d\n', ...
            H.iterations);
        fprintf('final rel change  = %.6e\n', ...
            H.finalRelChange);
        fprintf('tolerance met     = %d\n', ...
            H.toleranceMet);


    %% ====================================================================
    %  HYDICE URBAN
    % =====================================================================
    case 'URBAN'

        assert(~isempty(C.real.urbanMat), ...
            'Set C.real.urbanMat first.');

        assert(isfile(C.real.urbanMat), ...
            'Urban MAT file does not exist: %s', ...
            C.real.urbanMat);

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf(' Loading HYDICE Urban FINAL protocol\n');
        fprintf('============================================================\n');

        %% ---------------------------------------------------------------
        % 1. Load verified source
        %% ---------------------------------------------------------------

        S = load(C.real.urbanMat);

        assert(isfield(S,'Y'), ...
            'Urban file must contain variable Y.');

        assert(isfield(S,'nRow') && ...
               isfield(S,'nCol') && ...
               isfield(S,'nBand'), ...
            'Urban file must contain nRow, nCol, and nBand.');

        nRow  = double(S.nRow);
        nCol  = double(S.nCol);
        nBand = double(S.nBand);

        assert(isscalar(nRow) && ...
               isscalar(nCol) && ...
               isscalar(nBand), ...
            'Urban metadata must be scalar.');

        %% Verified source dimensions
        assert(nRow == 307 && ...
               nCol == 307 && ...
               nBand == 210, ...
            'Unexpected Urban metadata: %d x %d x %d.', ...
            nRow,nCol,nBand);

        Y2 = S.Y;

        assert(isnumeric(Y2) && ...
               ismatrix(Y2), ...
            'Urban Y must be a numeric 2-D matrix.');

        %% Verified layout:
        % 210 bands x 307*307 pixels
        assert(size(Y2,1) == 210 && ...
               size(Y2,2) == 307*307, ...
            ['Urban Y must have verified bands-by-pixels layout ', ...
             '210 x 94249.']);

        assert(all(isfinite(Y2(:))), ...
            'Urban Y contains nonfinite values.');

        fprintf('Stored Y size = %d x %d\n', ...
            size(Y2,1), ...
            size(Y2,2));

        fprintf('Metadata      = %d x %d x %d\n', ...
            nRow,nCol,nBand);

        fprintf('Verified layout = bands x pixels\n');

        %% ---------------------------------------------------------------
        % 2. Preserve source metadata before clearing S
        %% ---------------------------------------------------------------

        if isfield(S,'SlectBands')
            sourceSelectedBands = S.SlectBands;
        else
            sourceSelectedBands = [];
        end

        if isfield(S,'maxValue')
            sourceMaxValue = S.maxValue;
        else
            sourceMaxValue = [];
        end

        %% ---------------------------------------------------------------
        % 3. Reconstruct 307 x 307 x 210 cube
        %
        % IMPORTANT:
        % each row Y(k,:) is one vectorized spectral band.
        %
        % The orientation below has already been visually verified
        % using Urban bands 50 and 139.
        %% ---------------------------------------------------------------

        Xraw = zeros( ...
            nRow, ...
            nCol, ...
            nBand, ...
            'double');

        for k = 1:nBand

            Xraw(:,:,k) = reshape( ...
                double(Y2(k,:)), ...
                nRow, ...
                nCol);

        end

        clear Y2
        clear S

        fprintf('Reconstructed Urban cube = %d x %d x %d\n', ...
            size(Xraw,1), ...
            size(Xraw,2), ...
            size(Xraw,3));

        assert(isequal( ...
            size(Xraw), ...
            [307 307 210]), ...
            'Urban cube reconstruction failed.');

        %% ---------------------------------------------------------------
        % 4. Global min-max normalization BEFORE cropping
        %
        % This ordering is fixed by the manuscript protocol.
        %% ---------------------------------------------------------------

        xmin = min(Xraw(:));
        xmax = max(Xraw(:));

        assert(isfinite(xmin) && ...
               isfinite(xmax) && ...
               xmax > xmin, ...
            'Invalid Urban intensity range.');

        fprintf('Urban raw range = [%g, %g]\n', ...
            xmin,xmax);

        Xraw = ...
            (Xraw - xmin) ./ ...
            (xmax - xmin);

        fprintf('Urban normalized range = [%g, %g]\n', ...
            min(Xraw(:)), ...
            max(Xraw(:)));

        %% ---------------------------------------------------------------
        % 5. Paper-fixed crop:
        % rows 26:281
        % cols 26:281
        %% ---------------------------------------------------------------

        X = Xraw(26:281,26:281,:);

        clear Xraw

        fprintf('Working Urban cube = %d x %d x %d\n', ...
            size(X,1), ...
            size(X,2), ...
            size(X,3));

        assert(isequal( ...
            size(X), ...
            [256 256 210]), ...
            'Unexpected Urban working size.');

        %% ---------------------------------------------------------------
        % 6. Paper-fixed rank
        %% ---------------------------------------------------------------

        r = 25;

        fprintf('Urban fitted rank = %d\n',r);

        %% ---------------------------------------------------------------
        % 7. Urban-specific final HCL protocol
        %% ---------------------------------------------------------------

        Purban = C.hcl;

        %% Paper-fixed Urban budget
        Purban.maxIterations = 80;

        fprintf('Urban max iterations = %d\n', ...
            Purban.maxIterations);

        fprintf('Urban epsWarm = %.4g\n', ...
            Purban.epsWarm);

        fprintf('Urban epsSel  = %.4g\n', ...
            Purban.epsSel);

        fprintf('Urban epsRec  = %.4g\n', ...
            Purban.epsRec);

        %% ---------------------------------------------------------------
        % 8. Final frozen HCL-TRPCA
        %% ---------------------------------------------------------------

        H = hcl_trpca_final( ...
            X, ...
            r, ...
            Purban);

        %% ---------------------------------------------------------------
        % 9. Save complete final result
        %% ---------------------------------------------------------------

        R = struct();

        R.X = X;
        R.rank = r;
        R.H = H;

        R.protocol = Purban;

        R.sourceFile = C.real.urbanMat;

        R.sourceLayout = ...
            'bands_x_pixels_210x94249';

        R.reshapeConvention = ...
            ['each Y(k,:) reshaped to 307x307 ', ...
             'using MATLAB column-major order'];

        R.rawSize = ...
            [307 307 210];

        R.workingSize = ...
            [256 256 210];

        R.cropRows = ...
            26:281;

        R.cropCols = ...
            26:281;

        R.normalization = ...
            'global min-max before crop';

        R.rawMin = xmin;
        R.rawMax = xmax;

        %% These fields are stored ONLY as metadata.
        % They are NOT used to remove spectral bands.
        R.sourceSelectedBands = ...
            sourceSelectedBands;

        R.sourceMaxValue = ...
            sourceMaxValue;

        %% ---------------------------------------------------------------
        % 10. Save
        %% ---------------------------------------------------------------

        outFile = fullfile( ...
            C.outputRoot, ...
            'URBAN_FINAL.mat');

        save( ...
            outFile, ...
            'R', ...
            '-v7.3');

        fprintf('\n');
        fprintf('Urban final result saved:\n');
        fprintf('%s\n',outFile);

        %% ---------------------------------------------------------------
        % 11. Console diagnostics
        %% ---------------------------------------------------------------

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf(' Urban FINAL diagnostics\n');
        fprintf('============================================================\n');

        fprintf('working size      = %d x %d x %d\n', ...
            size(X,1), ...
            size(X,2), ...
            size(X,3));

        fprintf('rank              = %d\n', ...
            r);

        fprintf('epsRec            = %.6g\n', ...
            H.epsRec);

        fprintf('selected fraction = %.8f\n', ...
            H.selectedFraction);

        fprintf('recovery ESS      = %.8f\n', ...
            H.recoveryESS);

        fprintf('iterations        = %d\n', ...
            H.iterations);

        fprintf('final rel change  = %.8e\n', ...
            H.finalRelChange);

        fprintf('tolerance met     = %d\n', ...
            H.toleranceMet);

        fprintf('============================================================\n');


    %% ====================================================================
    %  INVALID MODE
    % =====================================================================
    otherwise

        error( ...
            ['Unknown real-data mode "%s". ', ...
             'Use PAVIA or URBAN. ', ...
             'SBI is handled by run_sbi_standardized.m.'], ...
             kind);

end

end


%% =========================================================================
%  HELPER: LARGEST NUMERIC 3-D ARRAY
% ==========================================================================
function X = largest_3d_numeric(file)

S = load(file);

fn = fieldnames(S);

bestN = 0;
X = [];

for i = 1:numel(fn)

    v = S.(fn{i});

    if isnumeric(v) && ...
       ndims(v) == 3 && ...
       numel(v) > bestN

        X = v;
        bestN = numel(v);

    end

end

assert(~isempty(X), ...
    'No numeric 3-D cube found in %s', ...
    file);

end


%% =========================================================================
%  HELPER: FOURIER-ENERGY RANK
% ==========================================================================
function r = energy_rank_fft(X,target,cap)

Xf = fft(X,[],3);

m = min( ...
    size(X,1), ...
    size(X,2));

E = zeros(m,1);

for k = 1:size(Xf,3)

    s = svd( ...
        Xf(:,:,k), ...
        'econ');

    E(1:numel(s)) = ...
        E(1:numel(s)) + ...
        abs(s(:)).^2;

end

totalEnergy = sum(E);

assert(isfinite(totalEnergy) && ...
       totalEnergy > 0, ...
    'Degenerate Fourier singular-value energy.');

c = cumsum(E) ./ totalEnergy;

r = find( ...
    c >= target, ...
    1, ...
    'first');

if isempty(r)
    r = m;
end

r = min(r,cap);

end