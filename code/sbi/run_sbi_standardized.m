function T = run_sbi_standardized()

% RUN_SBI_STANDARDIZED
%
% Final SBI runner for converted MAT files:
%
% Supported format:
%
%   Yuint8 : H x W x T uint8
%
% Optional:
%
%   BgtUint8 : H x W uint8
%
% The converted SBI files in this project contain:
%
%   SBI_Board.mat
%   SBI_CAVIAR1.mat
%   SBI_Hall_Monitor.mat
%   SBI_People_Foliage.mat
%
C = HCL_DSP_FINAL_config();
assert(~isempty(C.real.sbiRoot), ...
    'SBI root is empty.');
outDir = fullfile( ...
    C.outputRoot,...
    'SBI_FINAL');

if ~exist(outDir,'dir')
    mkdir(outDir);
end
rows = {};

rid = 0;

for i = 1:numel(C.real.sbiNames)
    name = C.real.sbiNames{i};
    f = fullfile( ...
        C.real.sbiRoot,...
        [name '.mat']);
    fprintf('\n');
    fprintf('================================================\n');
    fprintf(' Loading SBI sequence: %s\n',name);
    fprintf('================================================\n');
    assert(isfile(f), ...
        'Missing SBI file: %s',f);
    S = load(f);
    %% -------------------------------------------------
    % Load video tensor
    %% -------------------------------------------------
    if isfield(S,'Yuint8')
        Y = double(S.Yuint8)/255;
    elseif isfield(S,'Y')
        Y = double(S.Y);
        if max(Y(:)) > 1
            Y = Y/255;
        end
    else

        error( ...
            'SBI file %s does not contain Yuint8 or Y.', ...
            f);
    end
    fprintf('Video size = %d x %d x %d\n', ...
        size(Y,1), ...
        size(Y,2), ...
        size(Y,3));
    %% -------------------------------------------------
    % GT
    %% -------------------------------------------------
    hasGT=false;
    GT=[];
    if isfield(S,'BgtUint8')
        GT = double(S.BgtUint8)/255;
        hasGT=true;
    elseif isfield(S,'GT')
        GT = double(S.GT);
        if max(GT(:))>1
            GT=GT/255;
        end
        hasGT=true;
    end
    if hasGT
        fprintf('GT available: yes\n');
    else
        fprintf('GT available: no\n');
    end
    %% -------------------------------------------------
    % rank
    %% -------------------------------------------------

    r = energy_rank_fft( ...
        Y,...
        0.98,...
        15);
    fprintf('Rank = %d\n',r);



    %% -------------------------------------------------
    % HCL final
    %% -------------------------------------------------

    H = hcl_trpca_final( ...
        Y,...
        r,...
        C.hcl);



    %% -------------------------------------------------
    % Metrics
    %% -------------------------------------------------

    PSNR=NaN;
    SSIMv=NaN;
    AGE=NaN;


    if hasGT

        Bg = median(H.Xhat,3);


        Bg = min(max(Bg,0),1);


        mse = mean((Bg(:)-GT(:)).^2);


        PSNR = 10*log10(1/max(mse,eps));


        AGE = mean(abs(Bg(:)-GT(:)));


        if exist('ssim','file')

            SSIMv = ssim(Bg,GT);

        end

    end



    %% -------------------------------------------------
    % Save
    %% -------------------------------------------------

    R = struct();

    R.name=name;
    R.Y=Y;
    R.H=H;
    R.rank=r;

    R.hasGT=hasGT;

    if hasGT
        R.GT=GT;
    end


    save( ...
        fullfile(outDir,[name '_HCL_FINAL.mat']),...
        'R',...
        '-v7.3');


    rid=rid+1;


    rows(rid,:)={ ...
        name,...
        size(Y,1),...
        size(Y,2),...
        size(Y,3),...
        r,...
        H.selectedFraction,...
        H.recoveryESS,...
        H.iterations,...
        H.finalRelChange,...
        H.toleranceMet,...
        hasGT,...
        PSNR,...
        SSIMv,...
        AGE};



end



T = cell2table(rows,...
    'VariableNames',...
    {'Dataset',...
     'Height',...
     'Width',...
     'Frames',...
     'Rank',...
     'SelectedFraction',...
     'rESS',...
     'Iterations',...
     'FinalRelChange',...
     'ToleranceMet',...
     'HasGT',...
     'PSNR',...
     'SSIM',...
     'AGE'});


writetable( ...
    T,...
    fullfile(outDir,'sbi_final_summary.csv'));


disp(T);


end



function r = energy_rank_fft(X,target,cap)


Xf = fft(X,[],3);


m=min(size(X,1),size(X,2));


E=zeros(m,1);


for k=1:size(Xf,3)

    s=svd(Xf(:,:,k),'econ');

    E(1:numel(s))= ...
        E(1:numel(s))+abs(s(:)).^2;

end


c=cumsum(E)/sum(E);


r=find(c>=target,1,'first');


if isempty(r)
    r=m;
end


r=min(r,cap);


end