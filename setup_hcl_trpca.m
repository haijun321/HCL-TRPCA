function repoRoot = setup_hcl_trpca()
%SETUP_HCL_TRPCA Add the public HCL-TRPCA repository to the MATLAB path.

repoRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(repoRoot,'configs'));
addpath(genpath(fullfile(repoRoot,'code')));
addpath(fullfile(repoRoot,'scripts'));

fprintf('HCL-TRPCA repository added to MATLAB path:\n%s\n',repoRoot);
end
