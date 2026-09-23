function R = merge_confirmation_all_final()
% Merge original 30-seed confirmation and separately replayed frozen baselines.
% Original files are never overwritten.

C = HCL_DSP_FINAL_config();
outDir = fullfile(C.outputRoot,'E01_CONFIRMATION');
oldFile = fullfile(outDir,'confirmation_raw.csv');
baseFile = fullfile(outDir,'confirmation_baselines_final.csv');

assert(exist(oldFile,'file')==2,'Missing confirmation_raw.csv.');
assert(exist(baseFile,'file')==2,'Missing confirmation_baselines_final.csv.');

T0 = readtable(oldFile);
TB = readtable(baseFile);

% Normalize text columns
for name = ["Scenario","Method","InputHash","ReferenceHash"]
    if iscell(T0.(name)), T0.(name)=string(T0.(name)); else, T0.(name)=string(T0.(name)); end
    if iscell(TB.(name)), TB.(name)=string(TB.(name)); else, TB.(name)=string(TB.(name)); end
end

vars = {'Seed','Scenario','Method','NRE','DetectionF1','RecoveryMaskF1', ...
    'LevelMacroF1','rESS','SelectedFraction','InputHash','ReferenceHash','EpsRec'};
T0 = T0(:,vars);
TB = TB(:,vars);

% Verify every baseline row against the original pair hash
for i=1:height(TB)
    idx = T0.Seed==TB.Seed(i) & T0.Scenario==TB.Scenario(i);
    assert(any(idx),'No original pair for seed=%d scenario=%s.',TB.Seed(i),TB.Scenario(i));
    j = find(idx,1,'first');
    assert(T0.InputHash(j)==TB.InputHash(i),'InputHash mismatch at baseline row %d.',i);
    assert(T0.ReferenceHash(j)==TB.ReferenceHash(i),'ReferenceHash mismatch at baseline row %d.',i);
end

Tall = [T0;TB];

% Duplicate check
key = string(Tall.Seed)+"|"+Tall.Scenario+"|"+Tall.Method;
assert(numel(unique(key))==numel(key),'Duplicate Seed/Scenario/Method rows after merge.');

requiredMethods = ["tSVD-LRA","TNN-TRPCA","p-TRPCA","HCL-TRPCA","Oracle"];
scenarios = ["EB","ES","BS","EBS"];

for m = requiredMethods
    for sc = scenarios
        n = sum(Tall.Method==m & Tall.Scenario==sc);
        assert(n==30,'%s/%s count=%d; expected 30.',m,sc,n);
    end
end

writetable(Tall,fullfile(outDir,'confirmation_all_final.csv'));

% NRE summary
rows = {};
k=0;
for m = requiredMethods
    for sc = scenarios
        idx = Tall.Method==m & Tall.Scenario==sc;
        x = Tall.NRE(idx);
        k=k+1;
        rows(k,:) = {char(m),char(sc),numel(x),mean(x,'omitnan'),std(x,'omitnan')}; %#ok<AGROW>
    end
end
S = cell2table(rows,'VariableNames',{'Method','Scenario','N','MeanNRE','SDNRE'});
writetable(S,fullfile(outDir,'confirmation_all_final_summary.csv'));

% HCL diagnostics
TH = Tall(Tall.Method=="HCL-TRPCA",:);
drows = cell(numel(scenarios),8);
for i=1:numel(scenarios)
    sc=scenarios(i);
    idx=TH.Scenario==sc;
    drows(i,:)={char(sc),sum(idx), ...
        mean(TH.DetectionF1(idx),'omitnan'),std(TH.DetectionF1(idx),'omitnan'), ...
        mean(TH.RecoveryMaskF1(idx),'omitnan'),std(TH.RecoveryMaskF1(idx),'omitnan'), ...
        mean(TH.LevelMacroF1(idx),'omitnan'),std(TH.LevelMacroF1(idx),'omitnan')};
end
D = cell2table(drows,'VariableNames',{'Scenario','N', ...
    'MeanDetectionF1','SDDetectionF1','MeanRecoveryMaskF1','SDRecoveryMaskF1', ...
    'MeanLevelMacroF1','SDLevelMacroF1'});
writetable(D,fullfile(outDir,'confirmation_hcl_diagnostics_final.csv'));

save(fullfile(outDir,'confirmation_all_final_results.mat'),'Tall','S','D','-v7.3');

fprintf('\n===== FINAL 30-SEED NRE LATEX ROWS =====\n');
for m = requiredMethods
    fprintf('%s ',m);
    for sc = scenarios
        idx=Tall.Method==m & Tall.Scenario==sc;
        x=Tall.NRE(idx);
        fprintf('& $%.5f\\pm%.5f$ ',mean(x,'omitnan'),std(x,'omitnan'));
    end
    fprintf('\\\\\n');
end

fprintf('\n===== FINAL HCL DIAGNOSTIC LATEX ROWS =====\n');
for i=1:height(D)
    fprintf('%s & $%.4f\\pm%.4f$ & $%.4f\\pm%.4f$ & $%.4f\\pm%.4f$ \\\\\n', ...
        D.Scenario{i},D.MeanDetectionF1(i),D.SDDetectionF1(i), ...
        D.MeanRecoveryMaskF1(i),D.SDRecoveryMaskF1(i), ...
        D.MeanLevelMacroF1(i),D.SDLevelMacroF1(i));
end

R=struct('all',Tall,'summary',S,'diagnostics',D);
end
