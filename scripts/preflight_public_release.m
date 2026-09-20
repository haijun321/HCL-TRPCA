function report = preflight_public_release()
%PREFLIGHT_PUBLIC_RELEASE Fresh-clone dependency and configuration check.

report = struct();

repoRoot = setup_hcl_trpca();
C = HCL_DSP_FINAL_config();

report.repoRoot = repoRoot;
report.configOK = true;

report.hclCore = which('hcl_trpca_final');
report.syntheticEntry = which('run_HCL_DSP_FINAL');
report.paviaEntry = which('run_PaviaU_CONTROLLED_FINAL_V3_2');
report.urbanEntry = which('run_realdata_standardized');
report.sbiEntry = which('run_sbi_standardized');
report.sbiAuditEntry = which('audit_sbi_release_semantics');
report.tnnAdapter = which('my_TNN_adapter');
report.ptrpcaAdapter = which('my_pTRPCA_adapter');

report.tnnCore = which('hcldref.tnn_trpca_admm');
report.ptrpcaCore = which('ptrpca');

fprintf('\n=== HCL-TRPCA PUBLIC PREFLIGHT ===\n');
fprintf('Config: PASS\n');
fprintf('HCL core: %s\n',blank_if_empty(report.hclCore));
fprintf('Synthetic entry: %s\n',blank_if_empty(report.syntheticEntry));
fprintf('Pavia entry: %s\n',blank_if_empty(report.paviaEntry));
fprintf('Urban entry: %s\n',blank_if_empty(report.urbanEntry));
fprintf('SBI entry: %s\n',blank_if_empty(report.sbiEntry));
fprintf('SBI semantic audit: %s\n',blank_if_empty(report.sbiAuditEntry));
fprintf('TNN adapter: %s\n',blank_if_empty(report.tnnAdapter));
fprintf('p-TRPCA adapter: %s\n',blank_if_empty(report.ptrpcaAdapter));

if isempty(report.tnnCore)
    fprintf('TNN core: MISSING (hcldref.tnn_trpca_admm)\n');
else
    fprintf('TNN core: %s\n',report.tnnCore);
end

if isempty(report.ptrpcaCore)
    fprintf('p-TRPCA core: MISSING (ptrpca)\n');
else
    fprintf('p-TRPCA core: %s\n',report.ptrpcaCore);
end

report.passCore = ~isempty(report.hclCore) && ...
                  ~isempty(report.syntheticEntry) && ...
                  ~isempty(report.paviaEntry) && ...
                  ~isempty(report.urbanEntry) && ...
                  ~isempty(report.sbiEntry) && ...
                  ~isempty(report.sbiAuditEntry);

report.passBaselines = ~isempty(report.tnnCore) && ~isempty(report.ptrpcaCore);

fprintf('Core repository preflight: %d\n',report.passCore);
fprintf('Baseline dependency preflight: %d\n',report.passBaselines);

if ~isempty(report.sbiAuditEntry)
    report.sbiSemantics = audit_sbi_release_semantics();
else
    report.sbiSemantics = struct('pass',false);
end
end

function s = blank_if_empty(s)
if isempty(s)
    s = 'MISSING';
end
end
