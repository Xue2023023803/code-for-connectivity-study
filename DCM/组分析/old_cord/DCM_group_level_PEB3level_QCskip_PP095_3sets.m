
function DCM_group_level_PEB3level_QCskip_PP095_3sets()
% DCM_group_level_PEB3level_QCskip_PP095_3sets
%
% Outputs three 3-level (PEB-of-PEBs) results:
%   1) PEB3_connect     : connect-only (subject-level PEB over runs; then group PEB)
%   2) PEB3_unconnect   : unconnect-only (subject-level PEB over runs; then group PEB)
%   3) PEB3_condEffect  : within-subject condition effect (connect vs unconnect) using
%                         subject-level PEB with X2=[1, cond] then group PEB-of-PEBs
%
% QC rule (as requested): ONLY skip DCM files whose mean R^2 < 0.
% No other QC exclusions (no R^2 thresholding, no outlier rejection).
%
% Posterior probability threshold for exploratory reporting: PP > 0.95.
%
% Requirements:
%   - SPM12/SPM25 on MATLAB path (spm_dcm_peb, spm_dcm_peb_bmc, spm_BMS)
%   - run-wise DCM files already estimated
%   - QC mat (optional) from your QC script, used only to read per-file R2 if present
%
% Output:
%   root_prf/PEB_3level_QCskip_PP095_3sets/<ses-XX>/
%     metadata_DCM_index.mat
%     PEB2_connect_sub-XX.mat
%     PEB2_unconnect_sub-XX.mat
%     PEB2_condEffect_sub-XX.mat
%     PEB3_connect.mat,    BMA3_connect.mat
%     PEB3_unconnect.mat,  BMA3_unconnect.mat
%     PEB3_condEffect.mat, BMA3_condEffect.mat
%     report_*.txt
%
% Notes:
%   - PEB-of-PEBs is implemented exactly as in SPM PEB docs:
%       PEB1 = spm_dcm_peb(GCM1,X1);
%       PEB2 = spm_dcm_peb(GCM2,X2);
%       PEB3 = spm_dcm_peb({PEB1;PEB2},X3);
%
% Bao Xue / 2025-12-14

%% ---------------------- user settings ----------------------
root_prf      = '/home/xue/data/prf_result';
session_label = 'ses-01';
dcm_subdir    = 'DCM_manual_spm25_multiModel';

subject_list = { ...
    'sub-01','sub-02','sub-03','sub-04','sub-05', ...
    'sub-06','sub-07','sub-08','sub-09','sub-10'};

% Model names (optional, purely for printing)
modelNames = { ...
    'B1_globalGain', ...
    'B2_NPOhub', ...
    'B3_feedforward', ...
    'B4_feedback', ...
    'B5_parietalOnly', ...
    'B6_hierarchy', ...
    'B0_null'};

PP_THRESH = 0.95;   % exploratory threshold (unified)

% Which DCM parameter fields to take to PEB level (default: B only)
field = {'B'};   % change to {'A','B'} if you want both

% QC MAT (optional). Used ONLY to fetch per-file mean R2 quickly if present.
qc_mat = fullfile(root_prf, 'DCM_QC', sprintf('DCM_QC_runwise_BMS_%s.mat', session_label));

% Output dir
out_dir = fullfile(root_prf, 'PEB_3level_QCskip_PP095_3sets', session_label);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% ---------------------- sanity check: SPM ----------------------
if exist('spm', 'file') ~= 2 || exist('spm_dcm_peb', 'file') ~= 2
    error('SPM not found on path. Please add SPM12/SPM25 to MATLAB path first.');
end

%% ---------------------- load QC map (optional) ----------------------
qcR2Map = containers.Map('KeyType','char','ValueType','double');
if exist(qc_mat, 'file')
    tmp = load(qc_mat, 'QC');
    if isfield(tmp, 'QC') && ~isempty(tmp.QC) && isfield(tmp.QC, 'file') && isfield(tmp.QC, 'R2')
        for i = 1:numel(tmp.QC)
            f = tmp.QC(i).file;
            if ischar(f) && ~isempty(f)
                qcR2Map(normalize_path(f)) = tmp.QC(i).R2;
            end
        end
        fprintf('[QC] Loaded per-file R2 map from: %s (n=%d)\n', qc_mat, qcR2Map.Count);
    end
else
    fprintf('[QC] QC mat not found (OK): %s\n', qc_mat);
end

%% ---------------------- index DCM files and compute minimal metadata ----------------------
DCMidx = struct('subj',{},'cond',{},'run',{},'model',{},'file',{},'F',{},'meanR2',{},'skip',{},'skip_reason',{});
idx = 0;

for s = 1:numel(subject_list)
    subj = subject_list{s};
    subj_dcm_dir = fullfile(root_prf, subj, session_label, dcm_subdir);
    if ~isfolder(subj_dcm_dir)
        warning('[Index] Missing DCM dir, skip subject: %s', subj_dcm_dir);
        continue;
    end

    patt = sprintf('%s_%s_DCM_*_run*_spm25.mat', subj, session_label);
    files = dir(fullfile(subj_dcm_dir, patt));
    fprintf('[Index] %s: found %d DCM files in %s\n', subj, numel(files), subj_dcm_dir);

    for f = 1:numel(files)
        fn = files(f).name;
        full_fn = fullfile(files(f).folder, fn);

        % Parse: subj_ses-XX_DCM_<cond>_runNN_modelMM_spm25.mat
        expr = sprintf('%s_%s_DCM_(?<cond>.+?)_run(?<run>\\d+)(?:_model-?(?<model>\\d+))?_spm25\\.mat', subj, session_label);
        tok = regexp(fn, expr, 'names');
        if isempty(tok); continue; end

        condLabel = tok.cond;
        run_idx   = str2double(tok.run);
        if isfield(tok,'model') && ~isempty(tok.model)
            model_idx = str2double(tok.model);
        else
            model_idx = 1;
        end

        idx = idx + 1;
        DCMidx(idx).subj = subj;
        DCMidx(idx).cond = condLabel;
        DCMidx(idx).run  = run_idx;
        DCMidx(idx).model= model_idx;
        DCMidx(idx).file = full_fn;
        DCMidx(idx).F    = NaN;
        DCMidx(idx).meanR2 = NaN;
        DCMidx(idx).skip = false;
        DCMidx(idx).skip_reason = '';

        % Load DCM once to get F and (if needed) R2, and verify it is estimated.
        try
            S = load(full_fn, 'DCM');
            if ~isfield(S,'DCM'); error('No DCM variable'); end
            DCM = S.DCM;

            % Minimal validity check (not "QC", just required)
            if ~isfield(DCM,'Ep') || ~isfield(DCM,'Cp')
                DCMidx(idx).skip = true;
                DCMidx(idx).skip_reason = 'Invalid DCM: missing Ep/Cp (not estimated?)';
                continue;
            end

            if isfield(DCM,'F') && isfinite(DCM.F)
                DCMidx(idx).F = DCM.F;
            end

            % mean R2: prefer QC map if present; otherwise compute
            mR2 = NaN;
            key = normalize_path(full_fn);
            if qcR2Map.isKey(key)
                mR2 = qcR2Map(key);
                if isnan(mR2)
                    mR2 = compute_mean_R2_from_DCM(DCM);
                end
            else
                mR2 = compute_mean_R2_from_DCM(DCM);
            end
            DCMidx(idx).meanR2 = mR2;

            % Only exclusion rule requested: negative mean R2
            if ~isnan(mR2) && mR2 < 0
                DCMidx(idx).skip = true;
                DCMidx(idx).skip_reason = sprintf('QC: mean R2 < 0 (%.4f)', mR2);
            end

        catch ME
            DCMidx(idx).skip = true;
            DCMidx(idx).skip_reason = sprintf('Load/parse failed: %s', ME.message);
        end
    end
end

if isempty(DCMidx)
    error('No DCM files indexed. Check paths/naming.');
end

save(fullfile(out_dir,'metadata_DCM_index.mat'), 'DCMidx','root_prf','session_label','dcm_subdir','subject_list','modelNames','PP_THRESH','field');
fprintf('[Index] Saved DCM index: %s\n', fullfile(out_dir,'metadata_DCM_index.mat'));

%% ---------------------- Determine winning models (BMS) ----------------------
models_all = unique([DCMidx.model]);
models_all = models_all(:)';
nModel = max(models_all);

% Build log-evidence matrices for connect, unconnect, and combined
[lme_conn, subjs_conn] = build_lme_matrix(DCMidx, subject_list, 'connect', nModel);
[lme_unco, subjs_unco] = build_lme_matrix(DCMidx, subject_list, 'unconnect', nModel);
[lme_all , subjs_all ] = build_lme_matrix(DCMidx, subject_list, 'ALL', nModel);

win_model_connect   = choose_winning_model(lme_conn, 'connect', out_dir, modelNames);
win_model_unconnect = choose_winning_model(lme_unco, 'unconnect', out_dir, modelNames);
win_model_overall   = choose_winning_model(lme_all,  'combined', out_dir, modelNames);

fprintf('[BMS] Winning model (connect)   = %d (%s)\n', win_model_connect,   safe_model_name(win_model_connect, modelNames));
fprintf('[BMS] Winning model (unconnect) = %d (%s)\n', win_model_unconnect, safe_model_name(win_model_unconnect, modelNames));
fprintf('[BMS] Winning model (overall)   = %d (%s)\n', win_model_overall,   safe_model_name(win_model_overall, modelNames));

%% ---------------------- Run 3 analyses ----------------------
report_fn = fullfile(out_dir, sprintf('report_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
fid = fopen(report_fn,'w');
fprintf(fid, 'DCM 3-level PEB-of-PEBs (3 outputs) | QC: only skip mean R2 < 0 | PP>%.2f\n', PP_THRESH);
fprintf(fid, 'root_prf=%s\nsession=%s\ndcm_subdir=%s\nfield=%s\n\n', root_prf, session_label, dcm_subdir, strjoin(field,','));

% Section: QC skips
fprintf(fid, '==== QC skips (ONLY mean R2 < 0) and hard-invalid DCMs ====\n');
nSkip = 0;
for i = 1:numel(DCMidx)
    if DCMidx(i).skip
        nSkip = nSkip + 1;
        fprintf(fid, '[SKIP] %s | subj=%s cond=%s run=%02d model=%d | R2=%.4f | %s\n', ...
            DCMidx(i).file, DCMidx(i).subj, DCMidx(i).cond, DCMidx(i).run, DCMidx(i).model, ...
            DCMidx(i).meanR2, DCMidx(i).skip_reason);
    end
end
fprintf(fid, 'Total skipped files: %d / %d\n\n', nSkip, numel(DCMidx));

%% A) connect-only PEB3
fprintf(fid, '==== PEB3_connect (connect-only) ====\n');
res_connect = run_peb3_condition_only(DCMidx, subject_list, 'connect', win_model_connect, field, out_dir, PP_THRESH, modelNames);
write_result_summary(fid, res_connect);

%% B) unconnect-only PEB3
fprintf(fid, '\n==== PEB3_unconnect (unconnect-only) ====\n');
res_unconnect = run_peb3_condition_only(DCMidx, subject_list, 'unconnect', win_model_unconnect, field, out_dir, PP_THRESH, modelNames);
write_result_summary(fid, res_unconnect);

%% C) condEffect PEB3 (within-subject connect vs unconnect)
fprintf(fid, '\n==== PEB3_condEffect (within-subject condition effect: connect vs unconnect) ====\n');
res_cond = run_peb3_condition_effect(DCMidx, subject_list, win_model_overall, field, out_dir, PP_THRESH, modelNames);
write_result_summary(fid, res_cond);

fclose(fid);
fprintf('[Report] %s\n', report_fn);

end % main


%% =========================================================================
function res = run_peb3_condition_only(DCMidx, subject_list, condLabel, win_model, field, out_dir, PP_THRESH, modelNames)

% Build subject-level PEBs: within-subject over runs, intercept only
PEB2s = {};
PEB2_subj = {};
subj_used = {};
nDCM_used = [];
subj_fail = {};

for s = 1:numel(subject_list)
    subj = subject_list{s};

    files = select_dcm_files(DCMidx, subj, condLabel, win_model);
    if isempty(files); continue; end

    % X2: intercept only
    X2 = ones(numel(files),1);
    M2 = struct();
    M2.X = X2;
    M2.Xnames = {'mean'};
    M2.Q = 'all';

    try
        PEB2 = spm_dcm_peb(files, M2, field);
        PEB2s{end+1,1} = PEB2;
        PEB2_subj{end+1,1} = PEB2;
        subj_used{end+1,1} = subj;
        nDCM_used(end+1,1) = numel(files);

        save(fullfile(out_dir, sprintf('PEB2_%s_%s.mat', condLabel, subj)), 'PEB2', 'files', 'M2', 'field');

    catch ME
        subj_fail{end+1,1} = sprintf('%s: %s', subj, ME.message);
        continue;
    end
end

res = struct();
res.name = sprintf('PEB3_%s', condLabel);
res.cond = condLabel;
res.win_model = win_model;
res.win_model_name = safe_model_name(win_model, modelNames);
res.n_subjects = numel(subj_used);
res.subj_used = subj_used;
res.nDCM_used = nDCM_used;
res.subj_fail = subj_fail;

if res.n_subjects < 2
    res.ok = false;
    res.note = 'Not enough subjects for group PEB-of-PEBs (need >=2).';
    return;
end

% Level 3: group mean only
X3 = ones(res.n_subjects,1);
try
    PEB3 = spm_dcm_peb(PEB2s, X3);
    BMA3 = spm_dcm_peb_bmc(PEB3);

    res.ok = true;
    res.PEB3 = PEB3;
    res.BMA3 = BMA3;
    res.PP_THRESH = PP_THRESH;

    save(fullfile(out_dir, sprintf('%s.mat', res.name)), 'PEB3', 'X3', 'PEB2s', 'subj_used', 'nDCM_used', 'win_model', 'field');
    save(fullfile(out_dir, sprintf('BMA3_%s.mat', condLabel)), 'BMA3', 'PP_THRESH', 'win_model', 'condLabel', 'field');

catch ME
    res.ok = false;
    res.note = sprintf('Group PEB/BMA failed: %s', ME.message);
end

end


%% =========================================================================
function res = run_peb3_condition_effect(DCMidx, subject_list, win_model, field, out_dir, PP_THRESH, modelNames)

% Subject-level PEB: both conditions, X2=[1, cond], cond mean-centered per subject
PEB2s = {};
subj_used = {};
nDCM_used = [];
subj_fail = {};

for s = 1:numel(subject_list)
    subj = subject_list{s};

    files_conn = select_dcm_files(DCMidx, subj, 'connect', win_model);
    files_unco = select_dcm_files(DCMidx, subj, 'unconnect', win_model);

    if isempty(files_conn) || isempty(files_unco)
        % For estimating condition effect we require both conditions present
        continue;
    end

    files = [files_conn(:); files_unco(:)];

    condVec = [ ones(numel(files_conn),1); -ones(numel(files_unco),1) ];
    condVec = condVec - mean(condVec); % mean-center within subject

    X2 = [ones(numel(files),1), condVec];

    M2 = struct();
    M2.X = X2;
    M2.Xnames = {'mean','cond'};
    M2.Q = 'all';

    try
        PEB2 = spm_dcm_peb(files, M2, field);
        PEB2s{end+1,1} = PEB2;
        subj_used{end+1,1} = subj;
        nDCM_used(end+1,1) = numel(files);

        save(fullfile(out_dir, sprintf('PEB2_condEffect_%s.mat', subj)), 'PEB2', 'files', 'M2', 'field');

    catch ME
        subj_fail{end+1,1} = sprintf('%s: %s', subj, ME.message);
        continue;
    end
end

res = struct();
res.name = 'PEB3_condEffect';
res.cond = 'condEffect';
res.win_model = win_model;
res.win_model_name = safe_model_name(win_model, modelNames);
res.n_subjects = numel(subj_used);
res.subj_used = subj_used;
res.nDCM_used = nDCM_used;
res.subj_fail = subj_fail;

if res.n_subjects < 2
    res.ok = false;
    res.note = 'Not enough subjects with BOTH conditions for group PEB-of-PEBs (need >=2).';
    return;
end

% Level 3: group mean only (estimates group mean of "mean" and group mean of "cond")
X3 = ones(res.n_subjects,1);
try
    PEB3 = spm_dcm_peb(PEB2s, X3);
    BMA3 = spm_dcm_peb_bmc(PEB3);

    res.ok = true;
    res.PEB3 = PEB3;
    res.BMA3 = BMA3;
    res.PP_THRESH = PP_THRESH;

    save(fullfile(out_dir, sprintf('%s.mat', res.name)), 'PEB3', 'X3', 'PEB2s', 'subj_used', 'nDCM_used', 'win_model', 'field');
    save(fullfile(out_dir, 'BMA3_condEffect.mat'), 'BMA3', 'PP_THRESH', 'win_model', 'field');

catch ME
    res.ok = false;
    res.note = sprintf('Group PEB/BMA failed: %s', ME.message);
end

end


%% =========================================================================
function files = select_dcm_files(DCMidx, subj, condLabel, model_idx)
% Returns cell array of filenames for DCMs matching subj+cond+model and not skipped.
mask = strcmp({DCMidx.subj}, subj) & strcmp({DCMidx.cond}, condLabel) & ([DCMidx.model] == model_idx) & ~[DCMidx.skip];
files = {DCMidx(mask).file}';
end


%% =========================================================================
function [lme, subjs] = build_lme_matrix(DCMidx, subject_list, condLabel, nModel)
% lme: nSub x nModel, each entry is sum of DCM.F across all runs for that subject+cond+model (after skips)
nSub = numel(subject_list);
lme = nan(nSub, nModel);
subjs = subject_list(:);

for s = 1:nSub
    subj = subject_list{s};
    for m = 1:nModel
        if strcmpi(condLabel, 'ALL')
            mask = strcmp({DCMidx.subj}, subj) & ([DCMidx.model] == m) & ~[DCMidx.skip];
        else
            mask = strcmp({DCMidx.subj}, subj) & strcmp({DCMidx.cond}, condLabel) & ([DCMidx.model] == m) & ~[DCMidx.skip];
        end
        Fs = [DCMidx(mask).F];
        Fs = Fs(isfinite(Fs));
        if ~isempty(Fs)
            lme(s,m) = sum(Fs);
        end
    end
end
end


%% =========================================================================
function win_model = choose_winning_model(lme, tag, out_dir, modelNames)
% Choose winning model from log-evidence matrix (subjects x models).
% Tries RFX-BMS if enough complete rows; otherwise uses fixed-effect fallback.
win_model = NaN;

valid_rows = all(isfinite(lme),2);
nValid = sum(valid_rows);

if nValid >= 2 && exist('spm_BMS','file') == 2
    try
        [alpha,exp_r,xp,pxp,bor] = spm_BMS(lme(valid_rows,:));
        [~,win_model] = max(pxp);
        save(fullfile(out_dir, sprintf('BMS_%s.mat', tag)), 'lme','valid_rows','alpha','exp_r','xp','pxp','bor');
        return;
    catch
        % fallthrough to fallback
    end
end

% Fallback: fixed-effect on available evidence (nansum across subjects)
sumLME = nansum(lme,1);
[~,win_model] = max(sumLME);
save(fullfile(out_dir, sprintf('BMS_%s_FALLBACK.mat', tag)), 'lme','valid_rows','sumLME','win_model','modelNames');
end


%% =========================================================================
function write_result_summary(fid, res)
fprintf(fid, 'Winning model: %d (%s)\n', res.win_model, res.win_model_name);
fprintf(fid, 'Subjects included: %d\n', res.n_subjects);

if ~isempty(res.subj_used)
    for i = 1:numel(res.subj_used)
        fprintf(fid, '  - %s (nDCM=%d)\n', res.subj_used{i}, res.nDCM_used(i));
    end
end

if ~isempty(res.subj_fail)
    fprintf(fid, 'Subjects failed during PEB2:\n');
    for i = 1:numel(res.subj_fail)
        fprintf(fid, '  * %s\n', res.subj_fail{i});
    end
end

if ~res.ok
    fprintf(fid, 'STATUS: FAILED\nReason: %s\n', res.note);
    return;
end

fprintf(fid, 'STATUS: OK\n');

% Print PP>threshold parameters
try
    lines = summarize_highPP(res.BMA3, res.PP_THRESH, res.name);
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
catch
    fprintf(fid, '[Warn] Could not summarize PP>%.2f parameters (check BMA fields).\n', res.PP_THRESH);
end

end


%% =========================================================================
function lines = summarize_highPP(BMA, thr, tag)
% Summarize parameters with posterior probability > thr (exploratory PP threshold).
% Fixes common field-format issues:
%  - BMA.Pnames may be cell / string / char
%  - BMA.Ep may be sparse or non-vector; use spm_vec for robustness

lines = {};

if ~isstruct(BMA) || ~isfield(BMA,'Pp') || ~isfield(BMA,'Ep')
    lines{end+1} = sprintf('[%s] BMA missing Ep/Pp fields.', tag);
    return;
end

% --- Pp ---
try
    Pp = full(BMA.Pp(:));
catch
    Pp = BMA.Pp(:);
end

% --- Ep ---
try
    Ep = full(spm_vec(BMA.Ep));
catch
    try
        Ep = full(BMA.Ep(:));
    catch
        Ep = BMA.Ep(:);
    end
end
Ep = Ep(:);

% --- Pnames (robust) ---
Pnames = {};
if isfield(BMA,'Pnames') && ~isempty(BMA.Pnames)
    pn = BMA.Pnames;
    if iscell(pn)
        Pnames = pn(:);
    elseif isstring(pn)
        Pnames = cellstr(pn(:));
    elseif ischar(pn)
        Pnames = cellstr(pn);
    else
        try
            Pnames = cellstr(pn);
        catch
            Pnames = {};
        end
    end
end
if isempty(Pnames)
    Pnames = arrayfun(@(i) sprintf('param_%04d', i), 1:numel(Pp), 'UniformOutput', false)';
end

% Ensure each name is a char row
for i = 1:numel(Pnames)
    if isstring(Pnames{i}); Pnames{i} = char(Pnames{i}); end
end

% Align lengths defensively
n = min([numel(Pp), numel(Ep), numel(Pnames)]);
Pp = Pp(1:n);
Ep = Ep(1:n);
Pnames = Pnames(1:n);

% Optional: keep interpretation tidy
% - condEffect outputs typically want 'cond:' parameters
% - connect/unconnect outputs typically want 'mean:' parameters
prefix_mask = true(n,1);
try
    tag_l = lower(char(tag));
    is_cond = contains(tag_l, 'condeffect');
catch
    is_cond = ~isempty(strfind(lower(char(tag)), 'condeffect')); %#ok<STREMP>
end

% startsWith on cell can be version-dependent; use strncmp
has_cond = any(cellfun(@(s) ischar(s) && strncmp(s,'cond:',5), Pnames));
has_mean = any(cellfun(@(s) ischar(s) && strncmp(s,'mean:',5), Pnames));
if is_cond && has_cond
    prefix_mask = cellfun(@(s) ischar(s) && strncmp(s,'cond:',5), Pnames);
elseif ~is_cond && has_mean
    prefix_mask = cellfun(@(s) ischar(s) && strncmp(s,'mean:',5), Pnames);
end

% Threshold
mask = (Pp > thr) & isfinite(Ep) & prefix_mask;
idx  = find(mask);

lines{end+1} = sprintf('[%s] PP>%.2f parameters: %d / %d', tag, thr, numel(idx), sum(prefix_mask));
if isempty(idx)
    return;
end

% Sort by PP desc then |Ep| desc
[~,ord] = sortrows([Pp(idx), abs(Ep(idx))], [-1 -2]);
idx = idx(ord);

max_show = min(60, numel(idx));
for k = 1:max_show
    ii = idx(k);
    lines{end+1} = sprintf('  PP=%.3f | Ep=%+.4f | %s', Pp(ii), Ep(ii), Pnames{ii});
end
if numel(idx) > max_show
    lines{end+1} = sprintf('  ... (%d more above threshold)', numel(idx) - max_show);
end

end

Pp = BMA.Pp(:);
Ep = BMA.Ep(:);

if isfield(BMA,'Pnames')
    Pnames = cellstr(BMA.Pnames);
elseif isfield(BMA,'Pname')
    Pnames = cellstr(BMA.Pname);
else
    Pnames = arrayfun(@(i) sprintf('param_%04d', i), 1:numel(Pp), 'UniformOutput', false)';
end

mask = Pp > thr & isfinite(Ep);
idx = find(mask);

lines{end+1} = sprintf('[%s] PP>%.2f parameters: %d / %d', tag, thr, numel(idx), numel(Pp));
if isempty(idx); return; end

% Sort by PP desc then |Ep| desc
[~,ord] = sortrows([Pp(idx), abs(Ep(idx))], [-1 -2]);
idx = idx(ord);

max_show = min(60, numel(idx));
for k = 1:max_show
    i = idx(k);
    lines{end+1} = sprintf('  PP=%.3f | Ep=%+.4f | %s', Pp(i), Ep(i), Pnames{i});
end
if numel(idx) > max_show
    lines{end+1} = sprintf('  ... (%d more above threshold)', numel(idx) - max_show);
end

end


%% =========================================================================
function name = safe_model_name(m, modelNames)
if isnan(m) || m < 1 || m > numel(modelNames)
    name = 'unknown';
else
    name = modelNames{m};
end
end


%% =========================================================================
function p = normalize_path(p)
% Normalize path strings to improve matching across QC map
p = strrep(p, '\', '/');
end


%% =========================================================================
function meanR2 = compute_mean_R2_from_DCM(DCM)
% Prefer DCM.R2 if available; else compute from Y.y and DCM.y if possible.
meanR2 = NaN;
try
    if isfield(DCM,'R2') && ~isempty(DCM.R2)
        tmp = DCM.R2(:);
        tmp = tmp(~isnan(tmp));
        if ~isempty(tmp)
            meanR2 = mean(tmp);
            return;
        end
    end

    if isfield(DCM,'Y') && isfield(DCM.Y,'y') && isfield(DCM,'y') && ~isempty(DCM.Y.y) && ~isempty(DCM.y)
        Y_data = DCM.Y.y;
        Y_hat  = DCM.y;

        if iscell(Y_data); Y_data = cell2mat(Y_data); end
        if iscell(Y_hat);  Y_hat  = cell2mat(Y_hat);  end

        if ~isempty(Y_data) && ~isempty(Y_hat) && all(size(Y_data) == size(Y_hat))
            nReg = size(Y_data, 2);
            R2 = nan(1,nReg);
            for r = 1:nReg
                y  = Y_data(:,r);
                yh = Y_hat(:,r);
                if any(isnan(y)) || any(isnan(yh)); continue; end
                ymu = mean(y);
                SSE = sum((y - yh).^2);
                SST = sum((y - ymu).^2);
                if SST > 0
                    R2(r) = 1 - SSE./SST;
                end
            end
            tmp = R2(~isnan(R2));
            if ~isempty(tmp)
                meanR2 = mean(tmp);
            end
        end
    end
catch
    meanR2 = NaN;
end
end
