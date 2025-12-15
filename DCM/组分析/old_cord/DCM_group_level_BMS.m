%% DCM_group_level_BMS_spm25_compute.m
%  组水平 BMS + family-level BMS（unconnect）+ simple B 矩阵 +
%  统一 winning model 的 connection-wise 检验 +
%  每个条件使用各自 winning model 的 PEB/BMR（B 参数）
%  （本脚本只负责计算和保存，不画图）
%
%  适配 SPM25

clear; clc;

%% ====================== 路径 & 参数 ==========================
root_dir         = '/home/xue/data/prf_result';       % 数据根目录
analysis_dirname = 'DCM_manual_spm25_multiModel';     % 单被试 DCM 子目录
session_label    = 'ses-01';

% 被试列表
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09','sub-10'};
nSub     = numel(subjects);

% 条件列表
conds    = {'connect','unconnect'};

% 模型名字（顺序必须和 modelXX 的 XX 对齐）
model_names = {
    'B1_globalGain'   % model01
    'B2_NPOhub'       % model02
    'B3_feedforward'  % model03
    'B4_feedback'     % model04
    'B5_parietalOnly' % model05
    'B6_hierarchy'    % model06
    'B0_null'         % model07
};
N_MODEL = numel(model_names);

% 结果输出目录
out_dir = fullfile(root_dir, 'DCM_group_BMS');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('>>> root_dir          = %s\n', root_dir);
fprintf('>>> session_label     = %s\n', session_label);
fprintf('>>> analysis_dirname  = %s\n', analysis_dirname);
fprintf('>>> out_dir           = %s\n', out_dir);

% 检查 SPM 的 BMS 函数
if exist('spm_BMS','file') ~= 2
    warning(['>>> [Warn] 未检测到 spm_BMS，请先 addpath SPM12/SPM25 根目录，' ...
             '例如 addpath /path/to/spm12; spm(''defaults'',''fmri'');']);
else
    fprintf('>>> 已检测到 spm_BMS（将使用它做组水平 BMS）。\n');
end

%% ====================== 存放结果的结构体 ==========================
BMS_simple          = struct();
B_mean_all          = cell(1, numel(conds));   % 每个条件的 simple B_mean
best_model_idx_all  = nan(1, numel(conds));
B_diff_simple       = [];                      % simple 差异矩阵（unconnect - connect）

% 保存每个条件下「每个被试 × 每个模型」的最佳 DCM，用于后续统一模型分析 & PEB
DCMs_all = cell(1, numel(conds));

% 保存 PEB/BMR 结果
PEB_results = cell(1, numel(conds));   % 每个条件一个结构（PEB, BMA, PP 等）

%% ====================== 主循环：按条件跑 ==========================
for ic = 1:numel(conds)
    cond = conds{ic};
    fprintf('\n================= Condition: %s =================\n', cond);

    % DCMs{subIdx, modelIdx} 用来传给后续 BMS / 计算 B_mean / PEB
    DCMs  = cell(nSub, N_MODEL);
    % bestF(subIdx, modelIdx)：当前条件下，该被试该模型在所有 run 中的最大 F
    bestF = -Inf(nSub, N_MODEL);

    %% ---------- 逐被试扫描 DCM 目录 & 解析文件名 -----------------
    for is = 1:nSub
        subj = subjects{is};

        % DCM 目录：root/sub-XX/ses-01/DCM_manual_spm25_multiModel
        dcm_dir = fullfile(root_dir, subj, session_label, analysis_dirname);
        fprintf('  [Info] subj=%s, cond=%s\n', subj, cond);
        fprintf('        DCM 目录 = %s\n', dcm_dir);

        if ~exist(dcm_dir, 'dir')
            fprintf('        [Warn] 目录不存在 -> 这个被试这个条件直接跳过。\n');
            continue;
        end

        files = dir(fullfile(dcm_dir, '*.mat'));
        fprintf('        在该目录下找到 %d 个候选 .mat 文件。\n', numel(files));
        if isempty(files)
            continue;
        end

        % 文件名格式示例：
        % sub-03_ses-01_DCM_connect_run01_model01_spm25.mat
        for k = 1:numel(files)
            fname    = files(k).name;
            fullpath = fullfile(dcm_dir, fname);
            fprintf('        -> 解析文件: %s\n', fname);

            parts = strsplit(fname, '_');
            if numel(parts) < 7
                fprintf('           [Skip] 拆分后字段数 < 7 (只有 %d 段)，文件名格式不对，跳过。\n', numel(parts));
                continue;
            end

            subj_from_name = parts{1};  % sub-03
            ses_from_name  = parts{2};  % ses-01
            tag_from_name  = parts{3};  % DCM
            cond_from_name = parts{4};  % connect / unconnect
            run_str        = parts{5};  % run01
            model_str      = parts{6};  % model01

            % 只保留当前条件的 DCM
            if ~strcmpi(cond_from_name, cond)
                continue;
            end

            % sanity check（只是提示）
            if ~strcmpi(subj_from_name, subj)
                fprintf('           [Warn] 文件名里的 subj(%s) != 当前 subj(%s)\n', subj_from_name, subj);
            end
            if ~strcmpi(ses_from_name, session_label)
                fprintf('           [Warn] 文件名里的 ses(%s) != 当前 ses(%s)\n', ses_from_name, session_label);
            end
            if ~strcmpi(tag_from_name, 'DCM')
                fprintf('           [Warn] 第3段不是 "DCM" 而是 "%s"\n', tag_from_name);
            end

            % 解析 runXX 和 modelYY
            run_id   = NaN;
            model_id = NaN;

            if strncmp(run_str, 'run', 3)
                run_id = str2double(run_str(4:end));
            else
                fprintf('           [Warn] 第5段不是 run?? (%s)\n', run_str);
            end

            if strncmp(model_str, 'model', 5)
                model_id = str2double(model_str(6:end));
            else
                fprintf('           [Warn] 第6段不是 model?? (%s)\n', model_str);
            end

            if isnan(model_id) || model_id < 1 || model_id > N_MODEL
                fprintf('           [Skip] 无效模型编号 model_id=%g (N_MODEL=%d)，跳过。\n', model_id, N_MODEL);
                continue;
            end

            % 载入 DCM
            fprintf('           -> load %s\n', fullpath);
            S = load(fullpath, 'DCM');
            if ~isfield(S, 'DCM')
                fprintf('           [Skip] 文件中没有 DCM 变量，跳过。\n');
                continue;
            end
            DCM = S.DCM;

            % 读取 F（负自由能 ≈ log-evidence）
            if isfield(DCM, 'F') && isfinite(DCM.F)
                F = DCM.F;
            else
                fprintf('           [Warn] DCM 中缺少有效 F (log-evidence)，跳过此文件。\n');
                continue;
            end

            fprintf('           [Info] subj=%s, cond=%s, run=%g, model=%d(%s), F=%.3f\n', ...
                subj, cond, run_id, model_id, model_names{model_id}, F);

            % 同一 subj × model × cond，保留 F 最大的 run
            if F > bestF(is, model_id)
                bestF(is, model_id) = F;
                DCMs{is, model_id}  = DCM;
                fprintf('           [OK] 作为 subj=%s, model=%d(%s) 的当前最佳 DCM（run=%g, F=%.3f）。\n', ...
                    subj, model_id, model_names{model_id}, run_id, F);
            else
                fprintf('           [Info] 该 DCM 的 F=%.3f 未优于当前最佳 F=%.3f，忽略此 run。\n', ...
                    F, bestF(is, model_id));
            end
        end
    end

    %% ---------- 汇总：每个模型有多少有效被试 ---------------------
    for m = 1:N_MODEL
        n_valid = 0;
        for is = 1:nSub
            if ~isempty(DCMs{is, m})
                n_valid = n_valid + 1;
            end
        end
        fprintf('  [Summary] cond=%s, Model %d (%s): 有 %d 个有效 DCM。\n', ...
            cond, m, model_names{m}, n_valid);
    end

    %% 存一下这个条件下所有被试 × 模型的最佳 DCM，用于后续统一模型分析 & PEB
    DCMs_all{ic} = DCMs;

    %% ---------- 构造 log-evidence 矩阵并调用 spm_BMS -------------
    % lme: [nSub × N_MODEL]，每行为一个被试，每列一个模型
    lme = nan(nSub, N_MODEL);
    for is = 1:nSub
        for m = 1:N_MODEL
            if isempty(DCMs{is,m})
                continue;
            end
            D = DCMs{is,m};
            if isfield(D, 'F') && isfinite(D.F)
                lme(is,m) = D.F;
            else
                fprintf('  [Warn] subj=%s, model=%d(%s) 在 BMS 中缺少有效 F，将视为缺失。\n', ...
                    subjects{is}, m, model_names{m});
            end
        end
    end

    % 只保留“所有模型都有有效 F”的被试
    valid_rows = all(~isnan(lme), 2);
    n_validSub = sum(valid_rows);

    if n_validSub < 2
        fprintf('  [Summary] Condition %s: 有效被试(<2)不足以做 BMS，跳过。\n', cond);

        BMS_simple(ic).cond_name        = cond;
        BMS_simple(ic).posterior        = [];
        BMS_simple(ic).out              = [];
        BMS_simple(ic).best_model_idx   = NaN;
        BMS_simple(ic).best_model_name  = '';
        BMS_simple(ic).exp_r            = [];
        BMS_simple(ic).pxp              = [];
        BMS_simple(ic).B_mean           = [];
        BMS_simple(ic).subjects_in_BMS  = {};
        BMS_simple(ic).family           = [];
        BMS_simple(ic).model_rank_idx   = [];
        BMS_simple(ic).model_rank_pxp   = [];
        BMS_simple(ic).model_rank_exp_r = [];

        B_mean_all{ic}         = [];
        best_model_idx_all(ic) = NaN;
        continue;
    end

    fprintf('  [Info] Condition %s: 参与 BMS 的被试数 = %d\n', cond, n_validSub);
    lme_valid   = lme(valid_rows, :);
    subj_in_BMS = subjects(valid_rows);

    % 调用 spm_BMS 做 RFX-BMS（保护的 exceedance probability）
    posterior = [];
    out       = [];
    if exist('spm_BMS','file') ~= 2
        warning('  [Error] spm_BMS 不在路径中，无法执行 BMS。');
    else
        fprintf('  >>> 开始调用 spm_BMS 进行组水平 RFX-BMS...\n');
        try
            nsamps = 1e5;  % Monte Carlo 采样次数
            [alpha, exp_r, xp, pxp, bor] = spm_BMS(lme_valid, nsamps, 0, 0, 1);

            posterior = struct();
            posterior.alpha    = alpha;
            posterior.r        = exp_r(:);  % expected model frequencies
            posterior.xp       = xp(:);     % exceedance probabilities
            posterior.pxp      = pxp(:);    % protected exceedance probabilities
            posterior.bor      = bor;       % Bayes Omnibus Risk
            posterior.lme      = lme_valid;
            posterior.subjects = subj_in_BMS(:);

            out = struct();
            out.method      = 'RFX_spm_BMS';
            out.nsamps      = nsamps;
            out.lme         = lme_valid;
            out.subjects    = subj_in_BMS(:);
            out.model_names = model_names(:);

        catch ME
            warning('  [Error] 调用 spm_BMS 出错: %s', ME.message);
            posterior = [];
            out       = [];
        end
    end

    if isempty(posterior)
        fprintf('  [Summary] Condition %s: BMS 失败，posterior 为空。\n', cond);

        BMS_simple(ic).cond_name        = cond;
        BMS_simple(ic).posterior        = [];
        BMS_simple(ic).out              = [];
        BMS_simple(ic).best_model_idx   = NaN;
        BMS_simple(ic).best_model_name  = '';
        BMS_simple(ic).exp_r            = [];
        BMS_simple(ic).pxp              = [];
        BMS_simple(ic).B_mean           = [];
        BMS_simple(ic).subjects_in_BMS  = subj_in_BMS;
        BMS_simple(ic).family           = [];
        BMS_simple(ic).model_rank_idx   = [];
        BMS_simple(ic).model_rank_pxp   = [];
        BMS_simple(ic).model_rank_exp_r = [];

        B_mean_all{ic}         = [];
        best_model_idx_all(ic) = NaN;
        continue;
    end

    % 打印 expected model frequency
    exp_r = posterior.r(:);
    fprintf('  [BMS] Expected model frequencies (exp_r):\n');
    for m = 1:N_MODEL
        fprintf('         Model %d (%s): exp_r = %.4f\n', m, model_names{m}, exp_r(m));
    end

    % 打印 protected exceedance probability，并按 pxp 完整排序 7 个模型
    if isfield(posterior, 'pxp') && ~isempty(posterior.pxp)
        pxp = posterior.pxp(:);
        fprintf('  [BMS] Protected exceedance probability (pxp):\n');
        for m = 1:N_MODEL
            fprintf('         Model %d (%s): pxp = %.4f\n', m, model_names{m}, pxp(m));
        end

        [sorted_pxp, sort_idx] = sort(pxp, 'descend');
        fprintf('  [BMS] 模型按 pxp 从高到低完整排序 (cond = %s):\n', cond);
        for ii = 1:N_MODEL
            midx = sort_idx(ii);
            fprintf('         #%d: Model %d (%s): pxp = %.4f, exp_r = %.4f\n', ...
                ii, midx, model_names{midx}, sorted_pxp(ii), exp_r(midx));
        end

        % 保存排序信息到 posterior 结构
        posterior.sort_idx_by_pxp = sort_idx;
        posterior.sorted_pxp      = sorted_pxp;
        posterior.sorted_exp_r    = exp_r(sort_idx);
    else
        pxp = [];
        fprintf('  [Warn] posterior 中没有 pxp 字段或为空。\n');
    end

    % 最优模型：优先按 pxp 最大，其次按 exp_r 最大
    if ~isempty(pxp)
        [~, best_idx] = max(pxp);
        crit_str = 'pxp';
    else
        [~, best_idx] = max(exp_r);
        crit_str = 'exp_r';
    end
    best_name     = model_names{best_idx};
    fprintf('  [Result] Condition %s: Best model = %d (%s) (criterion = %s)\n', ...
        cond, best_idx, best_name, crit_str);

    %% ---------- 仅 unconnect 条件下做 family-level BMS ----------
    % Family 划分：
    %   Family 1: B1_globalGain, B2_NPOhub, B5_parietalOnly    -> 模型 1,2,5
    %   Family 2: B3_feedforward, B4_feedback, B6_hierarchy    -> 模型 3,4,6
    %   Family 3: B0_null                                      -> 模型 7
    family_struct = [];
    if strcmp(cond, 'unconnect')
        fprintf('  [Family-BMS] 对 unconnect 条件进行 family-level BMS...\n');

        if exist('spm_BMS','file') ~= 2
            fprintf('  [Family-BMS] spm_BMS 不在路径中，无法做 family-level BMS。\n');
            family_struct = [];
        else
            fam_names = { ...
                'Gain-like (B1,B2,B5)', ...
                'Directional/hierarchy (B3,B4,B6)', ...
                'Null (B0)'};
            % 每个模型所属 family 的编号（长度 = N_MODEL = 7）
            %   model1->1, model2->1, model3->2, model4->2, model5->1, model6->2, model7->3
            fam_index = [1 1 2 2 1 2 3];
            nFam      = numel(fam_names);

            % 打印 family 划分
            fprintf('  [Family-BMS] Family 划分如下：\n');
            for f = 1:nFam
                idx_models_f = find(fam_index == f);
                fprintf('         Family %d (%s): 包含模型 -> ', f, fam_names{f});
                for mm = idx_models_f
                    fprintf('%d(%s) ', mm, model_names{mm});
                end
                fprintf('\n');
            end

            % 计算每个 family 的 log-evidence（对 family 内模型做 log-sum-exp）
            lme_fam = nan(n_validSub, nFam);   % [n_validSub × nFam]
            for f = 1:nFam
                idx_models_f = find(fam_index == f);
                if isempty(idx_models_f)
                    continue;
                end
                evid_f = lme_valid(:, idx_models_f);   % [n_validSub × nModelsInFamily]
                for is2 = 1:n_validSub
                    row = evid_f(is2, :);
                    row = row(~isnan(row));
                    if isempty(row)
                        lme_fam(is2, f) = NaN;
                    else
                        mx = max(row);
                        lme_fam(is2, f) = mx + log(sum(exp(row - mx)));
                    end
                end
            end

            fam_valid_rows = all(~isnan(lme_fam), 2);
            n_validSub_fam = sum(fam_valid_rows);
            if n_validSub_fam < 2
                fprintf('  [Family-BMS] 有效被试(<2)不足，跳过 family-level BMS。\n');
                family_struct = [];
            else
                fprintf('  [Family-BMS] 参与 family-level BMS 的被试数 = %d\n', n_validSub_fam);
                try
                    nsamps_f = 1e5;
                    [alpha_f, exp_r_f, xp_f, pxp_f, bor_f] = spm_BMS(lme_fam(fam_valid_rows, :), nsamps_f, 0, 0, 1);

                    fprintf('  [Family-BMS] Expected family frequencies (exp_r_f):\n');
                    for f = 1:nFam
                        fprintf('         Family %d (%s): exp_r = %.4f\n', ...
                            f, fam_names{f}, exp_r_f(f));
                    end
                    fprintf('  [Family-BMS] Protected exceedance probability (pxp_f):\n');
                    for f = 1:nFam
                        fprintf('         Family %d (%s): pxp = %.4f\n', ...
                            f, fam_names{f}, pxp_f(f));
                    end

                    % 按 pxp_f 排序输出 summary
                    [~, fam_sort_idx] = sort(pxp_f, 'descend');
                    fprintf('  [Family-BMS] Families sorted by pxp_f (from high to low):\n');
                    for ii = 1:nFam
                        f = fam_sort_idx(ii);
                        fprintf('         #%d: Family %d (%s): pxp = %.4f, exp_r = %.4f\n', ...
                            ii, f, fam_names{f}, pxp_f(f), exp_r_f(f));
                    end

                    [~, best_f] = max(pxp_f);
                    fprintf('  [Family-BMS] Best family (by pxp_f) = %d (%s)\n', ...
                        best_f, fam_names{best_f});

                    % 在 winning family 内再做一次模型 BMS
                    idx_models_bestFam = find(fam_index == best_f);
                    fprintf('  [Family-BMS] 在 winning family 内再做一次模型级 BMS，包含模型: ');
                    for mm = idx_models_bestFam
                        fprintf('%d(%s) ', mm, model_names{mm});
                    end
                    fprintf('\n');

                    lme_bestFam = lme_valid(:, idx_models_bestFam);
                    % 再次保证每个被试对这些模型都有值
                    valid_rows_m = all(~isnan(lme_bestFam), 2);
                    n_validSub_m = sum(valid_rows_m);
                    if n_validSub_m < 2
                        fprintf('  [Family-BMS] winning family 内有效被试不足(<2)，跳过模型级 BMS。\n');
                        posterior_bestFam = [];
                    else
                        [alpha_m, exp_r_m, xp_m, pxp_m, bor_m] = spm_BMS(lme_bestFam(valid_rows_m,:), nsamps_f, 0, 0, 1);
                        fprintf('  [Family-BMS] winning family 内模型级 BMS 结果：\n');
                        for ii = 1:numel(idx_models_bestFam)
                            midx = idx_models_bestFam(ii);
                            fprintf('         Model %d (%s): exp_r = %.4f, pxp = %.4f\n', ...
                                midx, model_names{midx}, exp_r_m(ii), pxp_m(ii));
                        end
                        [~, best_idx_local] = max(exp_r_m);
                        best_model_in_family = idx_models_bestFam(best_idx_local);
                        fprintf('  [Family-BMS] winning family 内最优模型 = %d (%s)\n', ...
                            best_model_in_family, model_names{best_model_in_family});

                        posterior_bestFam = struct();
                        posterior_bestFam.alpha = alpha_m;
                        posterior_bestFam.r     = exp_r_m(:);
                        posterior_bestFam.xp    = xp_m(:);
                        posterior_bestFam.pxp   = pxp_m(:);
                        posterior_bestFam.bor   = bor_m;
                        posterior_bestFam.models = idx_models_bestFam(:);
                    end

                    fam_post = struct();
                    fam_post.alpha    = alpha_f;
                    fam_post.r        = exp_r_f(:);
                    fam_post.xp       = xp_f(:);
                    fam_post.pxp      = pxp_f(:);
                    fam_post.bor      = bor_f;
                    fam_post.lme      = lme_fam;
                    fam_post.subjects = subj_in_BMS(fam_valid_rows);

                    family_struct = struct();
                    family_struct.names               = fam_names(:);
                    family_struct.index_per_model     = fam_index(:);
                    family_struct.posterior           = fam_post;
                    family_struct.best_family_index   = best_f;
                    family_struct.best_family_name    = fam_names{best_f};
                    family_struct.BMS_within_family   = posterior_bestFam;

                catch ME
                    warning('  [Family-BMS] spm_BMS (family-level) 出错: %s', ME.message);
                    family_struct = [];
                end
            end
        end
    end

    %% ---------- 计算该条件下最佳模型的 B 矩阵 simple 平均 ---------------
    B_mean = [];
    fprintf('  >>> 计算 cond=%s, best_model=%d(%s) 的 B 矩阵简单平均...\n', ...
        cond, best_idx, best_name);

    nB = 0;
    for is = 1:nSub
        if isempty(DCMs{is, best_idx})
            continue;
        end

        D = DCMs{is, best_idx};

        % 优先 Ep.B，否则 D.B
        if isfield(D, 'Ep') && isfield(D.Ep, 'B') && ~isempty(D.Ep.B)
            B = D.Ep.B;
        elseif isfield(D, 'B') && ~isempty(D.B)
            B = D.B;
        else
            fprintf('      [Warn] subj=%s: DCM 中找不到 Ep.B 或 B 字段，跳过。\n', subjects{is});
            continue;
        end

        % 如果 B 是 cell（例如多调制），取第一个
        if iscell(B)
            if ~isempty(B{1})
                B = B{1};
            else
                fprintf('      [Warn] subj=%s: B 是 cell 但第一项为空，跳过。\n', subjects{is});
                continue;
            end
        end

        if isempty(B_mean)
            B_mean = zeros(size(B));
        end

        if ~isequal(size(B), size(B_mean))
            fprintf('      [Warn] subj=%s: B 尺寸(%s) 与当前平均尺寸(%s) 不一致，跳过该被试。\n', ...
                subjects{is}, mat2str(size(B)), mat2str(size(B_mean)));
            continue;
        end

        B_mean = B_mean + B;
        nB     = nB + 1;
    end

    if nB > 0
        B_mean = B_mean ./ nB;
        fprintf('  [OK] cond=%s: simple B_mean 计算完成 (聚合了 %d 个被试)。\n', cond, nB);
        fprintf('       B_mean 尺寸 = %s\n', mat2str(size(B_mean)));
    else
        B_mean = [];
        fprintf('  [Warn] cond=%s: 没有任何被试贡献 B 矩阵，B_mean 为空。\n', cond);
    end

    % 结果写入结构体
    BMS_simple(ic).cond_name        = cond;
    BMS_simple(ic).posterior        = posterior;
    BMS_simple(ic).out              = out;
    BMS_simple(ic).best_model_idx   = best_idx;
    BMS_simple(ic).best_model_name  = best_name;
    BMS_simple(ic).exp_r            = exp_r;
    BMS_simple(ic).pxp              = pxp;
    BMS_simple(ic).B_mean           = B_mean;
    BMS_simple(ic).subjects_in_BMS  = subj_in_BMS;
    BMS_simple(ic).family           = family_struct;
    if ~isempty(posterior) && isfield(posterior,'sort_idx_by_pxp')
        BMS_simple(ic).model_rank_idx   = posterior.sort_idx_by_pxp(:);
        BMS_simple(ic).model_rank_pxp   = posterior.sorted_pxp(:);
        BMS_simple(ic).model_rank_exp_r = posterior.sorted_exp_r(:);
    else
        BMS_simple(ic).model_rank_idx   = [];
        BMS_simple(ic).model_rank_pxp   = [];
        BMS_simple(ic).model_rank_exp_r = [];
    end

    B_mean_all{ic}         = B_mean;
    best_model_idx_all(ic) = best_idx;
end

%% ====================== 判定各条件 winning model 是否“独占鳌头” ==================
fprintf('\n================= Winning model dominance check =================\n');
% 判定规则：pxp(best) >= dom_thr_pxp 且 pxp(best) - pxp(second) >= gap_thr_pxp
dom_thr_pxp = 0.50;
gap_thr_pxp = 0.10;
fprintf('  判定规则: pxp(best) >= %.2f 且 pxp(best) - pxp(second) >= %.2f 视为“独占鳌头”。\n', ...
    dom_thr_pxp, gap_thr_pxp);

winning_model_info = struct('cond_name',{}, 'best_model_idx',{}, 'best_model_name',{}, ...
                            'best_pxp',{}, 'second_best_idx',{}, 'second_best_pxp',{}, ...
                            'is_dominant',{});
is_dominant = false(1, numel(conds));

for ic = 1:numel(conds)
    cond = conds{ic};
    info = struct();
    info.cond_name       = cond;
    info.best_model_idx  = best_model_idx_all(ic);

    if isnan(info.best_model_idx) || isempty(BMS_simple(ic).pxp)
        fprintf('  [Dominance] cond=%s: 缺少 pxp 或 best_model_idx 为 NaN，无法判定。\n', cond);
        info.best_model_name   = '';
        info.best_pxp          = NaN;
        info.second_best_idx   = NaN;
        info.second_best_pxp   = NaN;
        info.is_dominant       = false;
    else
        info.best_model_name   = model_names{info.best_model_idx};
        pxp = BMS_simple(ic).pxp(:);

        if numel(pxp) ~= N_MODEL
            fprintf('  [Dominance] cond=%s: pxp 长度(%d) != N_MODEL(%d)，无法判定。\n', ...
                cond, numel(pxp), N_MODEL);
            info.best_pxp        = NaN;
            info.second_best_idx = NaN;
            info.second_best_pxp = NaN;
            info.is_dominant     = false;
        else
            [sorted_pxp, sort_idx] = sort(pxp, 'descend');
            info.best_pxp = pxp(info.best_model_idx);
            if numel(sorted_pxp) >= 2
                info.second_best_idx = sort_idx(2);
                info.second_best_pxp = sorted_pxp(2);
            else
                info.second_best_idx = NaN;
                info.second_best_pxp = NaN;
            end

            info.is_dominant = ~isnan(info.best_pxp) && ~isnan(info.second_best_pxp) && ...
                               (info.best_pxp >= dom_thr_pxp) && ...
                               ((info.best_pxp - info.second_best_pxp) >= gap_thr_pxp);

            fprintf('  [Dominance] cond=%s: best=%d(%s), pxp(best)=%.3f, second=%d, pxp(second)=%.3f, dominant=%d\n', ...
                cond, info.best_model_idx, info.best_model_name, info.best_pxp, ...
                info.second_best_idx, info.second_best_pxp, info.is_dominant);
        end
    end

    winning_model_info(ic) = info;
    is_dominant(ic)        = info.is_dominant;
end

%% ====================== simple 差异矩阵（不画图） ==================
idx_connect   = find(strcmp(conds, 'connect'));
idx_unconnect = find(strcmp(conds, 'unconnect'));

if ~isempty(idx_connect) && ~isempty(idx_unconnect)
    Bc = B_mean_all{idx_connect};
    Bu = B_mean_all{idx_unconnect};

    if ~isempty(Bc) && ~isempty(Bu) && isequal(size(Bc), size(Bu))
        B_diff_simple = Bu - Bc;
        fprintf('\n[Diff] 已计算 simple B_mean(unconnect) - B_mean(connect)，尺寸 = %s\n', ...
            mat2str(size(B_diff_simple)));
    else
        fprintf('\n[Diff] connect/unconnect 的 simple B_mean 尺寸不一致或为空，无法计算差异。\n');
        B_diff_simple = [];
    end
else
    fprintf('\n[Diff] 未能在 conds 中同时找到 connect 和 unconnect 条件。\n');
    B_diff_simple = [];
end

%% ====================== 统一 winning model 的 connection-wise 检验 ==================
fprintf('\n================= Connection-wise tests (统一 winning model) =================\n');

connwise_results = struct([]);  % 可能有多个 unified model，对应多个结果块

if ~isempty(idx_connect) && ~isempty(idx_unconnect)
    DCMs_conn = DCMs_all{idx_connect};
    DCMs_unco = DCMs_all{idx_unconnect};

    if isempty(DCMs_conn) || isempty(DCMs_unco)
        fprintf('[ConnWise] connect 或 unconnect 条件下的 DCMs_all 为空，跳过统一模型分析。\n');
    else
        idx_dom = find(is_dominant);
        if isempty(idx_dom)
            fprintf('[ConnWise] 两个条件均无“独占鳌头”的 winning model，根据设定不进行 connection-wise tests。\n');
        else
            % 原始的“来源条件索引”和对应的 winning model 编号
            raw_models = best_model_idx_all(idx_dom);
            [unify_models_unique, ia_unique] = unique(raw_models, 'stable'); %#ok<NASGU>
            nJobs = numel(unify_models_unique);

            fprintf('[ConnWise] 将对 %d 个统一模型进行 connection-wise tests。\n', nJobs);
            job_count = 0;

            for uj = 1:nJobs
                unified_model_idx = unify_models_unique(uj);
                unified_model_name = model_names{unified_model_idx};

                % 找出哪些条件把这个模型视为 dominant winner
                src_cond_idx = idx_dom(raw_models == unified_model_idx);
                src_cond_names = conds(src_cond_idx);

                % 把来源条件名字拼成一个字符串
                src_str = '';
                for kk = 1:numel(src_cond_names)
                    if kk == 1
                        src_str = src_cond_names{kk};
                    else
                        src_str = [src_str ',' src_cond_names{kk}]; %#ok<AGROW>
                    end
                end

                job_count = job_count + 1;
                fprintf('\n[ConnWise] ---- Job #%d: Unified model = %d (%s), 来源条件 = %s ----\n', ...
                    job_count, unified_model_idx, unified_model_name, src_str);

                % 准备收集 B 矩阵
                B_conn_all   = [];
                B_unconn_all = [];
                subj_used    = {};
                used_count   = 0;

                % 先从一个 DCM 里读 ROI label（用于打印；如果没有 xY，退化成 ROI1,ROI2,...）
                roi_labels = {};
                DCM_example = [];
                for is = 1:nSub
                    if ~isempty(DCMs_conn{is, unified_model_idx})
                        DCM_example = DCMs_conn{is, unified_model_idx};
                        break;
                    end
                end
                if ~isempty(DCM_example) && isfield(DCM_example, 'xY')
                    try
                        nx = numel(DCM_example.xY);
                        roi_labels = cell(nx,1);
                        for ii = 1:nx
                            nm = DCM_example.xY(ii).name;  % 例如 '1 LNF'
                            if ischar(nm)
                                tok = regexp(nm, '^\s*\d+\s+(.*)$', 'tokens', 'once');
                                if ~isempty(tok)
                                    roi_labels{ii} = tok{1};
                                else
                                    roi_labels{ii} = nm;
                                end
                            else
                                roi_labels{ii} = sprintf('ROI%d', ii);
                            end
                        end
                    catch
                        roi_labels = {};
                    end
                end

                % 收集每个被试在统一模型下的 B_connect 和 B_unconnect
                for is = 1:nSub
                    D_c = DCMs_conn{is, unified_model_idx};
                    D_u = DCMs_unco{is, unified_model_idx};
                    if isempty(D_c) || isempty(D_u)
                        continue;
                    end

                    % 提取 connect 的 B
                    if isfield(D_c, 'Ep') && isfield(D_c.Ep, 'B') && ~isempty(D_c.Ep.B)
                        Bc = D_c.Ep.B;
                    elseif isfield(D_c, 'B') && ~isempty(D_c.B)
                        Bc = D_c.B;
                    else
                        continue;
                    end
                    if iscell(Bc)
                        if ~isempty(Bc{1})
                            Bc = Bc{1};
                        else
                            continue;
                        end
                    end

                    % 提取 unconnect 的 B
                    if isfield(D_u, 'Ep') && isfield(D_u.Ep, 'B') && ~isempty(D_u.Ep.B)
                        Bu = D_u.Ep.B;
                    elseif isfield(D_u, 'B') && ~isempty(D_u.B)
                        Bu = D_u.B;
                    else
                        continue;
                    end
                    if iscell(Bu)
                        if ~isempty(Bu{1})
                            Bu = Bu{1};
                        else
                            continue;
                        end
                    end

                    % 初始化 3D 数组
                    if isempty(B_conn_all)
                        [nR, nC] = size(Bc);
                        B_conn_all   = zeros(nR, nC, nSub);
                        B_unconn_all = zeros(nR, nC, nSub);
                        subj_used    = cell(nSub,1);
                        used_count   = 0;
                    end

                    if ~isequal(size(Bc), [nR, nC]) || ~isequal(size(Bu), [nR, nC])
                        fprintf('[ConnWise] subj=%s: B 尺寸不一致，跳过。\n', subjects{is});
                        continue;
                    end

                    used_count = used_count + 1;
                    B_conn_all(:,:,used_count)   = Bc;
                    B_unconn_all(:,:,used_count) = Bu;
                    subj_used{used_count}        = subjects{is};
                end

                if used_count < 2
                    fprintf('[ConnWise] Job #%d: 统一模型下同时有 connect/unconnect DCM 的被试不足(<2)，跳过。\n', job_count);
                    continue;
                end

                B_conn_all   = B_conn_all(:,:,1:used_count);
                B_unconn_all = B_unconn_all(:,:,1:used_count);
                subj_used    = subj_used(1:used_count);
                fprintf('[ConnWise] Job #%d: 统一模型下有效被试数 = %d\n', job_count, used_count);

                nR = size(B_conn_all,1);

                % 为每条边准备 p 值矩阵
                p_conn_vs0   = nan(nR, nR);
                p_uncon_vs0  = nan(nR, nR);
                p_diff       = nan(nR, nR);
                mean_conn    = nan(nR, nR);
                mean_uncon   = nan(nR, nR);
                mean_diff    = nan(nR, nR);

                % 对角线不分析
                for i = 1:nR
                    for j = 1:nR
                        if i == j
                            continue;
                        end

                        x_conn = squeeze(B_conn_all(i,j,:));
                        x_unco = squeeze(B_unconn_all(i,j,:));
                        x_diff = x_unco - x_conn;

                        mean_conn(i,j)  = mean(x_conn);
                        mean_uncon(i,j) = mean(x_unco);
                        mean_diff(i,j)  = mean(x_diff);

                        % connect vs 0
                        p_conn_vs0(i,j) = local_one_sample_p(x_conn);

                        % unconnect vs 0
                        p_uncon_vs0(i,j) = local_one_sample_p(x_unco);

                        % paired diff vs 0
                        p_diff(i,j) = local_one_sample_p(x_diff);
                    end
                end

                % FDR 校正（针对差值 p_diff）
                alpha_FDR = 0.05;
                mask_offdiag = true(nR, nR);
                mask_offdiag(1:nR+1:end) = false; % 去掉对角线
                p_vec  = p_diff(mask_offdiag);
                idx_ok = ~isnan(p_vec);
                p_vec_valid = p_vec(idx_ok);

                sig_mask_diff = false(nR, nR);
                if ~isempty(p_vec_valid)
                    [p_thresh, sig_vec_flat] = local_fdr_bh(p_vec_valid, alpha_FDR);
                    if ~isnan(p_thresh)
                        tmp_mask = false(size(p_vec));
                        tmp_mask(idx_ok) = sig_vec_flat;
                        sig_mask_diff(mask_offdiag) = tmp_mask;
                    end
                    fprintf('[ConnWise] Job #%d: FDR BH 阈值 p_thresh = %.4g (alpha=%.3f)\n', ...
                        job_count, p_thresh, alpha_FDR);
                else
                    fprintf('[ConnWise] Job #%d: 差值 p_diff 全为 NaN，无法做 FDR。\n');
                end

                % 打印差值显著的边
                fprintf('[ConnWise] Job #%d: 差值 (B_unconnect - B_connect) 经 BH-FDR 校正后显著的边：\n', job_count);
                any_sig = false;
                for i = 1:nR
                    for j = 1:nR
                        if i == j
                            continue;
                        end
                        if sig_mask_diff(i,j)
                            any_sig = true;
                            if ~isempty(roi_labels) && numel(roi_labels) >= nR
                                from_name = roi_labels{i};
                                to_name   = roi_labels{j};
                            else
                                from_name = sprintf('ROI%d', i);
                                to_name   = sprintf('ROI%d', j);
                            end
                            fprintf('   %s -> %s: mean_conn=%.4f, mean_uncon=%.4f, diff=%.4f, p_diff=%.4g\n', ...
                                from_name, to_name, mean_conn(i,j), mean_uncon(i,j), mean_diff(i,j), p_diff(i,j));
                        end
                    end
                end
                if ~any_sig
                    fprintf('   Job #%d: 无任何边在差值层面通过 FDR 校正。\n', job_count);
                end

                % 保存 connection-wise 结果
                connwise_results(job_count).unified_model_idx   = unified_model_idx;
                connwise_results(job_count).unified_model_name  = unified_model_name;
                connwise_results(job_count).source_cond_index   = src_cond_idx;
                connwise_results(job_count).source_cond_names   = src_cond_names;
                connwise_results(job_count).subjects            = subj_used;
                connwise_results(job_count).B_conn_all          = B_conn_all;
                connwise_results(job_count).B_unconn_all        = B_unconn_all;
                connwise_results(job_count).mean_conn           = mean_conn;
                connwise_results(job_count).mean_uncon          = mean_uncon;
                connwise_results(job_count).mean_diff           = mean_diff;
                connwise_results(job_count).p_conn_vs0          = p_conn_vs0;
                connwise_results(job_count).p_uncon_vs0         = p_uncon_vs0;
                connwise_results(job_count).p_diff              = p_diff;
                connwise_results(job_count).sig_mask_diff_FDR   = sig_mask_diff;
                connwise_results(job_count).roi_labels          = roi_labels;
            end
        end
    end
else
    fprintf('[ConnWise] 缺少 connect/unconnect 任一条件，跳过统一模型分析。\n');
end

%% ====================== PEB + BMR（对 B 参数，SPM25 版，每个条件用自己的 winning model） ==================
fprintf('\n================= PEB + BMR on B parameters (per condition, own winning model) =================\n');

PEB_results = cell(1, numel(conds));   % 重置 / 覆盖

for ic = 1:numel(conds)
    cond = conds{ic};
    fprintf('\n[PEB] Condition: %s\n', cond);

    % ---------- 使用该条件自己的 winning model ----------
    this_best_idx = best_model_idx_all(ic);
    if isnan(this_best_idx)
        fprintf('  [PEB] 该条件没有有效的 winning model，跳过。\n');
        PEB_results{ic} = [];
        continue;
    end
    unified_model_idx  = this_best_idx;
    unified_model_name = model_names{unified_model_idx};
    fprintf('  [PEB] 使用该条件的 winning model %d (%s) 进行 PEB。\n', ...
        unified_model_idx, unified_model_name);

    DCMs_cond = DCMs_all{ic};
    if isempty(DCMs_cond)
        fprintf('  [PEB] 该条件下 DCMs_all 为空。\n');
        PEB_results{ic} = [];
        continue;
    end

    % ---------- 收集该条件下统一模型的 DCM 列表（GCM） ----------
    GCM       = {};
    subj_list = {};
    for is = 1:nSub
        D = DCMs_cond{is, unified_model_idx};
        if isempty(D)
            continue;
        end
        GCM{end+1,1}        = D;            %#ok<AGROW>
        subj_list{end+1,1}  = subjects{is}; %#ok<AGROW>
    end

    n_condSub = numel(GCM);
    if n_condSub < 2
        fprintf('  [PEB] 该条件下统一模型的有效被试不足(<2)，跳过 PEB。\n');
        PEB_results{ic} = [];
        continue;
    end
    fprintf('  [PEB] 该条件下统一模型有效被试数 = %d\n', n_condSub);

    % ---------- 群组设计矩阵（group mean） ----------
    M        = struct();
    M.Q      = 'all';
    M.X      = ones(n_condSub, 1);
    M.Xnames = {'GroupMean'};

    % ---------- 调用 spm_dcm_peb & spm_dcm_peb_bmc ----------
    PEB_cond = [];
    BMA_cond = [];
    BMR_cond = [];   % 这里不单独使用 BMR，先留空占位

    try
        if exist('spm_dcm_peb','file') ~= 2
            error('spm_dcm_peb 不在路径中。');
        end
        PEB_cond = spm_dcm_peb(GCM, M, 'B');
        fprintf('  [PEB] spm_dcm_peb 完成。\n');

        if exist('spm_dcm_peb_bmc','file') ~= 2
            error('spm_dcm_peb_bmc 不在路径中。');
        end

        % post-hoc search + BMA（Occam window）
        BMA_cond = spm_dcm_peb_bmc(PEB_cond);
        fprintf('  [PEB] spm_dcm_peb_bmc 完成（post-hoc search + BMA）。\n');

    catch ME
        fprintf('  [PEB] 出错: %s\n', ME.message);
        PEB_results{ic} = [];
        continue;
    end

    % ---------- 计算 B 参数的后验概率 PP ----------
    Pp     = [];
    Pnames = {};
    Ep_vec = [];

    if ~isempty(BMA_cond) && isfield(BMA_cond, 'Ep') && isfield(BMA_cond, 'Cp')
        if exist('spm_vec','file') == 2 && exist('spm_Ncdf','file') == 2

            % 注意：这里全部转成 full，避免 sparse 导致 fprintf 报错
            Ep_vec = full(spm_vec(BMA_cond.Ep));   % 所有参数的后验均值向量
            Cp     = full(BMA_cond.Cp);            % 协方差矩阵

            sd     = sqrt(full(diag(Cp)));         % 标准差（对角线）
            Ep_vec = Ep_vec(:);
            sd     = sd(:);

            % 计算 PP：P(param ≠ 0 | Y)
            Pp = 1 - spm_Ncdf(0, abs(Ep_vec), sd);
            Pp = full(Pp(:));

            % 参数名字
            if isfield(BMA_cond, 'Pnames')
                Pnames = BMA_cond.Pnames(:);
            elseif isfield(PEB_cond, 'Pnames')
                Pnames = PEB_cond.Pnames(:);
            else
                Pnames = {};
            end
        else
            fprintf('  [PEB] 缺少 spm_vec 或 spm_Ncdf，无法计算 PP。\n');
        end
    else
        fprintf('  [PEB] BMA 缺少 Ep/Cp，无法计算 PP。\n');
    end

    % ---------- 打印 PP > 0.95 的 B 参数 ----------
    if ~isempty(Pp) && ~isempty(Pnames)
        fprintf('  [PEB/BMR] B 参数的后验概率 PP（>0.95 认为该连线在该条件下必要）：\n');
        thrPP  = 0.85;
        n_sigB = 0;
        for k = 1:numel(Pnames)
            name_k = Pnames{k};
            if isempty(name_k)
                continue;
            end
            % 只看以 'B' 开头的参数（即连接参数）
            if name_k(1) ~= 'B'
                continue;
            end
            if Pp(k) > thrPP
                fprintf('     %s: Ep=%.4f, PP=%.3f\n', name_k, Ep_vec(k), Pp(k));
                n_sigB = n_sigB + 1;
            end
        end
        if n_sigB == 0
            fprintf('     （无任何 B 参数的 PP 超过 %.2f）\n', thrPP);
        end
    else
        fprintf('  [PEB/BMR] 无法获得 B 参数的 PP 列表。\n');
    end

    % ---------- 保存结果 ----------
    PEB_results{ic} = struct();
    PEB_results{ic}.cond_name = cond;
    PEB_results{ic}.PEB       = PEB_cond;
    PEB_results{ic}.BMA       = BMA_cond;
    PEB_results{ic}.BMR       = BMR_cond;
    PEB_results{ic}.PP        = Pp;
    PEB_results{ic}.Pnames    = Pnames;
    PEB_results{ic}.Ep_vec    = Ep_vec;
    PEB_results{ic}.subjects  = subj_list;
end

%% ====================== 保存结果（供后续单独可视化脚本使用） ==================
out_file = fullfile(out_dir, sprintf('DCM_BMS_simple_%s_spm25.mat', session_label));
fprintf('\n[Save] 正在保存结果到: %s\n', out_file);
save(out_file, 'BMS_simple', 'subjects', 'conds', 'model_names', ...
     'B_mean_all', 'best_model_idx_all', 'B_diff_simple', ...
     'connwise_results', 'PEB_results', 'winning_model_info', ...
     'root_dir', 'analysis_dirname', 'session_label', 'out_dir', '-v7.3');
fprintf('[Save] 保存完成。\n');

%% ====================== 辅助函数：单样本检验 & FDR ==================
function p = local_one_sample_p(x)
% 对向量 x 做单样本检验 vs 0：
% 优先用 Lillietest 判断正态 -> ttest 或 signrank
x = x(:);
x = x(~isnan(x));
if numel(x) < 2
    p = NaN;
    return;
end
% 全零或几乎全相等时，差异为 0
if all(abs(x - x(1)) < 1e-12)
    p = 1;
    return;
end
if exist('lillietest','file') == 2
    try
        [h, ~] = lillietest(x);
        is_normal = (h == 0);
    catch
        is_normal = true;
    end
else
    is_normal = true;
end
if is_normal
    try
        [~, p] = ttest(x, 0);
    catch
        p = NaN;
    end
else
    try
        p = signrank(x, 0);
    catch
        p = NaN;
    end
end
end

function [p_thresh, sig_vec] = local_fdr_bh(p_vec, alpha)
% Benjamini-Hochberg FDR 控制
% 输入:
%   p_vec: 向量形式的 p 值
%   alpha: FDR 水平
% 输出:
%   p_thresh: 阈值；若无显著，则为 NaN
%   sig_vec : 逻辑向量，哪一些 p 通过 FDR（与 p_vec 同长度）

p_vec = p_vec(:);
m     = numel(p_vec);
[ps, idx] = sort(p_vec);
thr_idx = find(ps <= ( (1:m)' * alpha / m ), 1, 'last');
sig_vec = false(m,1);
if isempty(thr_idx)
    p_thresh = NaN;
else
    p_thresh = ps(thr_idx);
    sig_vec(idx(1:thr_idx)) = true;
end
end



