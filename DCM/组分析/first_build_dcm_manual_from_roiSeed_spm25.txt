%% ========================================================================
%  DCM run-wise QC + 简单 BMS 脚本（多模型版本）
%
%  假定你的 DCM 文件来自 build_dcm_manual_multiModel_spm25.m：
%    root_prf/sub-XX/ses-01/DCM_manual_spm25_multiModel/
%      sub-06_ses-01_DCM_connect_run01_model01_spm25.mat
%      sub-06_ses-01_DCM_connect_run01_model02_spm25.mat
%      ...
%      sub-06_ses-01_DCM_unconnect_run01_model07_spm25.mat
%
%  步骤：
%    1) 扫描每个被试的 DCM 目录，找到所有 DCM 文件
%    2) 从文件名中解析 cond、run、model
%    3) 读取每个 DCM 的 F 和 mean R²
%    4) 对每个 subj × run 做“单被试 BMS”
%         - 将该 subj + run 的所有有效 DCM 当作候选模型 family
%         - 用 F 作为 log evidence，计算 posterior 概率，选 winning model
%    5) 用 winning model 的 F & R² 做 run-wise QC（R² 阈值 + 极端 outlier）
%    6) 打印建议剔除的 run，以及每个 model 的估计失败率
%    7) 将 QC 结果保存到 /home/xue/data/prf_result/DCM_QC/ 下
%
%  备注：
%    - 对单个数据集（一个 run）来说，BMS 在数学上就是选择 F 最大的模型，
%      这里额外算了一下 posterior 方便你后续看不确定性。
% ========================================================================

clear; clc;

%% ---------------------- 基本参数设置（真实路径） ------------------------

% 1) 和 DCM 构建脚本保持一致
root_prf   = '/home/xue/data/prf_result';      % build_dcm_manual_multiModel_spm25 里的 root_prf
prf_ses    = 'ses-01';                         % session 名
dcm_subdir = 'DCM_manual_spm25_multiModel';    % 构建脚本 out_dir 的最后一级

% 2) 被试列表（你自己控制）
subject_list = { 'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09','sub-10' };

% 3) QC 阈值（可以按需要调）
R2_min   = 0.05;   % winning DCM 的 mean R² 低于这个认为 run 质量较差
z_thresh = 3;      % F / R² 的 z 分数 outlier 阈值（>3 视为极端）

% 4) 模型名字（可选，用于输出更友好）
modelNames = { ...
    'B1_globalGain', ...
    'B2_NPOhub', ...
    'B3_feedforward', ...
    'B4_feedback', ...
    'B5_parietalOnly', ...
    'B6_hierarchy', ...
    'B0_null'};

%% ---------------------- 存放结果的结构体 --------------------------------

QC = struct('subj', {}, ...         % subj ID
            'run',  {}, ...         % run 编号
            'model', {}, ...        % model 编号（从文件名或 DCM 中解析）
            'model_name', {}, ...   % 模型名字（如果能匹配）
            'cond', {}, ...         % 条件名（connect / unconnect 等）
            'F',    {}, ...         % DCM.F
            'R2',   {}, ...         % mean R²
            'is_valid', {}, ...     % 该 DCM 是否估计有效
            'note', {}, ...         % 错误 / 说明
            'file', {});            % 完整文件名

RunSummary = struct('subj', {}, 'run', {}, ...
                    'win_model', [], 'win_model_name', {}, ...
                    'win_cond', {}, ...
                    'F', [], 'R2', [], ...
                    'post_prob', [], ...  % winning model 的 posterior 概率
                    'qc_pass', [], 'reasons', {});

idx = 0;

%% ---------------------- 遍历被试 & 自动识别 run / model -----------------

for s = 1:numel(subject_list)
    subj = subject_list{s};
    
    % 对应 DCM 目录：root_prf/subj/ses-01/DCM_manual_spm25_multiModel
    subj_dcm_dir = fullfile(root_prf, subj, prf_ses, dcm_subdir);
    if ~isfolder(subj_dcm_dir)
        warning('找不到 DCM 目录: %s (跳过 %s)', subj_dcm_dir, subj);
        continue;
    end
    
    % 列出该目录下所有符合 pattern 的 DCM 文件
    % 例如:
    %   sub-06_ses-01_DCM_connect_run01_spm25.mat
    %   sub-06_ses-01_DCM_connect_run01_model01_spm25.mat
    patt = sprintf('%s_%s_DCM_*_run*_spm25.mat', subj, prf_ses);
    dcm_files = dir(fullfile(subj_dcm_dir, patt));
    
    if isempty(dcm_files)
        warning('在 %s 下找不到任何 DCM 文件 (跳过 %s)', subj_dcm_dir, subj);
        continue;
    end
    
    fprintf('\n[SUBJ] %s: 在 %s 下找到 %d 个 DCM 文件\n', ...
        subj, subj_dcm_dir, numel(dcm_files));
    
    for f = 1:numel(dcm_files)
        fn      = dcm_files(f).name;
        full_fn = fullfile(dcm_files(f).folder, fn);
        
        % 从文件名中解析 cond、run、model
        % 允许形式：
        %   sub-06_ses-01_DCM_connect_run01_spm25.mat
        %   sub-06_ses-01_DCM_connect_run01_model01_spm25.mat
        %   或旧版：sub-06_ses-01_DCM_connect_run01_model-1_spm25.mat
        expr = sprintf( ...
            '%s_%s_DCM_(?<cond>.+?)_run(?<run>\\d+)(?:_model-?(?<model>\\d+))?_spm25\\.mat', ...
            subj, prf_ses);
        tokens = regexp(fn, expr, 'names');
        
        if isempty(tokens)
            warning('  [Skip] 文件名不符合预期模式，无法解析 run/model: %s', fn);
            continue;
        end
        
        condLabel = tokens.cond;
        run_idx   = str2double(tokens.run);
        if isfield(tokens,'model') && ~isempty(tokens.model)
            model_idx = str2double(tokens.model);
        else
            model_idx = 1;  % 没有写 model 的就当作 1
        end
        
        % ==== 初始化当前 DCM 的 QC 记录 ====
        idx = idx + 1;
        QC(idx).subj       = subj;
        QC(idx).run        = run_idx;
        QC(idx).model      = model_idx;
        QC(idx).model_name = '';   % 后面尝试补上
        QC(idx).cond       = condLabel;
        QC(idx).F          = NaN;
        QC(idx).R2         = NaN;
        QC(idx).is_valid   = false;
        QC(idx).note       = '';
        QC(idx).file       = full_fn;
        
        % 尝试从 modelNames 映射一个名字（如果合理）
        if model_idx >= 1 && model_idx <= numel(modelNames)
            QC(idx).model_name = modelNames{model_idx};
        end
        
        % ==== 尝试加载 DCM，并读取 F / R² ====
        try
            S = load(full_fn, 'DCM');
            if ~isfield(S, 'DCM')
                QC(idx).note = '没有 DCM 变量';
                continue;
            end
            DCM = S.DCM;
            
            % 如果 DCM 里自带 model_id / model_name，可以用它覆盖
            if isfield(DCM, 'model_id') && isnumeric(DCM.model_id)
                QC(idx).model = DCM.model_id;
                model_idx     = DCM.model_id;
            end
            if isfield(DCM, 'model_name') && ischar(DCM.model_name)
                QC(idx).model_name = DCM.model_name;
            end
            
            % ------- 1) F -------
            if isfield(DCM, 'F') && isfinite(DCM.F)
                QC(idx).F = DCM.F;
            else
                QC(idx).note = 'F 缺失或非有限';
            end
            
            % ------- 2) mean R² -------
            meanR2 = NaN;
            
            % 优先使用 DCM.R2
            if isfield(DCM, 'R2') && ~isempty(DCM.R2)
                tmp = DCM.R2(:);
                tmp = tmp(~isnan(tmp));
                if ~isempty(tmp)
                    meanR2 = mean(tmp);
                end
            else
                % 没有 DCM.R2，就尝试用 Y.y vs DCM.y 计算
                if isfield(DCM, 'Y') && isfield(DCM.Y, 'y') && isfield(DCM, 'y')
                    try
                        Y_data = DCM.Y.y; % 真实 BOLD
                        Y_hat  = DCM.y;   % 预测 BOLD
                        
                        if iscell(Y_data)
                            Y_data = cell2mat(Y_data);
                        end
                        if iscell(Y_hat)
                            Y_hat = cell2mat(Y_hat);
                        end
                        
                        if ~isempty(Y_data) && ~isempty(Y_hat) && all(size(Y_data) == size(Y_hat))
                            nReg = size(Y_data, 2);
                            R2_per_reg = nan(1, nReg);
                            for rr = 1:nReg
                                y  = Y_data(:, rr);
                                yh = Y_hat(:, rr);
                                
                                y  = y(:);
                                yh = yh(:);
                                
                                if any(isnan(y)) || any(isnan(yh))
                                    continue;
                                end
                                
                                y_mean = mean(y);
                                SSE = sum((y - yh).^2);
                                SST = sum((y - y_mean).^2);
                                
                                if SST > 0
                                    R2_per_reg(rr) = 1 - SSE / SST;
                                end
                            end
                            
                            tmp = R2_per_reg(~isnan(R2_per_reg));
                            if ~isempty(tmp)
                                meanR2 = mean(tmp);
                            end
                        else
                            QC(idx).note = 'Y.y 与 DCM.y 维度不匹配';
                        end
                    catch ME
                        QC(idx).note = sprintf('计算 R² 失败: %s', ME.message);
                    end
                else
                    QC(idx).note = '无 R2 字段，也无法根据 Y.y/DCM.y 计算';
                end
            end
            
            QC(idx).R2 = meanR2;
            
            % ------- 3) 判定该 DCM 是否“估计有效” -------
            if isfinite(QC(idx).F) && ~isnan(QC(idx).R2)
                QC(idx).is_valid = true;
            else
                if isempty(QC(idx).note)
                    QC(idx).note = 'F 或 R² 非有效';
                end
            end
            
        catch ME
            QC(idx).note     = sprintf('加载/解析 DCM 失败: %s', ME.message);
            QC(idx).is_valid = false;
        end
        
    end % dcm_files
end % subjects

if isempty(QC)
    error('没有任何 DCM 被成功读取，请检查路径和文件命名。');
end

fprintf('\n========== 第一步：逐 DCM 估计情况汇总完成 ==========\n');

%% ---------------------- run 级别：BMS + winning DCM QC -------------------

% 这里的 “run” 是 subj 内的 run 编号，不区分 cond。
% 即：对同一 subj + run 的所有 DCM（不同 cond & model）一起做 BMS。

keys      = arrayfun(@(q) sprintf('%s__run%02d', q.subj, q.run), QC, 'UniformOutput', false);
uniq_keys = unique(keys);

RunSummary = struct('subj', {}, 'run', {}, ...
                    'win_model', [], 'win_model_name', {}, ...
                    'win_cond', {}, ...
                    'F', [], 'R2', [], ...
                    'post_prob', [], ...
                    'qc_pass', [], 'reasons', {});

rs_idx = 0;

for u = 1:numel(uniq_keys)
    key   = uniq_keys{u};
    parts = split(key, '__run');
    subj  = parts{1};
    run   = str2double(parts{2});
    
    mask  = strcmp({QC.subj}, subj) & [QC.run] == run;
    these = QC(mask);
    
    % 只看估计有效的 DCM
    valid_mask = [these.is_valid];
    if ~any(valid_mask)
        % 这个 run 所有 DCM 都挂了 → 直接认为 run 有问题
        rs_idx = rs_idx + 1;
        RunSummary(rs_idx).subj           = subj;
        RunSummary(rs_idx).run            = run;
        RunSummary(rs_idx).win_model      = NaN;
        RunSummary(rs_idx).win_model_name = '';
        RunSummary(rs_idx).win_cond       = '';
        RunSummary(rs_idx).F              = NaN;
        RunSummary(rs_idx).R2             = NaN;
        RunSummary(rs_idx).post_prob      = NaN;
        RunSummary(rs_idx).qc_pass        = false;
        RunSummary(rs_idx).reasons        = { ...
            '该 run 的所有 DCM 都估计失败（F/R² 无效），可能设计或代码有严重问题'};
        continue;
    end
    
    valid_these = these(valid_mask);
    allF        = [valid_these.F];
    
    % ---- “单被试 BMS”：用 F 作为 log evidence 计算 posterior ----
    % 为数值稳定，用 F - max(F)
    F_max = max(allF);
    dF    = allF - F_max;
    post  = exp(dF);
    post  = post / sum(post);
    
    % 选 posterior 最大（等价于 F 最大）的 DCM 作为 winning DCM
    [~, k] = max(post);
    win    = valid_these(k);
    
    rs_idx = rs_idx + 1;
    RunSummary(rs_idx).subj           = subj;
    RunSummary(rs_idx).run            = run;
    RunSummary(rs_idx).win_model      = win.model;
    RunSummary(rs_idx).win_model_name = win.model_name;
    RunSummary(rs_idx).win_cond       = win.cond;
    RunSummary(rs_idx).F              = win.F;
    RunSummary(rs_idx).R2             = win.R2;
    RunSummary(rs_idx).post_prob      = post(k);
    RunSummary(rs_idx).qc_pass        = true;   % 先默认通过
    RunSummary(rs_idx).reasons        = {};
end

% 计算 winning DCM 的 F / R² 的 z-score，用于 outlier 检测
allF_win  = [RunSummary.F];
allR2_win = [RunSummary.R2];

valid_F  = isfinite(allF_win);
valid_R2 = ~isnan(allR2_win);

zF  = nan(size(allF_win));
zR2 = nan(size(allR2_win));

if any(valid_F)
    muF = mean(allF_win(valid_F));
    sdF = std(allF_win(valid_F));
    if sdF > 0
        zF(valid_F) = (allF_win(valid_F) - muF) / sdF;
    else
        zF(valid_F) = 0;  % 所有 F 一样就都不算 outlier
    end
end
if any(valid_R2)
    muR2 = mean(allR2_win(valid_R2));
    sdR2 = std(allR2_win(valid_R2));
    if sdR2 > 0
        zR2(valid_R2) = (allR2_win(valid_R2) - muR2) / sdR2;
    else
        zR2(valid_R2) = 0;
    end
end

% 根据阈值更新 qc_pass & reasons
for i = 1:numel(RunSummary)
    rs = RunSummary(i);
    reasons = {};
    pass    = true;
    
    % 1) F 是否有限
    if ~isfinite(rs.F)
        pass = false;
        reasons{end+1} = 'winning DCM 的 F 非有限';
    end
    
    % 2) R² 是否太低
    if isnan(rs.R2)
        pass = false;
        reasons{end+1} = 'winning DCM 的 R² = NaN';
    elseif rs.R2 < R2_min
        pass = false;
        reasons{end+1} = sprintf('winning DCM 的 mean R² (%.3f) < 阈值 %.3f', rs.R2, R2_min);
    end
    
    % 3) F / R² 是否为极端 outlier（可选）
    if ~isnan(zF(i)) && abs(zF(i)) > z_thresh
        pass = false;
        reasons{end+1} = sprintf('winning DCM 的 F 是极端 outlier (z=%.2f)', zF(i));
    end
    if ~isnan(zR2(i)) && abs(zR2(i)) > z_thresh
        pass = false;
        reasons{end+1} = sprintf('winning DCM 的 R² 是极端 outlier (z=%.2f)', zR2(i));
    end
    
    RunSummary(i).qc_pass = pass;
    RunSummary(i).reasons = reasons;
end

fprintf('\n========== 第二步：按 subj × run 的 BMS + QC 结果 ==========\n\n');

% 打印建议剔除的 run
bad_idx = find(~[RunSummary.qc_pass]);
if isempty(bad_idx)
    fprintf('>> 所有 run 的 winning DCM 在当前阈值下都通过 QC ✅\n');
else
    fprintf('>> 建议在后续 PEB / BMS 中剔除（或至少人工检查）的 run：\n');
    for k = bad_idx
        rs = RunSummary(k);
        if isempty(rs.win_model_name)
            mstr = sprintf('model=%d', rs.win_model);
        else
            mstr = sprintf('model=%d (%s)', rs.win_model, rs.win_model_name);
        end
        fprintf('  - %s, run %02d, winning DCM: cond=%s, %s, F=%.3f, R²=%.3f, post=%.3f\n', ...
            rs.subj, rs.run, rs.win_cond, mstr, rs.F, rs.R2, rs.post_prob);
        for jj = 1:numel(rs.reasons)
            fprintf('      * %s\n', rs.reasons{jj});
        end
    end
end

%% ---------------------- model 级别：看哪个 model 经常挂 -----------

all_model = [QC.model];
all_valid = [QC.is_valid];

uniq_models = unique(all_model(~isnan(all_model)));

fprintf('\n========== 第三步：按 model 的估计稳定性检查 ==========\n\n');

for mm = 1:numel(uniq_models)
    m = uniq_models(mm);
    mask_m = (all_model == m);
    total_m   = sum(mask_m);
    valid_m   = sum(mask_m & all_valid);
    invalid_m = total_m - valid_m;
    ratio_bad = invalid_m / max(total_m,1);
    
    if m >= 1 && m <= numel(modelNames)
        mname = modelNames{m};
    else
        mname = sprintf('model-%d', m);
    end
    
    fprintf('Model %d (%s): 总计 %d 个 DCM，估计有效 %d 个，失败 %d 个 (失败率 %.1f%%)\n', ...
        m, mname, total_m, valid_m, invalid_m, 100*ratio_bad);
    
    if ratio_bad > 0.5
        fprintf('  >> 警告：Model %d 的估计失败率 > 50%%，可能该模型的设计/脚本有问题，请重点检查。\n', m);
    end
end

%% ---------------------- 保存 QC 结果到 group 目录 -----------------------

qc_dir = fullfile(root_prf, 'DCM_QC');
if ~exist(qc_dir, 'dir')
    mkdir(qc_dir);
end

out_mat = fullfile(qc_dir, sprintf('DCM_QC_runwise_BMS_%s.mat', prf_ses));
save(out_mat, 'QC', 'RunSummary', 'root_prf', 'prf_ses', 'dcm_subdir', ...
              'subject_list', 'R2_min', 'z_thresh');

fprintf('\nQC 结果已保存到: %s\n', out_mat);
fprintf('脚本结束。\n');

