%% DCM_group_PEB_condContrast_spm25.m
%  组水平 PEB（跨条件）—— B 参数，带条件回归量：
%    M.X = [Mean, Cond(connect_vs_unconnect)]
%
%  核心逻辑：
%    - 对每个被试、每个条件（connect/unconnect），从所有 run 中选 F 最大的 DCM
%    - 选定一个统一的模型（unified_model_idx），确保所有被试在两条件下都是同一模型结构
%    - 构造 GCM：每个被试 2 行（connect, unconnect）
%    - 设计矩阵 M.X：
%          行 1: Mean=1, Cond=+1 (connect)
%          行 2: Mean=1, Cond=-1 (unconnect)
%          如此循环每个被试
%    - spm_dcm_peb(GCM, M, 'B') + spm_dcm_peb_bmc(PEB)
%    - 从 BMA 中分离出：
%          Ep_mean_B(i,j), Pp_mean_B(i,j)  -- 平均信息流
%          Ep_cond_B(i,j), Pp_cond_B(i,j)  -- 条件差异（connect vs unconnect）
%    - 重建每个条件的 group-level B：
%          B_connect_hat   = Ep_mean_B + Ep_cond_B
%          B_unconnect_hat = Ep_mean_B - Ep_cond_B
%
%  注意：
%    - 本脚本只进行估计和文本输出，不画网络图。
%    - 画图可以在单独的可视化脚本中基于本脚本输出的 Ep/Pp 矩阵来完成。

clear; clc; close all;

%% ====================== 路径 & 参数 ==========================
root_dir         = '/home/xue/data/prf_result';       % 数据根目录
analysis_dirname = 'DCM_manual_spm25_multiModel';     % 单被试 DCM 子目录
session_label    = 'ses-01';

% 被试列表（需要确保与单被试分析一致）
subjects = {'sub-03','sub-04','sub-05','sub-06','sub-07','sub-09','sub-10'};
nSub     = numel(subjects);

% 条件列表（必须是这两个名字，对应文件名中的 cond 字段）
conds    = {'connect','unconnect'};
if numel(conds) ~= 2
    error('本脚本假定只有两个条件（connect/unconnect），当前 conds 数量 = %d。', numel(conds));
end

% 模型名字（顺序必须和文件名中的 modelXX 对齐）
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

% ========= 统一使用哪一个模型做 PEB？（请按你的 BMS 结果手动设置） =========
% 例如：根据之前的 BMS / family-BMS 结果，选择你认为最合理的网络结构
% 比如这里统一使用 model 02 (B2_NPOhub)
unified_model_idx = 2;
if unified_model_idx < 1 || unified_model_idx > N_MODEL
    error('unified_model_idx=%d 无效（1-%d）。', unified_model_idx, N_MODEL);
end
unified_model_name = model_names{unified_model_idx};

fprintf('>>> 组水平 PEB（条件对比）将使用统一模型 %d (%s)\n', ...
    unified_model_idx, unified_model_name);
fprintf('>>> root_dir          = %s\n', root_dir);
fprintf('>>> session_label     = %s\n', session_label);
fprintf('>>> analysis_dirname  = %s\n', analysis_dirname);

% 输出目录
out_dir = fullfile(root_dir, 'DCM_group_PEB_condContrast');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% 检查 SPM 函数
if exist('spm_dcm_peb','file') ~= 2 || exist('spm_dcm_peb_bmc','file') ~= 2
    error('未检测到 spm_dcm_peb 或 spm_dcm_peb_bmc，请先 addpath SPM 并设置路径。');
end
if exist('spm_vec','file') ~= 2 || exist('spm_Ncdf','file') ~= 2
    error('缺少 spm_vec 或 spm_Ncdf，请确认 SPM 在 MATLAB path 中。');
end

%% ====================== 读取 DCM：每被试 × 每条件 × 指定模型 ======================
% 为 connect/unconnect 分别记录 best-F DCM
DCM_conn   = cell(nSub,1);
DCM_unconn = cell(nSub,1);
bestF_conn   = -Inf(nSub,1);
bestF_unconn = -Inf(nSub,1);

for is = 1:nSub
    subj    = subjects{is};
    dcm_dir = fullfile(root_dir, subj, session_label, analysis_dirname);
    fprintf('\n[Scan] subj=%s, DCM 目录: %s\n', subj, dcm_dir);

    if ~exist(dcm_dir, 'dir')
        fprintf('  [Warn] 目录不存在，跳过该被试。\n');
        continue;
    end

    files = dir(fullfile(dcm_dir, '*.mat'));
    fprintf('  [Info] 找到 %d 个 .mat 文件。\n', numel(files));
    if isempty(files)
        continue;
    end

    for k = 1:numel(files)
        fname    = files(k).name;
        fullpath = fullfile(dcm_dir, fname);

        % 文件名格式示例：
        %   sub-03_ses-01_DCM_connect_run01_model02_spm25.mat
        parts = strsplit(fname, '_');
        if numel(parts) < 6
            fprintf('  [Skip] 文件名字段数 (%d) < 6: %s\n', numel(parts), fname);
            continue;
        end

        subj_from_name = parts{1};  %#ok<NASGU>
        ses_from_name  = parts{2};  %#ok<NASGU>
        tag_from_name  = parts{3};  %#ok<NASGU>
        cond_from_name = parts{4};  % connect / unconnect
        run_str        = parts{5};  %#ok<NASGU>
        model_str      = parts{6};  % modelXX

        % 条件必须是 connect/unconnect 之一
        if ~strcmpi(cond_from_name, 'connect') && ~strcmpi(cond_from_name, 'unconnect')
            continue;
        end

        % 模型编号
        model_id = NaN;
        if strncmp(model_str, 'model', 5)
            model_id = str2double(model_str(6:end));
        end
        if isnan(model_id) || model_id ~= unified_model_idx
            % 只保留我们想用的统一模型
            continue;
        end

        % 加载 DCM
        S = load(fullpath, 'DCM');
        if ~isfield(S, 'DCM')
            fprintf('  [Skip] %s 中无 DCM 变量。\n', fname);
            continue;
        end
        DCM = S.DCM;

        if ~isfield(DCM, 'F') || ~isfinite(DCM.F)
            fprintf('  [Skip] %s 中缺少有效 F（log-evidence）。\n', fname);
            continue;
        end
        F = DCM.F;

        % 依据 cond_from_name 更新 best-F DCM
        if strcmpi(cond_from_name, 'connect')
            if F > bestF_conn(is)
                bestF_conn(is) = F;
                DCM_conn{is}   = DCM;
                fprintf('  [OK] subj=%s, connect: 发现更优 F=%.3f 的 DCM（%s）。\n', subj, F, fname);
            end
        elseif strcmpi(cond_from_name, 'unconnect')
            if F > bestF_unconn(is)
                bestF_unconn(is) = F;
                DCM_unconn{is}   = DCM;
                fprintf('  [OK] subj=%s, unconnect: 发现更优 F=%.3f 的 DCM（%s）。\n', subj, F, fname);
            end
        end
    end
end

%% ====================== 检查哪些被试两条件都有 DCM ======================
has_conn   = ~cellfun(@isempty, DCM_conn);
has_unconn = ~cellfun(@isempty, DCM_unconn);
use_mask   = has_conn & has_unconn;

if ~any(use_mask)
    error('没有任何被试同时具备 connect 和 unconnect 条件的 DCM（模型 %d）。', unified_model_idx);
end

subjects_used = subjects(use_mask);
nSub_used     = numel(subjects_used);
fprintf('\n[Summary] 统一模型 %d(%s) 下，有效被试数（两条件都有 DCM） = %d\n', ...
    unified_model_idx, unified_model_name, nSub_used);
for ii = 1:nSub_used
    fprintf('   #%d: %s\n', ii, subjects_used{ii});
end

%% ====================== 构造 GCM & 设计矩阵 M.X ======================
% GCM: 每被试两行：connect (+1) 和 unconnect (–1)
GCM      = {};
design_X = [];
id_subj  = {};
id_cond  = {};

for is = 1:nSub
    if ~use_mask(is)
        continue;
    end

    % connect
    GCM{end+1,1}      = DCM_conn{is};           %#ok<AGROW>
    design_X(end+1,:) = [1, +1];                %#ok<AGROW>
    id_subj{end+1,1}  = subjects{is};          %#ok<AGROW>
    id_cond{end+1,1}  = 'connect';             %#ok<AGROW>

    % unconnect
    GCM{end+1,1}      = DCM_unconn{is};        %#ok<AGROW>
    design_X(end+1,:) = [1, -1];               %#ok<AGROW>
    id_subj{end+1,1}  = subjects{is};          %#ok<AGROW>
    id_cond{end+1,1}  = 'unconnect';           %#ok<AGROW>
end

nRows = size(design_X,1);
fprintf('\n[Design] GCM 条目数（被试×条件） = %d\n', nRows);

M        = struct();
M.X      = design_X;
M.Xnames = {'Mean', 'Cond(connect_vs_unconnect)'};  % 对应列 1&2
M.Q      = 'all';

%% ====================== 运行 PEB + BMR/BMA ============================
fprintf('\n[PEB] 开始在 B 参数上拟合 PEB（统一模型 %d:%s, N_sub=%d）...\n', ...
    unified_model_idx, unified_model_name, nSub_used);

PEB = spm_dcm_peb(GCM, M, 'B');
fprintf('  [PEB] spm_dcm_peb 完成。\n');

BMA = spm_dcm_peb_bmc(PEB);
fprintf('  [PEB] spm_dcm_peb_bmc 完成（post-hoc search + BMA）。\n');

%% ====================== 从 BMA 中计算 Ep、PP 向量 ====================
Pnames = {};
if isfield(BMA, 'Pnames') && ~isempty(BMA.Pnames)
    Pnames = BMA.Pnames(:);
elseif isfield(PEB, 'Pnames')
    Pnames = PEB.Pnames(:);
else
    error('PEB/BMA 结果中找不到 Pnames 字段。');
end

% 所有 group-level 参数的后验均值（可能包含比自由参数更多的条目）
Ep_vec_full = full(spm_vec(BMA.Ep));
Ep_vec_full = Ep_vec_full(:);

% 协方差矩阵只对应真正的自由参数
Cp     = full(BMA.Cp);
sd_vec = sqrt(full(diag(Cp)));
sd_vec = sd_vec(:);

n_cov = numel(sd_vec);  % 自由参数个数（以 Cp 为准）

% 以 Cp 的长度为准截取 Ep_vec 和 Pnames
if numel(Ep_vec_full) < n_cov
    warning('Ep_vec_full 长度 (%d) 小于 Cp 对角线长度 (%d)，仅使用 Ep_vec_full 全部条目。', ...
        numel(Ep_vec_full), n_cov);
    n_use = numel(Ep_vec_full);
    Ep_vec = Ep_vec_full;
    sd_vec = sd_vec(1:n_use);
else
    n_use = n_cov;
    Ep_vec = Ep_vec_full(1:n_use);
end

if numel(Pnames) < n_use
    warning('Pnames 长度 (%d) 小于需要的参数个数 (%d)，将 Ep_vec/sd_vec 也截到 Pnames 长度。', ...
        numel(Pnames), n_use);
    n_use  = numel(Pnames);
    Ep_vec = Ep_vec(1:n_use);
    sd_vec = sd_vec(1:n_use);
elseif numel(Pnames) > n_use
    % Pnames 比 Cp/Ep_vec 长，直接截断多出来的
    Pnames = Pnames(1:n_use);
end

% 现在三者长度一致
if ~(numel(Ep_vec) == numel(sd_vec) && numel(Ep_vec) == numel(Pnames))
    error('长度对齐失败：Ep_vec=%d, sd_vec=%d, Pnames=%d。', ...
        numel(Ep_vec), numel(sd_vec), numel(Pnames));
end

% 计算每个参数的后验概率 P(|β|>0)
Pp_vec = 1 - spm_Ncdf(0, abs(Ep_vec), sd_vec);
Pp_vec = Pp_vec(:);

%% ====================== 将 B 参数映射回 ROI×ROI×(Mean/Cond) ==========
% 假定 B 参数名形如：
%   'B(i,j,k)' 其中 i=target, j=source, k=设计矩阵列（1=Mean, 2=Cond）
% 若缺少第三个索引，则默认 k=1（Mean）

% 先确定 ROI 数量
max_idx = 0;
for kk = 1:numel(Pnames)
    name_k = Pnames{kk};
    if isempty(name_k) || name_k(1) ~= 'B'
        continue;
    end
    tok = regexp(name_k, '^B\((\d+),(\d+)(?:,(\d+))?\)$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    ii = str2double(tok{1});
    jj = str2double(tok{2});
    if ~isnan(ii), max_idx = max(max_idx, ii); end
    if ~isnan(jj), max_idx = max(max_idx, jj); end
end

if max_idx == 0
    error('Pnames 中未检测到任何 B(i,j,*) 类型的参数。');
end

nROI = max_idx;
fprintf('\n[Map] 推断 ROI 数量 nROI = %d\n', nROI);

% 初始化矩阵：Mean 和 Cond 两个回归量上的 B(i,j)
Ep_mean_B = zeros(nROI, nROI);  % β_mean(i,j)
Pp_mean_B = zeros(nROI, nROI);
Ep_cond_B = zeros(nROI, nROI);  % β_cond(i,j) 对应 connect vs unconnect
Pp_cond_B = zeros(nROI, nROI);

nReg = size(M.X,2);
if nReg < 2
    warning('设计矩阵列数 < 2，Cond 回归量可能不存在。');
end

for kk = 1:numel(Pnames)
    name_k = Pnames{kk};
    if isempty(name_k) || name_k(1) ~= 'B'
        continue;
    end

    tok = regexp(name_k, '^B\((\d+),(\d+)(?:,(\d+))?\)$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end

    ii    = str2double(tok{1});  % target
    jj    = str2double(tok{2});  % source
    if numel(tok) >= 3 && ~isempty(tok{3})
        kReg = str2double(tok{3});
        if isnan(kReg), kReg = 1; end
    else
        % 若无第三个索引，默认属于第 1 列 Mean
        kReg = 1;
    end

    if isnan(ii) || isnan(jj) || ii < 1 || ii > nROI || jj < 1 || jj > nROI
        continue;
    end

    beta = Ep_vec(kk);
    pp   = Pp_vec(kk);

    % 若同一 (i,j,reg) 有多参数（例如 BMR 后的组合），简单保留 |Ep| 更大的那个
    if kReg == 1
        if abs(beta) > abs(Ep_mean_B(ii,jj))
            Ep_mean_B(ii,jj) = beta;
            Pp_mean_B(ii,jj) = pp;
        end
    elseif kReg == 2
        if abs(beta) > abs(Ep_cond_B(ii,jj))
            Ep_cond_B(ii,jj) = beta;
            Pp_cond_B(ii,jj) = pp;
        end
    else
        % 如果有更多回归量（>2），此处暂时忽略
        continue;
    end
end

%% ====================== 估计每个条件的 group-level B ==================
% 设计矩阵行模式:
%   connect:  Mean=1, Cond=+1
%   unconn :  Mean=1, Cond=–1
%
% 因此:
%   B_connect_hat   = β_mean + β_cond
%   B_unconnect_hat = β_mean - β_cond

B_connect_hat   = Ep_mean_B + Ep_cond_B;
B_unconnect_hat = Ep_mean_B - Ep_cond_B;

%% ====================== 从 DCM 中读取 ROI label（若失败用 ROI1..n） ====
roi_labels = {};
for is = 1:nSub
    if ~use_mask(is)
        continue;
    end
    D = DCM_conn{is};
    if isempty(D)
        D = DCM_unconn{is};
    end
    if isempty(D) || ~isfield(D, 'xY') || isempty(D.xY)
        continue;
    end
    try
        nx = numel(D.xY);
        tmp = cell(nx,1);
        for ii = 1:nx
            nm = D.xY(ii).name;   % 例如 '1 LNPC1'
            if ischar(nm)
                tok = regexp(nm, '^\s*\d+\s+(.*)$', 'tokens', 'once');
                if ~isempty(tok)
                    tmp{ii} = strtrim(tok{1});
                else
                    tmp{ii} = strtrim(nm);
                end
            else
                tmp{ii} = sprintf('ROI%d', ii);
            end
        end
        roi_labels = tmp;
        break;
    catch
        roi_labels = {};
    end
end

if isempty(roi_labels)
    fprintf('[ROI] 无法从 DCM 中读取 ROI 名称，使用默认 ROI1..ROI%d。\n', nROI);
    roi_labels = arrayfun(@(ii) sprintf('ROI%d', ii), 1:nROI, 'UniformOutput', false);
else
    % 保证长度与 nROI 一致
    if numel(roi_labels) < nROI
        roi_labels(numel(roi_labels)+1:nROI) = ...
            arrayfun(@(ii) sprintf('ROI%d', ii), numel(roi_labels)+1:nROI, 'UniformOutput', false);
    elseif numel(roi_labels) > nROI
        roi_labels = roi_labels(1:nROI);
    end
end

%% ====================== 文本输出：方向 & 条件差异 ======================
thrPP_main = 0.85;   % 平均效应阈值（可根据需要调整，例如 0.90）
thrPP_cond = 0.85;   % 条件差异阈值

fprintf('\n================= 结果摘要 =================\n');
fprintf('统一模型: %d (%s), ROI 数: %d, 有效被试: %d\n', ...
    unified_model_idx, unified_model_name, nROI, nSub_used);

% 1) 平均信息流（Mean 列）
fprintf('\n[Mean] 平均条件下（connect+unconnect）PP > %.2f 的 B 连接（方向: source -> target）：\n', thrPP_main);
any_mean = false;
for ii = 1:nROI
    for jj = 1:nROI
        if ii == jj, continue; end
        if Pp_mean_B(ii,jj) > thrPP_main
            any_mean = true;
            fprintf('  %s <- %s : Ep_mean = %.4f, PP = %.3f\n', ...
                roi_labels{ii}, roi_labels{jj}, Ep_mean_B(ii,jj), Pp_mean_B(ii,jj));
        end
    end
end
if ~any_mean
    fprintf('  (无任何 B(i,j) 在 Mean 上 PP 超过 %.2f)\n', thrPP_main);
end

% 2) 条件差异（Cond 列：connect_vs_unconnect）
fprintf('\n[Cond] 条件效应（connect vs unconnect）PP > %.2f 的 B 连接（方向: source -> target）：\n', thrPP_cond);
fprintf('       Ep_cond>0 -> connect > unconnect ; Ep_cond<0 -> unconnect > connect\n');
any_cond = false;
for ii = 1:nROI
    for jj = 1:nROI
        if ii == jj, continue; end
        if Pp_cond_B(ii,jj) > thrPP_cond
            any_cond = true;
            if Ep_cond_B(ii,jj) >= 0
                dir_str = 'connect > unconnect';
            else
                dir_str = 'unconnect > connect';
            end
            fprintf('  %s <- %s : Ep_cond = %.4f, PP = %.3f  (%s)\n', ...
                roi_labels{ii}, roi_labels{jj}, Ep_cond_B(ii,jj), Pp_cond_B(ii,jj), dir_str);
        end
    end
end
if ~any_cond
    fprintf('  (无任何 B(i,j) 在 Cond 上 PP 超过 %.2f)\n', thrPP_cond);
end

%% ====================== 保存结果 ===========================
out_file = fullfile(out_dir, sprintf('DCM_group_PEB_condContrast_%s_model%02d_spm25.mat', ...
    session_label, unified_model_idx));

fprintf('\n[Save] 正在保存结果到: %s\n', out_file);

save(out_file, ...
    'GCM', 'M', ...
    'root_dir', 'analysis_dirname', 'session_label', ...
    'subjects', 'subjects_used', 'conds', ...
    'model_names', 'unified_model_idx', 'unified_model_name', ...
    'PEB', 'BMA', ...
    'Pnames', 'Ep_vec', 'Pp_vec', ...
    'Ep_mean_B', 'Pp_mean_B', ...
    'Ep_cond_B', 'Pp_cond_B', ...
    'B_connect_hat', 'B_unconnect_hat', ...
    'roi_labels', ...
    'id_subj', 'id_cond', ...
    '-v7.3');

fprintf('[Save] 保存完成。\n');
