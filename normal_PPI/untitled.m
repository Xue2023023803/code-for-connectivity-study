%% ppi_graph_analysis_group.m
% 基于 bilateral PPI 结果 (results_fix.beta_PPI_bilat) 做简单图论分析
% - 节点强度（strength）
% - 度数（degree，基于绝对值归一化后的阈值）
% - 聚类系数（clustering coefficient）
% + 全局强度（Global strength / mean connectivity = mean(|W|)）
% 并在 group-level 上做平均，并画出每个条件的网络图

clear; clc; close all;

%% ---------- 基本设置（按需修改） ----------
root_bids = '/home/xue/data/BIDS_prep';  % <--- 改成你的 BIDS_prep 根目录

% 参与 PPI 分析的被试列表（按你的实际情况修改）
subjects = {  'sub-01','sub-02','sub-03', 'sub-04','sub-05', 'sub-06', 'sub-07','sub-08','sub-09', 'sub-10'};

condLevels = {'connect','unconnect'};  % 要分析的条件

% 图论阈值：对每个被试、每个条件，将 |W| 归一化到 [0,1] 再以此阈值生成二值图
thr_norm = 0.3;  % 可以之后改成 0.2 / 0.4 看鲁棒性

%% ---------- 预分配（等第一次 load 确定 ROI 数） ----------
nSub = numel(subjects);
baseROI_list = {};
nBase = [];

% 我们想要保存的 group-level 指标：
% strength(ROI × cond × subj), degree(ROI × cond × subj), clustering(ROI × cond × subj)
strength_all   = [];
degree_all     = [];
cluster_all    = [];

% 仅保留一个 global-level 指标：Global strength / mean connectivity（cond × subj）
global_meanConn_all = [];

for s = 1:nSub
    subj = subjects{s};
    fprintf('\n===== 处理被试 %s =====\n', subj);

    % FFX bilateral 结果文件路径
    ffx_file = fullfile(root_bids, 'derivatives', 'ppi_roi_runwise', subj, ...
        sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral.mat', subj));

    if ~exist(ffx_file, 'file')
        warning('找不到 FFX 文件: %s，跳过该被试。', ffx_file);
        continue;
    end

    S = load(ffx_file, 'results_fix');
    if ~isfield(S, 'results_fix')
        warning('文件 %s 中没有 results_fix，跳过该被试。', ffx_file);
        continue;
    end
    results_fix = S.results_fix;

    if isempty(baseROI_list)
        % 第一个有效被试，初始化 ROI 信息和存储矩阵
        baseROI_list = results_fix.baseROI_list;
        nBase = numel(baseROI_list);

        nCond = numel(condLevels);
        strength_all = nan(nBase, nCond, nSub);
        degree_all   = nan(nBase, nCond, nSub);
        cluster_all  = nan(nBase, nCond, nSub);

        % 初始化 global-level 存储
        global_meanConn_all = nan(nCond, nSub);

        fprintf('    [Info] 检测到 base ROI 数 = %d\n', nBase);
        disp(baseROI_list(:)');
    end

    for ci = 1:numel(condLevels)
        condLabel = condLevels{ci};

        if ~isfield(results_fix, condLabel) || ...
           ~isfield(results_fix.(condLabel), 'beta_PPI_bilat') || ...
           isempty(results_fix.(condLabel).beta_PPI_bilat)

            warning('    [Subj %s] 条件 %s 没有 beta_PPI_bilat，填 NaN。', ...
                subj, condLabel);
            continue;
        end

        W = results_fix.(condLabel).beta_PPI_bilat;   % [nBase × nBase]

        % 对称化：这里先简单取 (W + W')/2，保留符号
        W_sym = (W + W.') / 2;

        % 主对角线置 0（不考虑自连）
        W_sym(1:nBase+1:end) = 0;

        % --- 计算图论指标 ---
        % 返回：节点级 strength/degree/clustering + global mean|W|
        [node_strength, node_degree, node_cluster, ~, ~, global_meanConn] = ...
            compute_graph_measures(W_sym, thr_norm);

        strength_all(:, ci, s) = node_strength;
        degree_all(:,   ci, s) = node_degree;
        cluster_all(:,  ci, s) = node_cluster;

        % 保存 global-level 指标（cond × subj）
        global_meanConn_all(ci, s) = global_meanConn;

        fprintf('    [Subj %s, Cond %s] strength/degree/clustering 已计算。\n', ...
            subj, condLabel);
        fprintf('    [Subj %s, Cond %s] global mean|W| 已计算。\n', ...
            subj, condLabel);
    end
end

%% ---------- Group-level 统计：均值 ± SEM 并打印 ----------
fprintf('\n================ Group-level 图论指标 (mean ± SEM) ================\n');

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    % 取在该条件下不全是 NaN 的被试
    data_strength = squeeze(strength_all(:, ci, :));   % [ROI × subj]
    data_degree   = squeeze(degree_all(:,   ci, :));
    data_cluster  = squeeze(cluster_all(:,  ci, :));

    valid_subj = ~all(isnan(data_strength), 1);
    nValid = sum(valid_subj);

    if nValid == 0
        fprintf('\n[Cond %s] 无有效被试，跳过。\n', condLabel);
        continue;
    end

    data_strength = data_strength(:, valid_subj);
    data_degree   = data_degree(:,   valid_subj);
    data_cluster  = data_cluster(:,  valid_subj);

    m_strength = mean(data_strength, 2, 'omitnan');
    m_degree   = mean(data_degree,   2, 'omitnan');
    m_cluster  = mean(data_cluster,  2, 'omitnan');

    sem_strength = std(data_strength, 0, 2, 'omitnan') ./ sqrt(nValid);
    sem_degree   = std(data_degree,   0, 2, 'omitnan') ./ sqrt(nValid);
    sem_cluster  = std(data_cluster,  0, 2, 'omitnan') ./ sqrt(nValid);

    % global mean|W|：按同一 valid_subj 做 mean ± SEM
    g_meanConn = global_meanConn_all(ci, valid_subj);
    m_g_meanConn   = mean(g_meanConn, 'omitnan');
    sem_g_meanConn = std(g_meanConn, 0, 'omitnan') ./ sqrt(nValid);

    fprintf('\n===== 条件 %s：N = %d =====\n', condLabel, nValid);

    % 仅打印一个全局指标
    fprintf('GLOBAL | mean|W| = %6.3f ± %5.3f\n', m_g_meanConn, sem_g_meanConn);

    % 原有 ROI-level 指标输出保持不变
    for r = 1:nBase
        fprintf('ROI %-6s | strength = %6.3f ± %5.3f | degree = %4.2f ± %4.2f | C = %5.3f ± %5.3f\n', ...
            baseROI_list{r}, ...
            m_strength(r), sem_strength(r), ...
            m_degree(r),   sem_degree(r), ...
            m_cluster(r),  sem_cluster(r));
    end
end

%% ---------- 可视化：条件层面的平均网络图 ----------
fprintf('\n================ 绘制每个条件的 group-level 网络图 ================\n');

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    % 取该条件下所有被试的 W_sym（我们重新从文件读取一遍）
    W_group = zeros(nBase, nBase);
    count_W = zeros(nBase, nBase);

    for s = 1:nSub
        subj = subjects{s};
        ffx_file = fullfile(root_bids, 'derivatives', 'ppi_roi_runwise', subj, ...
            sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral.mat', subj));
        if ~exist(ffx_file, 'file')
            continue;
        end
        S = load(ffx_file, 'results_fix');
        if ~isfield(S, 'results_fix') || ~isfield(S.results_fix, condLabel)
            continue;
        end
        res = S.results_fix.(condLabel);
        if ~isfield(res, 'beta_PPI_bilat') || isempty(res.beta_PPI_bilat)
            continue;
        end

        W = res.beta_PPI_bilat;
        W_sym = (W + W.') / 2;
        W_sym(1:nBase+1:end) = 0;

        mask_nonNan = ~isnan(W_sym);
        W_group(mask_nonNan) = W_group(mask_nonNan) + W_sym(mask_nonNan);
        count_W(mask_nonNan) = count_W(mask_nonNan) + 1;
    end

    W_mean = W_group ./ max(count_W, 1);  % 避免除 0
    W_mean(count_W == 0) = NaN;

    % 用同样的 thr_norm 生成二值图（基于 |W_mean| 的归一化）
    [node_strength_g, node_degree_g, ~, A_bin, W_abs_norm] = ...
        compute_graph_measures(W_mean, thr_norm);

    % 圆形布局
    theta = linspace(0, 2*pi, nBase+1);
    theta = theta(1:nBase);
    R = 1.0;
    x = R * cos(theta);
    y = R * sin(theta);

    % 节点大小 ~ group-level strength
    max_strength = max(node_strength_g);
    if max_strength <= 0
        node_size = 300 * ones(nBase,1);
    else
        node_size = 200 + 400 * (node_strength_g / max_strength); % [200, 600]
    end

    % 画图
    figure('Name', sprintf('PPI Graph - %s (group mean)', condLabel));
    set(gcf, 'Color', [1 1 1]);
    hold on; axis equal off;

    % 画边（上三角）
    for i = 1:nBase-1
        for j = i+1:nBase
            if isnan(W_mean(i,j))
                continue;
            end

            % 只画通过二值阈值的边
            if ~A_bin(i,j)
                continue;
            end

            wij = W_mean(i,j);
            % 颜色：正 β 用红色，负 β 用蓝色，中性用灰色
            if wij > 0
                edge_color = [1 0 0];      % 红
            elseif wij < 0
                edge_color = [0 0.4 1];    % 蓝
            else
                edge_color = [0.7 0.7 0.7];
            end

            % 线宽 ~ |W_mean| (归一化后)
            wij_norm = W_abs_norm(i,j);
            lw = 1 + 4 * wij_norm;        % [1,5]

            line([x(i) x(j)], [y(i) y(j)], 'Color', edge_color, 'LineWidth', lw);
        end
    end

    % 画节点
    for i = 1:nBase
        scatter(x(i), y(i), node_size(i), ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5);
        text(x(i)*1.15, y(i)*1.15, baseROI_list{i}, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end

    title(sprintf('PPI Graph (%s) - group mean (thr=%.2f)', condLabel, thr_norm));
end

fprintf('\n>>> 图论分析和可视化完成，可以查看命令行输出和图像窗口。\n');

%% ---------- 子函数：计算图论指标 ----------
function [node_strength, node_degree, node_cluster, A_bin, W_abs_norm, global_meanConn] = ...
    compute_graph_measures(W_sym, thr_norm)
% 输入:
%   W_sym    : n×n 对称加权矩阵（可以有 NaN；主对角线会被置 0）
%   thr_norm : 阈值（基于 |W| 归一化后 [0,1]；例 thr_norm = 0.3）
%
% 输出:
%   node_strength  : 每个节点的加权强度 sum(|w_ij|)
%   node_degree    : 在二值图下的度数
%   node_cluster   : 在二值图下的聚类系数
%   A_bin          : 阈值后的二值邻接矩阵
%   W_abs_norm     : 归一化后的 |W|，供画图使用
%   global_meanConn: 全网 mean(|w_ij|)（上三角、去对角、忽略 NaN）

    n = size(W_sym,1);

    % 主对角线归零
    W_sym(1:n+1:end) = 0;

    % 绝对值加权
    W_abs = abs(W_sym);
    mask_valid = ~isnan(W_abs);

    max_w = max(W_abs(mask_valid));
    if isempty(max_w) || max_w == 0
        W_abs_norm = zeros(n);
    else
        W_abs_norm = W_abs / max_w;  % 归一化到 [0,1]
    end

    % 二值邻接矩阵
    A_bin = W_abs_norm > thr_norm;
    A_bin = triu(A_bin,1);
    A_bin = A_bin + A_bin.';  % 对称化
    A_bin(1:n+1:end) = 0;

    % 节点强度：sum |w_ij|
    node_strength = sum(W_abs, 2, 'omitnan');

    % 节点度数：在二值图中有多少邻居
    node_degree = sum(A_bin, 2);

    % 聚类系数（无向二值图）
    node_cluster = zeros(n,1);
    for i = 1:n
        neigh = find(A_bin(i,:));
        k = numel(neigh);
        if k < 2
            node_cluster(i) = 0;
            continue;
        end
        % 子图中的边数
        subA = A_bin(neigh, neigh);
        e = sum(subA(:)) / 2;  % 每条边算了两次
        node_cluster(i) = 2*e / (k*(k-1));
    end

    % 全局 mean|W|：上三角边（去对角），忽略 NaN
    ut = triu(true(n), 1);
    idx = ut & ~isnan(W_abs);
    if any(idx(:))
        global_meanConn = mean(W_abs(idx), 'omitnan');
    else
        global_meanConn = 0;
    end
end

