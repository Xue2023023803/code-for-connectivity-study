%% ppi_graph_analysis_group.m
% 基于 bilateral PPI 结果 (results_fix.(cond).beta_PPI_bilat) 做简单图论分析
% - 节点强度（strength, sum(|w_ij|)）
% - 度数（degree，基于 |W| 归一化后的阈值二值图）
% - 聚类系数（clustering coefficient，二值无向图）
% 并在 group-level 上做平均，画每个条件的网络图
% 额外输出：
% - ROI 分开的三线表：connect vs unconnect 的 paired t-test + 显著性标注
% - 不区分ROI的全局检验：先对ROI折叠（每被试一个值）再 paired t-test

clear; clc; close all;

%% ---------- 基本设置（按需修改） ----------
root_bids = '/home/xue/data/BIDS_prep';  % <--- 改成你的 BIDS_prep 根目录

% 参与 PPI 分析的被试列表（按你的实际情况修改）
subjects = {'sub-01','sub-03','sub-04','sub-05', ...
            'sub-06','sub-07','sub-09','sub-10'};

condLevels = {'connect','unconnect'};  % 要分析的条件（顺序无所谓）

% 图论阈值：对每个被试、每个条件，将 |W| 归一化到 [0,1] 再以此阈值生成二值图
thr_norm = 0.3;  % 可以之后改成 0.2 / 0.4 看鲁棒性

%% ---------- 预分配（等第一次 load 确定 ROI 数） ----------
nSub = numel(subjects);
baseROI_list = {};
nBase = [];

% strength(ROI × cond × subj), degree(ROI × cond × subj), clustering(ROI × cond × subj)
strength_all   = [];
degree_all     = [];
cluster_all    = [];

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
        baseROI_list = results_fix.baseROI_list;
        nBase = numel(baseROI_list);

        nCond = numel(condLevels);
        strength_all = nan(nBase, nCond, nSub);
        degree_all   = nan(nBase, nCond, nSub);
        cluster_all  = nan(nBase, nCond, nSub);

        fprintf('    [Info] 检测到 base ROI 数 = %d\n', nBase);
        disp(baseROI_list(:)');
    end

    for ci = 1:numel(condLevels)
        condLabel = condLevels{ci};

        if ~isfield(results_fix, condLabel) || ...
           ~isfield(results_fix.(condLabel), 'beta_PPI_bilat') || ...
           isempty(results_fix.(condLabel).beta_PPI_bilat)
            warning('    [Subj %s] 条件 %s 没有 beta_PPI_bilat，填 NaN。', subj, condLabel);
            continue;
        end

        W = results_fix.(condLabel).beta_PPI_bilat;   % [nBase × nBase]

        % 对称化：取 (W + W')/2，保留符号
        W_sym = (W + W.') / 2;

        % 主对角线置 0（不考虑自连）
        W_sym(1:nBase+1:end) = 0;

        % --- 计算图论指标 ---
        [node_strength, node_degree, node_cluster] = compute_graph_measures(W_sym, thr_norm);

        strength_all(:, ci, s) = node_strength;
        degree_all(:,   ci, s) = node_degree;
        cluster_all(:,  ci, s) = node_cluster;

        fprintf('    [Subj %s, Cond %s] strength/degree/clustering 已计算。\n', subj, condLabel);
    end
end

%% ---------- Group-level 统计：均值 ± SEM 并打印 ----------
fprintf('\n================ Group-level 图论指标 (mean ± SEM) ================\n');

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

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

    fprintf('\n===== 条件 %s：N = %d =====\n', condLabel, nValid);
    for r = 1:nBase
        fprintf('ROI %-6s | strength = %6.3f ± %5.3f | degree = %4.2f ± %4.2f | C = %5.3f ± %5.3f\n', ...
            baseROI_list{r}, ...
            m_strength(r), sem_strength(r), ...
            m_degree(r),   sem_degree(r), ...
            m_cluster(r),  sem_cluster(r));
    end
end

%% ---------- 条件差异三线表：paired t-test (connect vs unconnect) ----------
idxC = find(strcmp(condLevels, 'connect'));
idxU = find(strcmp(condLevels, 'unconnect'));
if isempty(idxC) || isempty(idxU)
    error('condLevels 中必须包含 connect 和 unconnect 才能做差异检验。');
end

out_dir = fullfile(root_bids, 'derivatives', 'ppi_roi_runwise', 'group_graph_metrics');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

T_strength = make_three_line_table_paired(strength_all, baseROI_list, idxC, idxU, 'strength');
T_degree   = make_three_line_table_paired(degree_all,   baseROI_list, idxC, idxU, 'degree');
T_cluster  = make_three_line_table_paired(cluster_all,  baseROI_list, idxC, idxU, 'clustering');

fprintf('\n================ 条件差异（paired t-test）：connect vs unconnect ================\n');
disp(T_strength);
disp(T_degree);
disp(T_cluster);

writetable(T_strength, fullfile(out_dir, 'threeLine_strength_connect_vs_unconnect.csv'));
writetable(T_degree,   fullfile(out_dir, 'threeLine_degree_connect_vs_unconnect.csv'));
writetable(T_cluster,  fullfile(out_dir, 'threeLine_clustering_connect_vs_unconnect.csv'));

T_all = [T_strength; T_degree; T_cluster];
writetable(T_all, fullfile(out_dir, 'threeLine_ALLmetrics_connect_vs_unconnect_long.csv'));
fprintf('[Saved] 三线表 CSV 已保存到: %s\n', out_dir);

%% ---------- 不区分ROI：先对ROI折叠(每被试一个值) -> paired t-test ----------
% mode 可选："mean" 或 "median"
Tg_strength = global_test_collapseROI(strength_all, idxC, idxU, 'strength', "mean");
Tg_degree   = global_test_collapseROI(degree_all,   idxC, idxU, 'degree',   "mean");
Tg_cluster  = global_test_collapseROI(cluster_all,  idxC, idxU, 'clustering',"mean");

fprintf('\n================ 不区分ROI（ROI折叠后）paired t-test：connect vs unconnect ================\n');
disp(Tg_strength);
disp(Tg_degree);
disp(Tg_cluster);

writetable([Tg_strength; Tg_degree; Tg_cluster], fullfile(out_dir,'GLOBAL_noROI_pairedT.csv'));

%% ---------- 可视化：条件层面的平均网络图 ----------
fprintf('\n================ 绘制每个条件的 group-level 网络图 ================\n');

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

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
        W_group(mask_nonNan)  = W_group(mask_nonNan) + W_sym(mask_nonNan);
        count_W(mask_nonNan)  = count_W(mask_nonNan) + 1;
    end

    W_mean = W_group ./ max(count_W, 1);
    W_mean(count_W == 0) = NaN;

    [node_strength_g, ~, ~, A_bin, W_abs_norm] = compute_graph_measures(W_mean, thr_norm);

    % 圆形布局
    theta = linspace(0, 2*pi, nBase+1);
    theta = theta(1:nBase);
    R = 1.0;
    x = R * cos(theta);
    y = R * sin(theta);

    % 节点大小 ~ group-level strength
    max_strength = max(node_strength_g);
    if isempty(max_strength) || max_strength <= 0
        node_size = 300 * ones(nBase,1);
    else
        node_size = 200 + 400 * (node_strength_g / max_strength); % [200, 600]
    end

    figure('Name', sprintf('PPI Graph - %s (group mean)', condLabel));
    set(gcf, 'Color', [1 1 1]);
    hold on; axis equal off;

    % 画边（上三角）
    for i = 1:nBase-1
        for j = i+1:nBase
            if isnan(W_mean(i,j)); continue; end
            if ~A_bin(i,j); continue; end

            wij = W_mean(i,j);

            % 颜色：正 β 用红色，负 β 用蓝色，0 用灰色
            if wij > 0
                edge_color = [1 0 0];
            elseif wij < 0
                edge_color = [0 0.4 1];
            else
                edge_color = [0.7 0.7 0.7];
            end

            % 线宽 ~ |W_mean| 归一化
            wij_norm = W_abs_norm(i,j);
            lw = 1 + 4 * wij_norm;  % [1,5]
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

%% =====================================================================
%% ========================= Local functions ============================
%% =====================================================================

function [node_strength, node_degree, node_cluster, A_bin, W_abs_norm] = compute_graph_measures(W_sym, thr_norm)
% 输入:
%   W_sym    : n×n 对称加权矩阵（允许 NaN；主对角线会被置 0）
%   thr_norm : 阈值（基于 |W| 归一化到 [0,1] 后）
% 输出:
%   node_strength : sum(|w_ij|)
%   node_degree   : 二值图下度数
%   node_cluster  : 二值无向图聚类系数
%   A_bin         : 阈值后的二值邻接矩阵
%   W_abs_norm    : 归一化后的 |W|

    n = size(W_sym,1);

    W_sym(1:n+1:end) = 0;

    W_abs = abs(W_sym);
    mask_valid = ~isnan(W_abs);

    max_w = max(W_abs(mask_valid));
    if isempty(max_w) || max_w == 0
        W_abs_norm = zeros(n);
    else
        W_abs_norm = W_abs / max_w;
    end

    A_bin = W_abs_norm > thr_norm;
    A_bin = triu(A_bin,1);
    A_bin = A_bin + A_bin.';
    A_bin(1:n+1:end) = 0;

    node_strength = sum(W_abs, 2, 'omitnan');
    node_degree   = sum(A_bin, 2);

    node_cluster = zeros(n,1);
    for i = 1:n
        neigh = find(A_bin(i,:));
        k = numel(neigh);
        if k < 2
            node_cluster(i) = 0;
            continue;
        end
        subA = A_bin(neigh, neigh);
        e = sum(subA(:)) / 2;
        node_cluster(i) = 2*e / (k*(k-1));
    end
end

function T = make_three_line_table_paired(metric_all, roi_list, idxC, idxU, metric_name)
% metric_all: [nROI x nCond x nSub]
% 输出：每 ROI 一行（paired t-test），并包含 mean±SEM、Δ、t/df/p、显著性标注

    nROI = size(metric_all,1);
    nSub = size(metric_all,3);

    ROI  = strings(nROI,1);
    Npair = nan(nROI,1);

    mC = nan(nROI,1); semC = nan(nROI,1);
    mU = nan(nROI,1); semU = nan(nROI,1);
    dMean = nan(nROI,1);

    tval = nan(nROI,1); df = nan(nROI,1); pval = nan(nROI,1);
    Sig  = strings(nROI,1);

    connect_meanSEM   = strings(nROI,1);
    unconnect_meanSEM = strings(nROI,1);

    for r = 1:nROI
        ROI(r) = string(roi_list{r});

        a = squeeze(metric_all(r, idxC, :)); % connect
        b = squeeze(metric_all(r, idxU, :)); % unconnect

        ok = isfinite(a) & isfinite(b);
        np = sum(ok);
        Npair(r) = np;

        if np < 2
            Sig(r) = "NA";
            connect_meanSEM(r)   = "NA";
            unconnect_meanSEM(r) = "NA";
            continue;
        end

        a = a(ok); b = b(ok);

        mC(r) = mean(a);
        mU(r) = mean(b);
        semC(r) = std(a,0) / sqrt(np);
        semU(r) = std(b,0) / sqrt(np);

        d = a - b;
        dMean(r) = mean(d);

        [~, p, ~, stats] = ttest(a, b);
        tval(r) = stats.tstat;
        df(r)   = stats.df;
        pval(r) = p;

        Sig(r) = sig_star(p);

        connect_meanSEM(r)   = sprintf('%.3f ± %.3f', mC(r), semC(r));
        unconnect_meanSEM(r) = sprintf('%.3f ± %.3f', mU(r), semU(r));
    end

    Metric = repmat(string(metric_name), nROI, 1);

    T = table(Metric, ROI, Npair, connect_meanSEM, unconnect_meanSEM, dMean, tval, df, pval, Sig, ...
        'VariableNames', {'Metric','ROI','N','connect_meanSEM','unconnect_meanSEM','Delta_CminusU','t','df','p','Sig'});
end

function T = global_test_collapseROI(metric_all, idxC, idxU, metric_name, mode)
% 不区分 ROI：每个被试先把 ROI 折叠成一个值，再做 paired t-test
% mode: "mean" 或 "median"

    nSub = size(metric_all,3);
    xC = nan(nSub,1);
    xU = nan(nSub,1);

    if ischar(mode); mode = string(mode); end

    for s = 1:nSub
        a = squeeze(metric_all(:, idxC, s)); % [nROI x 1]
        b = squeeze(metric_all(:, idxU, s));
        ok = isfinite(a) & isfinite(b);      % 两条件都有效的 ROI

        if sum(ok) < 1
            continue;
        end

        if mode == "median"
            xC(s) = median(a(ok));
            xU(s) = median(b(ok));
        else
            xC(s) = mean(a(ok));
            xU(s) = mean(b(ok));
        end
    end

    okS = isfinite(xC) & isfinite(xU);
    N   = sum(okS);

    if N < 2
        warning('%s: 有效被试数 < 2，无法做 t-test。', metric_name);
        T = table(string(metric_name), N, "NA", "NA", NaN, NaN, NaN, NaN, "NA", ...
            'VariableNames', {'Metric','N','connect_meanSEM','unconnect_meanSEM','Delta_CminusU','t','df','p','Sig'});
        return;
    end

    xC = xC(okS); xU = xU(okS);

    mC = mean(xC); mU = mean(xU);
    semC = std(xC,0) / sqrt(N);
    semU = std(xU,0) / sqrt(N);
    dMean = mean(xC - xU);

    [~, p, ~, stats] = ttest(xC, xU);
    tval = stats.tstat;
    df   = stats.df;

    Sig = sig_star(p);

    T = table(string(metric_name), N, ...
        string(sprintf('%.3f ± %.3f', mC, semC)), ...
        string(sprintf('%.3f ± %.3f', mU, semU)), ...
        dMean, tval, df, p, string(Sig), ...
        'VariableNames', {'Metric','N','connect_meanSEM','unconnect_meanSEM','Delta_CminusU','t','df','p','Sig'});
end

function s = sig_star(p)
% 常用显著性标注：*** <0.001, ** <0.01, * <0.05, ns >=0.05
% 注意：此函数在本脚本中只定义一次，避免“已在此作用域内声明”的报错
    if ~isfinite(p)
        s = "NA";
    elseif p < 0.001
        s = "***";
    elseif p < 0.01
        s = "**";
    elseif p < 0.05
        s = "*";
    else
        s = "ns";
    end
end
