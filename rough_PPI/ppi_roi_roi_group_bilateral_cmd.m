%% ppi_roi_roi_group_bilateral_cmd.m
% 组水平：基于 bilateral ROI-ROI PPI 的 group-level 分析
%
% 假定每个被试已有：
%   <root_bids>/derivatives/ppi_roi/<subj>/<subj>_ppi_roi_roi_bilateral_allconds.mat
%
% 该文件中包含变量：
%   - results_bilat.unconnect
%   - results_bilat.connect
%
% 每个条件下：
%   results_bilat.(cond).(seed).(target).beta_ppi  (被试内 L/R 已固定效应合并)
%
% 本脚本步骤：
%   Step 1) 对每个条件，找到所有被试共同的 bilateral ROI 集合
%   Step 2) 对每个条件、每条有向边(seed->target)，跨被试做 t 检验 (beta vs 0)
%   Step 3) 对每个条件构建对称 beta/p 矩阵、FDR 校正并可视化 group-level 网络
%   Step 4) 对 unconnect vs connect 做组水平差异检验 (within-subject Δbeta)，FDR+可视化

clear; clc;
fprintf('\n================ Group-level bilateral ROI-ROI PPI ================\n');

%% 0. 配置

% 被试列表（按需修改）
subj_list = { 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07','sub-09', 'sub-10'};

nSubj = numel(subj_list);

% 单被试 PPI 结果所在根目录（与单被试脚本 cfg.root_bids 保持一致）
root_bids = '/home/xue/data/BIDS_prep';

% group-level 输出目录
out_dir_group = fullfile(root_bids, 'derivatives', 'ppi_roi_group');
if ~exist(out_dir_group, 'dir')
    fprintf('  [Info] 创建 group 输出目录: %s\n', out_dir_group);
    mkdir(out_dir_group);
end

% 统计相关参数
alpha_unc = 0.05;   % 未校正显著性阈值
q_fdr     = 0.05;   % FDR 控制水平

fprintf('  [Config] nSubj = %d\n', nSubj);
fprintf('  [Config] root_bids = %s\n', root_bids);

%% 1. 第一次遍历被试：获取条件列表 & 每个条件共同的 bilateral ROI 集

fprintf('\n[Step 1] 统计所有被试共同的 bilateral ROI 集合...\n');

cond_levels = {};
roi_common  = struct();

for si = 1:nSubj
    subj = subj_list{si};
    fprintf('  [Subj %s] 读取 bilateral PPI 文件...\n', subj);

    subj_ppi_dir  = fullfile(root_bids, 'derivatives', 'ppi_roi', subj);
    subj_ppi_file = fullfile(subj_ppi_dir, ...
        sprintf('%s_ppi_roi_roi_bilateral_allconds.mat', subj));

    assert(exist(subj_ppi_file, 'file') == 2, ...
        '找不到被试 %s 的 bilateral PPI 文件: %s', subj, subj_ppi_file);

    S = load(subj_ppi_file, 'results_bilat');
    assert(isfield(S, 'results_bilat'), ...
        '文件 %s 中不包含 results_bilat', subj_ppi_file);

    res_bilat = S.results_bilat;

    if si == 1
        % 第一位被试时记录条件名和 ROI 列表
        cond_levels = fieldnames(res_bilat);
        fprintf('    条件列表: %s\n', strjoin(cond_levels', ', '));

        for ci = 1:numel(cond_levels)
            cond_label = cond_levels{ci};
            roi_names  = fieldnames(res_bilat.(cond_label));
            roi_common.(cond_label) = roi_names(:);  % 初始为该被试的所有 ROI
        end
    else
        % 从第二位被试开始，对每个条件取 ROI 名的交集
        for ci = 1:numel(cond_levels)
            cond_label = cond_levels{ci};
            if ~isfield(res_bilat, cond_label)
                error('被试 %s 的 results_bilat 中不存在条件 %s', ...
                    subj, cond_label);
            end
            roi_names_subj = fieldnames(res_bilat.(cond_label));
            roi_names_subj = roi_names_subj(:);

            roi_common.(cond_label) = intersect(roi_common.(cond_label), ...
                                                roi_names_subj);
        end
    end
end

% 打印共同 ROI 信息
for ci = 1:numel(cond_levels)
    cond_label   = cond_levels{ci};
    roi_list_cond = roi_common.(cond_label);
    fprintf('  [Info] 条件 %s：共同 bilateral ROI 数 = %d，ROI = %s\n', ...
        cond_label, numel(roi_list_cond), strjoin(roi_list_cond', ', '));
end

%% 2. 第二次遍历被试：按条件构建 beta_all_subj (subj x ROI x ROI)，并做组统计

fprintf('\n[Step 2] 构建组水平 beta 矩阵并进行单条件 t 检验 (beta vs 0)...\n');

group_stats = struct();

for ci = 1:numel(cond_levels)
    cond_label = cond_levels{ci};
    roi_list   = roi_common.(cond_label);
    nROI       = numel(roi_list);

    if nROI < 2
        warning('  [Warning] 条件 %s 的共同 ROI 数 < 2，跳过该条件的组分析。', cond_label);
        continue;
    end

    fprintf('\n  === 条件 %s ===\n', cond_label);
    fprintf('      使用 %d 个共同 bilateral ROI: %s\n', ...
        nROI, strjoin(roi_list', ', '));

    % beta_all_subj: [nSubj x nROI x nROI]，对角线保持 NaN
    beta_all_subj = NaN(nSubj, nROI, nROI);

    % 再次遍历被试，填充 beta_all_subj
    for si = 1:nSubj
        subj = subj_list{si};
        subj_ppi_dir  = fullfile(root_bids, 'derivatives', 'ppi_roi', subj);
        subj_ppi_file = fullfile(subj_ppi_dir, ...
            sprintf('%s_ppi_roi_roi_bilateral_allconds.mat', subj));

        S = load(subj_ppi_file, 'results_bilat');
        res_bilat = S.results_bilat;

        assert(isfield(res_bilat, cond_label), ...
            '被试 %s 的 results_bilat 中不存在条件 %s', subj, cond_label);

        res_cond = res_bilat.(cond_label);

        % 遍历 ROI 作为 seed/target
        for i = 1:nROI
            seed_name = roi_list{i};
            if ~isfield(res_cond, seed_name), continue; end

            res_seed = res_cond.(seed_name);

            for j = 1:nROI
                if i == j, continue; end

                target_name = roi_list{j};
                if ~isfield(res_seed, target_name), continue; end

                rr = res_seed.(target_name);
                if isfield(rr, 'beta_ppi') && ~isnan(rr.beta_ppi)
                    beta_all_subj(si, i, j) = rr.beta_ppi;
                end
            end
        end
    end

    % 组水平：对每条有向边 beta_all_subj(:, i, j) 做 t 检验 (vs 0)
    beta_mean_dir = NaN(nROI, nROI);
    t_dir         = NaN(nROI, nROI);
    p_dir         = NaN(nROI, nROI);
    n_dir         = NaN(nROI, nROI);

    fprintf('      对每条有向边 (seed->target) 做组水平 t 检验...\n');

    for i = 1:nROI
        for j = 1:nROI
            if i == j, continue; end

            b_vec = squeeze(beta_all_subj(:, i, j));   % [nSubj x 1]
            valid = ~isnan(b_vec);
            n     = sum(valid);

            if n < 2
                % 样本数太少，无法做 t 检验
                continue;
            end

            b_valid = b_vec(valid);
            m  = mean(b_valid);
            sd = std(b_valid, 0);    % N-1 归一化
            if sd == 0
                beta_mean_dir(i,j) = m;
                t_dir(i,j)         = Inf;
                p_dir(i,j)         = 0;
                n_dir(i,j)         = n;
            else
                se = sd / sqrt(n);
                t  = m / se;
                df = n - 1;
                p  = 2 * tcdf(-abs(t), df);

                beta_mean_dir(i,j) = m;
                t_dir(i,j)         = t;
                p_dir(i,j)         = p;
                n_dir(i,j)         = n;
            end
        end
    end

    % 构建对称 beta/p 矩阵用于可视化：双向平均 beta，p 取最小值
    beta_sym = NaN(nROI, nROI);
    p_sym    = NaN(nROI, nROI);

    for i = 1:nROI
        for j = i+1:nROI
            pair_beta = [beta_mean_dir(i,j), beta_mean_dir(j,i)];
            pair_p    = [p_dir(i,j),         p_dir(j,i)];

            valid_beta = pair_beta(~isnan(pair_beta));
            valid_p    = pair_p(~isnan(pair_p));

            if ~isempty(valid_beta)
                beta_sym(i,j) = mean(valid_beta);
                beta_sym(j,i) = beta_sym(i,j);
            end
            if ~isempty(valid_p)
                p_sym(i,j) = min(valid_p);
                p_sym(j,i) = p_sym(i,j);
            end
        end
        beta_sym(i,i) = 0;
        p_sym(i,i)    = NaN;
    end

    % ---------- FDR 校正（对称 p_sym 的上三角） ----------
    mask_ut = triu(~isnan(p_sym), 1);
    p_vec   = p_sym(mask_ut);
    m_edges = numel(p_vec);

    sig_fdr_mat = false(size(p_sym));
    p_crit = NaN;

    if m_edges > 0
        [h_fdr_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_mat(mask_ut) = h_fdr_vec;
        sig_fdr_mat = sig_fdr_mat | sig_fdr_mat.';  % 对称化

        fprintf('      条件 %s：FDR q=%.3f，%d/%d 条边通过 FDR (p<=%.3g)\n', ...
            cond_label, q_fdr, sum(h_fdr_vec), m_edges, p_crit);
    else
        fprintf('      条件 %s：没有可用于 FDR 的边。\n', cond_label);
    end

    % 未校正显著性矩阵
    sig_unc_mat = (p_sym < alpha_unc) & ~isnan(p_sym);
    sig_unc_mat = sig_unc_mat | sig_unc_mat.';

    % 保存组统计结果
    gs = struct();
    gs.cond_label    = cond_label;
    gs.roi_list      = roi_list;
    gs.beta_all_subj = beta_all_subj;
    gs.beta_mean_dir = beta_mean_dir;
    gs.t_dir         = t_dir;
    gs.p_dir         = p_dir;
    gs.n_dir         = n_dir;
    gs.beta_sym      = beta_sym;
    gs.p_sym         = p_sym;
    gs.sig_fdr_mat   = sig_fdr_mat;
    gs.sig_unc_mat   = sig_unc_mat;
    gs.alpha_unc     = alpha_unc;
    gs.q_fdr         = q_fdr;
    gs.p_crit_fdr    = p_crit;

    group_stats.(cond_label) = gs;
end

% 保存 group-level 统计结果
out_group_file = fullfile(out_dir_group, ...
    'ppi_roi_roi_group_bilateral_singlecond.mat');
save(out_group_file, 'subj_list', 'group_stats', 'alpha_unc', 'q_fdr', '-v7.3');
fprintf('\n[Step 2] 单条件组统计结果已保存到: %s\n', out_group_file);

%% 3. group-level 脑网络可视化（每个条件单独绘图）

fprintf('\n[Step 3] 可视化每个条件的 group-level bilateral PPI 网络...\n');

% 颜色定义（与单被试版本一致）
col_red    = [1   0     0  ];  % FDR 正向显著
col_blue   = [0   0.6   1  ];  % FDR 负向显著
col_orange = [1   0.5   0  ];  % 未校正正向显著但 FDR 不显著
col_green  = [0   1     0  ];  % 未校正负向显著但 FDR 不显著
col_white  = [1   1     1  ];

cond_levels_group = fieldnames(group_stats);

for ci = 1:numel(cond_levels_group)
    cond_label = cond_levels_group{ci};
    gs = group_stats.(cond_label);

    roi_all     = gs.roi_list;
    nROI        = numel(roi_all);
    beta_sym    = gs.beta_sym;
    p_sym       = gs.p_sym;
    sig_fdr_mat = gs.sig_fdr_mat;
    sig_unc_mat = gs.sig_unc_mat;

    if nROI < 2
        warning('  [Warning] 条件 %s 的 ROI 数 < 2，跳过绘图。', cond_label);
        continue;
    end

    fprintf('  [Network] 条件 %s：绘制 group-level bilateral PPI 网络...\n', cond_label);

    % 节点坐标：均匀放在一个圆上
    theta = linspace(0, 2*pi, nROI+1);
    theta = theta(1:nROI);
    R     = 1;

    x = R * cos(theta);
    y = R * sin(theta);

    % 计算最大 |beta| 用于设置线宽
    valid_abs_beta = abs(beta_sym(~isnan(beta_sym)));
    max_abs_beta   = max(valid_abs_beta);
    if isempty(max_abs_beta) || max_abs_beta == 0
        max_abs_beta = 1;
    end

    % 画图
    figure('Name', sprintf('Group - %s PPI network (bilateral, FDR)', cond_label));
    set(gcf, 'Color', [0 0 0]);
    ax = gca;
    set(ax, 'Color', [0 0 0]);
    hold on; axis equal off;

    % 画边（只画上三角，避免重复）
    for i = 1:nROI-1
        for j = i+1:nROI
            if isnan(beta_sym(i,j))
                continue;
            end

            b_ij = beta_sym(i,j);
            p_ij = p_sym(i,j);

            is_fdr = sig_fdr_mat(i,j);
            is_unc = sig_unc_mat(i,j);

            edge_color = col_white;

            if is_fdr
                % FDR 显著
                if b_ij > 0
                    edge_color = col_red;
                elseif b_ij < 0
                    edge_color = col_blue;
                end
            elseif is_unc
                % 未校正显著但 FDR 不显著
                if b_ij > 0
                    edge_color = col_orange;
                elseif b_ij < 0
                    edge_color = col_green;
                end
            else
                edge_color = col_white;
            end

            % 线宽随 |beta| 变化
            w = 1.5 + 3 * (abs(b_ij) / max_abs_beta);

            % 画边
            line([x(i) x(j)], [y(i) y(j)], ...
                'Color', edge_color, 'LineWidth', w);

            % 在线中点标注 β 值
            xm = (x(i) + x(j)) / 2;
            ym = (y(i) + y(j)) / 2;
            text(xm, ym, sprintf('%.3f', b_ij), ...
                'Color', edge_color, ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'none');
        end
    end

    % 画节点
    node_size = 300;
    for i = 1:nROI
        scatter(x(i), y(i), node_size, ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5);
        text(x(i)*1.2, y(i)*1.2, roi_all{i}, ...
            'Color', [1 1 1], ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end

    title(sprintf('Group: %s (bilateral PPI, FDR)', cond_label), ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold');
end

%% 4. 条件差异：unconnect vs connect 的组水平 Δbeta 分析＋网络图

fprintf('\n[Step 4] 组水平条件差异分析：unconnect vs connect (Δbeta = uncon - con)...\n');

% 检查条件是否包含 unconnect / connect
conds_avail = fieldnames(group_stats);
assert(ismember('unconnect', conds_avail) && ismember('connect', conds_avail), ...
    'group_stats 中必须同时包含条件 unconnect 和 connect 才能做差异分析。');

gs_un = group_stats.unconnect;
gs_co = group_stats.connect;

% 共同 ROI 集合，用于差异分析
roi_un = gs_un.roi_list;
roi_co = gs_co.roi_list;
roi_diff = intersect(roi_un, roi_co);
nROI_diff = numel(roi_diff);

fprintf('  [Diff] 共同 bilateral ROI 数 = %d，ROI = %s\n', ...
    nROI_diff, strjoin(roi_diff', ', '));

if nROI_diff < 2
    warning('  [Diff] 共同 ROI 数 < 2，无法进行条件差异网络分析。');
else
    % 在各自条件中的索引
    [~, idx_un] = ismember(roi_diff, roi_un);
    [~, idx_co] = ismember(roi_diff, roi_co);

    if any(idx_un == 0) || any(idx_co == 0)
        error('  [Diff] ROI 映射出现问题，请检查 roi_list。');
    end

    % 拿出两个条件的 beta_all_subj，只保留共同 ROI 子集
    beta_un_all = gs_un.beta_all_subj(:, idx_un, idx_un);  % [nSubj x nROI_diff x nROI_diff]
    beta_co_all = gs_co.beta_all_subj(:, idx_co, idx_co);  % 同维度

    % 对每条有向边构建 Δbeta(subj, i, j) = beta_un - beta_co
    delta_all_subj = NaN(nSubj, nROI_diff, nROI_diff);
    delta_mean_dir = NaN(nROI_diff, nROI_diff);
    t_diff_dir     = NaN(nROI_diff, nROI_diff);
    p_diff_dir     = NaN(nROI_diff, nROI_diff);
    n_diff_dir     = NaN(nROI_diff, nROI_diff);

    fprintf('  [Diff] 对每条有向边 (seed->target) 做组内配对 t 检验...\n');

    for i = 1:nROI_diff
        for j = 1:nROI_diff
            if i == j, continue; end

            b_un = squeeze(beta_un_all(:, i, j));  % [nSubj x 1]
            b_co = squeeze(beta_co_all(:, i, j));  % [nSubj x 1]

            valid = ~isnan(b_un) & ~isnan(b_co);
            n     = sum(valid);
            if n < 2
                continue;
            end

            b_un_v = b_un(valid);
            b_co_v = b_co(valid);
            d_vec  = b_un_v - b_co_v;   % Δbeta = uncon - con

            delta_all_subj(valid, i, j) = d_vec;

            m  = mean(d_vec);
            sd = std(d_vec, 0);

            if sd == 0
                delta_mean_dir(i,j) = m;
                t_diff_dir(i,j)     = Inf;
                p_diff_dir(i,j)     = 0;
                n_diff_dir(i,j)     = n;
            else
                se = sd / sqrt(n);
                t  = m / se;
                df = n - 1;
                p  = 2 * tcdf(-abs(t), df);

                delta_mean_dir(i,j) = m;
                t_diff_dir(i,j)     = t;
                p_diff_dir(i,j)     = p;
                n_diff_dir(i,j)     = n;
            end
        end
    end

    % 构建对称 Δbeta / p 矩阵：双向平均 Δbeta，p 取最小值
    delta_sym = NaN(nROI_diff, nROI_diff);
    p_sym_diff = NaN(nROI_diff, nROI_diff);

    for i = 1:nROI_diff
        for j = i+1:nROI_diff
            pair_delta = [delta_mean_dir(i,j), delta_mean_dir(j,i)];
            pair_p     = [p_diff_dir(i,j),     p_diff_dir(j,i)];

            valid_d = pair_delta(~isnan(pair_delta));
            valid_p = pair_p(~isnan(pair_p));

            if ~isempty(valid_d)
                delta_sym(i,j) = mean(valid_d);
                delta_sym(j,i) = delta_sym(i,j);
            end
            if ~isempty(valid_p)
                p_sym_diff(i,j) = min(valid_p);
                p_sym_diff(j,i) = p_sym_diff(i,j);
            end
        end
        delta_sym(i,i)  = 0;
        p_sym_diff(i,i) = NaN;
    end

    % FDR (Δbeta 的 p 值)
    mask_ut = triu(~isnan(p_sym_diff), 1);
    p_vec   = p_sym_diff(mask_ut);
    m_edges = numel(p_vec);

    sig_fdr_diff = false(size(p_sym_diff));
    p_crit_diff  = NaN;

    if m_edges > 0
        [h_fdr_vec, p_crit_diff] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_diff(mask_ut) = h_fdr_vec;
        sig_fdr_diff = sig_fdr_diff | sig_fdr_diff.';   % 对称化

        fprintf('  [Diff] FDR q=%.3f，%d/%d 条边在 unconnect vs connect 间存在显著差异 (p<=%.3g)\n', ...
            q_fdr, sum(h_fdr_vec), m_edges, p_crit_diff);
    else
        fprintf('  [Diff] 无可用于差异 FDR 的边。\n');
    end

    % 未校正显著
    sig_unc_diff = (p_sym_diff < alpha_unc) & ~isnan(p_sym_diff);
    sig_unc_diff = sig_unc_diff | sig_unc_diff.';

    % 保存差异结果
    diff_stats = struct();
    diff_stats.cond_pair      = {'unconnect','connect'};
    diff_stats.roi_list       = roi_diff;
    diff_stats.delta_all_subj = delta_all_subj;
    diff_stats.delta_mean_dir = delta_mean_dir;
    diff_stats.t_diff_dir     = t_diff_dir;
    diff_stats.p_diff_dir     = p_diff_dir;
    diff_stats.n_diff_dir     = n_diff_dir;
    diff_stats.delta_sym      = delta_sym;
    diff_stats.p_sym_diff     = p_sym_diff;
    diff_stats.sig_fdr_diff   = sig_fdr_diff;
    diff_stats.sig_unc_diff   = sig_unc_diff;
    diff_stats.alpha_unc      = alpha_unc;
    diff_stats.q_fdr          = q_fdr;
    diff_stats.p_crit_diff    = p_crit_diff;

    % 保存 group + diff 结果
    out_group_file2 = fullfile(out_dir_group, ...
        'ppi_roi_roi_group_bilateral_with_diff.mat');
    save(out_group_file2, 'subj_list', 'group_stats', 'diff_stats', ...
        'alpha_unc', 'q_fdr', '-v7.3');
    fprintf('  [Diff] 组差异结果已保存到: %s\n', out_group_file2);

    % 绘制差异网络图：Δbeta = uncon - con
    fprintf('  [Diff] 绘制 group-level 差异网络图 (Δbeta = uncon - con)...\n');

    roi_all = roi_diff;
    nROI    = nROI_diff;

    % 节点坐标
    theta = linspace(0, 2*pi, nROI+1);
    theta = theta(1:nROI);
    R     = 1;
    x = R * cos(theta);
    y = R * sin(theta);

    % 最大 |Δbeta|
    valid_abs_delta = abs(delta_sym(~isnan(delta_sym)));
    max_abs_delta   = max(valid_abs_delta);
    if isempty(max_abs_delta) || max_abs_delta == 0
        max_abs_delta = 1;
    end

    figure('Name', 'Group diff (unconnect - connect) PPI network (bilateral, FDR)');
    set(gcf, 'Color', [0 0 0]);
    ax = gca;
    set(ax, 'Color', [0 0 0]);
    hold on; axis equal off;

    for i = 1:nROI-1
        for j = i+1:nROI
            if isnan(delta_sym(i,j))
                continue;
            end

            d_ij = delta_sym(i,j);       % mean Δbeta
            p_ij = p_sym_diff(i,j);
            is_fdr = sig_fdr_diff(i,j);
            is_unc = sig_unc_diff(i,j);

            % 颜色逻辑：
            %   Δbeta>0: unconnect > connect
            %   Δbeta<0: connect > unconnect
            edge_color = col_white;

            if is_fdr
                if d_ij > 0
                    edge_color = col_red;   % FDR 显著，uncon > con
                elseif d_ij < 0
                    edge_color = col_blue;  % FDR 显著，con > uncon
                end
            elseif is_unc
                if d_ij > 0
                    edge_color = col_orange; % 未校正显著，uncon > con
                elseif d_ij < 0
                    edge_color = col_green;  % 未校正显著，con > uncon
                end
            else
                edge_color = col_white;
            end

            % 线宽随 |Δbeta| 变化
            w = 1.5 + 3 * (abs(d_ij) / max_abs_delta);

            % 画边
            line([x(i) x(j)], [y(i) y(j)], ...
                'Color', edge_color, 'LineWidth', w);

            % 在线中点标注 Δbeta 值
            xm = (x(i) + x(j)) / 2;
            ym = (y(i) + y(j)) / 2;
            text(xm, ym, sprintf('%.3f', d_ij), ...
                'Color', edge_color, ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'none');
        end
    end

    % 画节点
    node_size = 300;
    for i = 1:nROI
        scatter(x(i), y(i), node_size, ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5);
        text(x(i)*1.2, y(i)*1.2, roi_all{i}, ...
            'Color', [1 1 1], ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end

    title('Group diff: unconnect - connect (bilateral PPI, FDR)', ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('\n================ Group-level PPI 分析结束 ================\n');

%% 本地 FDR 函数（Benjamini-Hochberg）
function [h, p_crit] = fdr_bh_local(p_vals, q)
    % p_vals: 列向量或行向量的 p 值（可含 NaN）
    % q     : FDR 控制水平，例如 0.05
    %
    % 返回：
    %   h      : 与 p_vals 同维度的逻辑向量，true 表示通过 FDR
    %   p_crit : 阈值 p*，即最大的 p，使得 p(k) <= (k/m)*q

    p = p_vals(:);
    h = false(size(p));

    valid = ~isnan(p);
    pv    = p(valid);
    if isempty(pv)
        p_crit = NaN;
        h(~valid) = false;
        h = reshape(h, size(p_vals));
        return;
    end

    [pv_sorted, ~] = sort(pv);
    m_valid = numel(pv_sorted);
    crit    = (1:m_valid)' * (q / m_valid);

    cmp = pv_sorted <= crit;
    if any(cmp)
        k_max   = find(cmp, 1, 'last');
        p_crit  = pv_sorted(k_max);
        sig_idx = p <= p_crit & valid;
        h(sig_idx) = true;
    else
        p_crit = 0;
    end

    h = reshape(h, size(p_vals));
end

