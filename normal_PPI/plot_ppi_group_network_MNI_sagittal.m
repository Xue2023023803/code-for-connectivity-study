%% ppi_group_network_MNI_sagittal.m
% Group-level ROI–ROI PPI 网络，可视化在 SPM avg152T1 的矢状切片上
% - 使用 Harvey & Dumoulin 2017 提供的 NTO / NPO / NPC1-3 / NF 的 MNI 坐标
% - 对每个条件 (connect, unconnect) 做 group-level t 检验
% - 额外画一个 connect − unconnect 的差异网络
% - 边颜色和阈值逻辑与你之前的 FDR 可视化保持一致

clear; clc;

%% 1. 基本配置：被试列表和数据路径
subj_list = { 'sub-01','sub-02','sub-03', 'sub-04','sub-05', 'sub-06', 'sub-07','sub-08','sub-09', 'sub-10'};

nSub = numel(subj_list);

data_root = '/home/xue/data/BIDS_prep/derivatives/ppi_roi_runwise';
condLevels = {'connect','unconnect'};   % 只关心这两个条件

%% 2. 读入每个被试的 FFX 结果，收集 beta_PPI_bilat
beta_all = struct();   % beta_all.(cond) : [nBase × nBase × nSub]
baseROI_list = [];
for ci = 1:numel(condLevels)
    beta_all.(condLevels{ci}) = [];  % 先占位
end

for s = 1:nSub
    subj = subj_list{s};
    ffx_file = fullfile(data_root, subj, ...
        sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral.mat', subj));
    assert(exist(ffx_file,'file')==2, '找不到 FFX 文件: %s', ffx_file);

    S = load(ffx_file, 'results_fix');
    results_fix = S.results_fix;

    if s == 1
        baseROI_list = results_fix.baseROI_list;
        nBase = numel(baseROI_list);

        for ci = 1:numel(condLevels)
            condLabel = condLevels{ci};
            beta_all.(condLabel) = nan(nBase, nBase, nSub);
        end
    else
        % 简单检查 baseROI_list 是否一致
        assert(isequal(baseROI_list, results_fix.baseROI_list), ...
            '被试 %s 的 baseROI_list 与第一个被试不一致', subj);
    end

    for ci = 1:numel(condLevels)
        condLabel = condLevels{ci};
        if isfield(results_fix, condLabel) && ...
           isfield(results_fix.(condLabel), 'beta_PPI_bilat')
            beta_all.(condLabel)(:, :, s) = results_fix.(condLabel).beta_PPI_bilat;
        else
            warning('被试 %s 缺少条件 %s 的 beta_PPI_bilat，保持 NaN', subj, condLabel);
        end
    end
end

fprintf('>>> 已加载 %d 名被试, base ROI 数 = %d\n', nSub, nBase);
disp('    baseROI_list:');
disp(baseROI_list(:)');

%% 3. 对每个条件做 group-level t 检验
%    逻辑：先在被试内对 i→j / j→i 取平均，得到每个被试在 {i,j} 上的“无向连接强度”，
%    然后对 N 个被试的平均值做 one-sample t 检验 vs 0，再做 FDR。
alpha_unc = 0.05;
q_fdr     = 0.05;

stats = struct();

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};
    beta_sub  = beta_all.(condLabel);       % [nBase × nBase × nSub]

    beta_sym  = nan(nBase, nBase);  % 存每条无向边的 group-level β 均值
    p_sym     = nan(nBase, nBase);  % 存对应的 p 值

    for i = 1:nBase
        for j = i+1:nBase
            % 取出两个方向
            vals_ij = squeeze(beta_sub(i, j, :));  % nSub×1
            vals_ji = squeeze(beta_sub(j, i, :));  % nSub×1

            % 被试内方向平均：对每个被试，把 i→j / j→i 在被试内取平均
            m_sub = nan(nSub, 1);
            for s = 1:nSub
                tmp = [vals_ij(s), vals_ji(s)];
                tmp = tmp(~isnan(tmp));  % 去掉 NaN（有一边缺失时用另一边）
                if ~isempty(tmp)
                    m_sub(s) = mean(tmp);
                else
                    m_sub(s) = NaN;
                end
            end

            valid = ~isnan(m_sub);
            n_valid = sum(valid);

            if n_valid >= 2
                [~, p] = ttest(m_sub(valid), 0);
                bmean  = mean(m_sub(valid));
            elseif n_valid == 1
                % 只有一个有效被试，p 不具统计意义，给一个非显著值
                p     = 1;
                bmean = m_sub(valid);
            else
                p     = NaN;
                bmean = NaN;
            end

            beta_sym(i, j) = bmean;
            beta_sym(j, i) = bmean;
            p_sym(i, j)    = p;
            p_sym(j, i)    = p;
        end
        beta_sym(i, i) = 0;
        p_sym(i, i)    = NaN;
    end

    % FDR（只对上三角做）
    mask_ut = triu(~isnan(p_sym), 1);
    p_vec   = p_sym(mask_ut);
    sig_fdr = false(size(p_sym));
    p_crit  = NaN;

    if ~isempty(p_vec)
        [h_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr(mask_ut) = h_vec;
        sig_fdr = sig_fdr | sig_fdr.';  % 对称化
        fprintf('条件 %s: FDR(q=%.3f) 通过边数 = %d / %d (p_crit = %.3g)\n', ...
            condLabel, q_fdr, sum(h_vec), numel(h_vec), p_crit);
    else
        fprintf('条件 %s: 无有效 p 值用于 FDR\n', condLabel);
    end

    stats.(condLabel).beta_sym = beta_sym;
    stats.(condLabel).p_sym    = p_sym;
    stats.(condLabel).sig_fdr  = sig_fdr;
    stats.(condLabel).p_crit   = p_crit;
end

%% 4. 计算 connect − unconnect 的差异网络
%    同样先在被试内做方向平均，然后对 N 个被试的差异值做 one-sample t 检验 vs 0
if all(isfield(beta_all, {'connect','unconnect'}))
    beta_diff_sub = beta_all.connect - beta_all.unconnect;  % [nBase × nBase × nSub]

    beta_diff_sym = nan(nBase, nBase);
    p_diff_sym    = nan(nBase, nBase);

    for i = 1:nBase
        for j = i+1:nBase
            vals_ij = squeeze(beta_diff_sub(i, j, :));  % nSub×1
            vals_ji = squeeze(beta_diff_sub(j, i, :));  % nSub×1

            % 被试内方向平均：每个被试在 {i,j} 上的差异 = (i→j 与 j→i 的平均差异)
            m_sub = nan(nSub, 1);
            for s = 1:nSub
                tmp = [vals_ij(s), vals_ji(s)];
                tmp = tmp(~isnan(tmp));
                if ~isempty(tmp)
                    m_sub(s) = mean(tmp);
                else
                    m_sub(s) = NaN;
                end
            end

            valid = ~isnan(m_sub);
            n_valid = sum(valid);

            if n_valid >= 2
                [~, p] = ttest(m_sub(valid), 0);   % 差异 vs 0
                bmean  = mean(m_sub(valid));
            elseif n_valid == 1
                p     = 1;
                bmean = m_sub(valid);
            else
                p     = NaN;
                bmean = NaN;
            end

            beta_diff_sym(i, j) = bmean;
            beta_diff_sym(j, i) = bmean;
            p_diff_sym(i, j)    = p;
            p_diff_sym(j, i)    = p;
        end
        beta_diff_sym(i, i) = 0;
        p_diff_sym(i, i)    = NaN;
    end

    % FDR for difference
    mask_ut = triu(~isnan(p_diff_sym), 1);
    p_vec   = p_diff_sym(mask_ut);
    sig_fdr_diff = false(size(p_diff_sym));
    p_crit_diff  = NaN;

    if ~isempty(p_vec)
        [h_vec, p_crit_diff] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_diff(mask_ut) = h_vec;
        sig_fdr_diff = sig_fdr_diff | sig_fdr_diff.';
        fprintf('差异 (connect − unconnect): FDR(q=%.3f) 通过边数 = %d / %d (p_crit = %.3g)\n', ...
            q_fdr, sum(h_vec), numel(h_vec), p_crit_diff);
    else
        fprintf('差异 (connect − unconnect): 无有效 p 值用于 FDR\n');
    end

    stats.diff.beta_sym = beta_diff_sym;
    stats.diff.p_sym    = p_diff_sym;
    stats.diff.sig_fdr  = sig_fdr_diff;
    stats.diff.p_crit   = p_crit_diff;
else
    warning('beta_all 中缺少 connect/unconnect，无法计算差异网络');
end

%% 5. 可视化：在 SPM avg152T1 的矢状切片上画 3 张图
% 颜色约定：
%   FDR 显著: 正向 = 红色, 负向 = 蓝色
%   未经 FDR 但 p<0.05: 正向 = 橙色, 负向 = 绿色
%   其余 = 白色细线
%
% 差异图中：
%   β_diff > 0 表示 connect 的 PPI > unconnect
%   β_diff < 0 表示 unconnect 的 PPI > connect

% 找到 SPM canonical 模板路径
if exist('spm','file') ~= 2
    error('找不到 SPM 函数 spm，请确认 SPM 已加入 MATLAB 路径');
end
spm_dir = spm('Dir');
template_nii = fullfile(spm_dir, 'canonical', 'avg152T1.nii');
assert(exist(template_nii,'file')==2, '找不到 avg152T1.nii: %s', template_nii);

% 加载模板
V = spm_vol(template_nii);
Y = spm_read_vols(V);     % 3D: [Nx × Ny × Nz]
dim = V.dim;

% 选一个中线附近的矢状切片 (固定 x index)，也可以改成你想要的 x
ix = round(dim(1)/2);     % 中线附近

% 提前计算每个 baseROI 的 (j,k) 像素坐标（MNI y,z → 模板 voxel）
node_xy = compute_roi_pixel_positions(baseROI_list, V);

% 绘图（connect, unconnect, connect−unconnect）
plot_ppi_mni_network(V, Y, ix, baseROI_list, node_xy, stats.connect, ...
    'Connect: ROI–ROI PPI (numerosity vs baseline)', alpha_unc);

plot_ppi_mni_network(V, Y, ix, baseROI_list, node_xy, stats.unconnect, ...
    'Unconnect: ROI–ROI PPI (numerosity vs baseline)', alpha_unc);

if isfield(stats, 'diff')
    plot_ppi_mni_network(V, Y, ix, baseROI_list, node_xy, stats.diff, ...
        'Connect − Unconnect: difference in PPI (β_{conn} − β_{unconn})', alpha_unc);
end

fprintf('\n>>> Group-level MNI 矢状切片网络可视化完成。\n');

%% ========== 子函数 1: MNI → 模板像素坐标 (只用 y,z) ==========
function node_xy = compute_roi_pixel_positions(baseROI_list, V)
% 使用 Harvey & Dumoulin 2017 中的 MNI 坐标（左右平均）来给节点放置位置
% 文献中给出的 (x,y,z) mm，如：
%   NTO 右:  44, -75, -4; 左: -42, -77, -3
%   NPO 右:  25, -82, 34; 左: -23, -80, 32
%   NPC1右:  22, -61, 60; 左: -22, -59, 61
%   NPC2右:  33, -40, 52; 左: -38, -43, 48
%   NPC3右:  45, -30, 40; 左: -48, -29, 34
%   NF  右:  24, -11, 52; 左: -22, -11, 50
%
% 我们对左右 y,z 取平均，只保留 y,z，用来表示拓扑上的相对位置。

harvey_yz = struct();

% NTO
harvey_yz.NTO  = [mean([-75, -77]), mean([-4, -3])];    % [y, z]
% NPO
harvey_yz.NPO  = [mean([-82, -80]), mean([34, 32])];
% NPC1
harvey_yz.NPC1 = [mean([-61, -59]), mean([60, 61])];
% NPC2
harvey_yz.NPC2 = [mean([-40, -43]), mean([52, 48])];
% NPC3
harvey_yz.NPC3 = [mean([-30, -29]), mean([40, 34])];
% NF
harvey_yz.NF   = [mean([-11, -11]), mean([52, 50])];

nBase = numel(baseROI_list);
node_xy = nan(nBase, 2);

% 从 affine 矩阵中取 y,z 的线性变换参数
M = V.mat;
a_y = M(2, 2); b_y = M(2, 4);     % y_world = a_y * j + b_y
a_z = M(3, 3); b_z = M(3, 4);     % z_world = a_z * k + b_z

for b = 1:nBase
    roi = baseROI_list{b};
    if isfield(harvey_yz, roi)
        yz = harvey_yz.(roi);      % [y_world, z_world]
        y_w = yz(1);
        z_w = yz(2);

        % world → voxel index（注意这里是浮点，可以用于 scatter）
        j = (y_w - b_y) / a_y;
        k = (z_w - b_z) / a_z;

        node_xy(b, :) = [j, k];
    else
        % 如果有其他 ROI 名（不在 6 个 numerosity map 里），先暂时留 NaN，
        % 后面绘图阶段会 fallback 到一个环形布局。
        node_xy(b, :) = [NaN, NaN];
    end
end

end

%% ========== 子函数 2: 在 MNI 模板矢状切片上画网络 ==========
function plot_ppi_mni_network(V, Y, ix, baseROI_list, node_xy, stat_struct, ...
    fig_title, alpha_unc)

beta_sym = stat_struct.beta_sym;
p_sym    = stat_struct.p_sym;
sig_fdr  = stat_struct.sig_fdr;

nBase = numel(baseROI_list);
dim   = V.dim;

% 提取指定矢状切片: [Ny × Nz]
slice = squeeze(Y(ix, :, :));   % x 固定，y,z 变化

% 计算 |β| 的最大值用于归一化线宽
mask_beta = ~isnan(beta_sym) & triu(true(size(beta_sym)), 1);
valid_abs = abs(beta_sym(mask_beta));
if isempty(valid_abs)
    max_abs_beta = 1;
else
    max_abs_beta = max(valid_abs);
    if max_abs_beta == 0
        max_abs_beta = 1;
    end
end

% FDR & 未校正显著 (p < alpha_unc)
mask_ut = triu(~isnan(p_sym), 1);
p_vec   = p_sym(mask_ut);
sig_unc = (p_sym < alpha_unc) & ~isnan(p_sym);
sig_unc = sig_unc | sig_unc.';

% 颜色定义
col_red    = [1.0, 0.0, 0.0];    % FDR 正向
col_blue   = [0.0, 0.6, 1.0];    % FDR 负向
col_orange = [1.0, 0.5, 0.0];    % uncorrected 正向
col_green  = [0.0, 1.0, 0.0];    % uncorrected 负向
col_white  = [1.0, 1.0, 1.0];    % 不显著

% 如果有些 ROI 没有 Harvey 坐标，则在环上给它们一个 fallback 位置
fallback_idx = find(any(isnan(node_xy), 2));
if ~isempty(fallback_idx)
    n_fb = numel(fallback_idx);
    theta = linspace(0, 2*pi, n_fb+1);
    theta = theta(1:n_fb);
    R = min(dim(2:3)) * 0.3;   % 环的半径（相对尺寸）

    cx = dim(2) * 0.5;
    cy = dim(3) * 0.5;

    for k = 1:n_fb
        bi = fallback_idx(k);
        node_xy(bi, 1) = cx + R * cos(theta(k));   % j
        node_xy(bi, 2) = cy + R * sin(theta(k));   % k
    end
end

% 新建 figure
fig = figure('Name', fig_title, 'Color', [0 0 0]);
ax  = axes('Parent', fig);
hold(ax, 'on');
axis(ax, 'image');
set(ax, 'YDir', 'normal');   % 避免上下颠倒
set(ax, 'Color', [0 0 0]);

% 显示模板切片
% 注意：imagesc(slice') 后，x 轴对应原来的 j，y 轴对应原来的 k
lo = prctile(slice(:), 5);
hi = prctile(slice(:), 95);
imagesc(slice', 'Parent', ax, [lo, hi]);
colormap(ax, gray);

xlim(ax, [1 dim(2)]);
ylim(ax, [1 dim(3)]);

% 画边（只画一次上三角即可）
for i = 1:nBase-1
    for j = i+1:nBase
        if isnan(beta_sym(i, j))
            continue;
        end
        b_ij = beta_sym(i, j);
        p_ij = p_sym(i, j);

        % 判定显著性
        is_fdr = sig_fdr(i, j);
        is_unc = sig_unc(i, j);

        % 默认颜色：白色（不显著）
        edge_color = col_white;

        if is_fdr
            if b_ij > 0
                edge_color = col_red;
            elseif b_ij < 0
                edge_color = col_blue;
            end
        elseif is_unc
            if b_ij > 0
                edge_color = col_orange;
            elseif b_ij < 0
                edge_color = col_green;
            end
        else
            % 不显著就保持白色细线
            edge_color = col_white;
        end

        % 线宽随 |β| 变化
        w = 1.0 + 3.0 * (abs(b_ij) / max_abs_beta);

        % 节点坐标 (j,k) → 显示中的 (x=j, y=k)
        x1 = node_xy(i, 1);
        y1 = node_xy(i, 2);
        x2 = node_xy(j, 1);
        y2 = node_xy(j, 2);

        if any(isnan([x1,y1,x2,y2]))
            continue;
        end

        line(ax, [x1 x2], [y1 y2], 'Color', edge_color, 'LineWidth', w);

        % 在边中点标 β 值
        xm = (x1 + x2) / 2;
        ym = (y1 + y2) / 2;
        text(ax, xm, ym, sprintf('%.3f', b_ij), ...
            'Color', edge_color, ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end
end

% 画节点
node_size = 120;
for b = 1:nBase
    x = node_xy(b, 1);
    y = node_xy(b, 2);
    if any(isnan([x,y]))
        continue;
    end
    scatter(ax, x, y, node_size, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 1.5);
    text(ax, x, y - 4, baseROI_list{b}, ...
        'Color', [1 1 1], ...
        'FontSize', 11, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');
end

title(ax, fig_title, ...
    'Color', [1 1 1], ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

end

%% ========== 子函数 3: 本地 FDR (Benjamini–Hochberg) ==========
function [h, p_crit] = fdr_bh_local(p_vals, q)
% p_vals: 向量形式的 p 值（可含 NaN）
% q     : FDR 水平，例如 0.05
%
% 返回：
%   h      : 与 p_vals 同维度的逻辑向量，true 表示通过 FDR
%   p_crit : 阈值 p*，即最大的 p，使得 p(k) <= (k/m)*q

p_flat = p_vals(:);
h      = false(size(p_flat));

valid  = ~isnan(p_flat);
pv     = p_flat(valid);

if isempty(pv)
    p_crit = NaN;
    h      = reshape(h, size(p_vals));
    return;
end

[p_sorted, sort_idx] = sort(pv);
m_valid = numel(p_sorted);

crit = (1:m_valid)' * (q / m_valid);
cmp  = (p_sorted <= crit);

if any(cmp)
    k_max  = find(cmp, 1, 'last');
    p_crit = p_sorted(k_max);
    sig_mask = (pv <= p_crit);
    h(valid) = sig_mask;
else
    p_crit = 0;
end

h = reshape(h, size(p_vals));

end
