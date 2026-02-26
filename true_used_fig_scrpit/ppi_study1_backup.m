%% ppi_group_network_BrainNet_undirected_USEFDR_whiteNonSig.m
% Group-level ROI–ROI PPI network visualization.
% 本次可视化：circle graph（节点沿圆周排布，显著边用“直线”连接）
%
% 你的要求（已实现）：
%   1) 仅可视化部分改成 circle graph（不再用弦带/BrainNet）
%   2) 只画显著边（FDR 或非校正阈值），不显著不画
%   3) 线宽 = |beta|
%   4) 线色 = 正/负（红/蓝）
%   5) 连线改为直线
%   6) 节点顺序（顺时针，从下方开始）：NTO NPO NPC1 NPC2 NPC3 NF
%   7) 连线“雾化/高级”：颜色做淡化（与白混合）+ 透明度随强度变化（EdgeAlpha）
%   8) 标题上移

clear; clc; close all;

%% ===================== USER SETTINGS =====================
subj_list  = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09','sub-10'};
condLevels = {'connect','unconnect'};
data_root  = '/home/xue/data/BIDS_prep/derivatives/ppi_roi_runwise';

DEBUG_ROILIST_PRINT_EACH_SUBJ = true;
STOP_ON_ROILIST_MISMATCH      = true;

USE_FDR   = true;
alpha_unc = 0.05;
q_fdr     = 0.05;

SHOW_TEXTBOX_INFO    = false;
MAX_EDGES_IN_TEXTBOX = 50;

% 线宽缩放（越大越粗）
EDGE_SIZE_RATIO = 0.25;

% 线色（基色）：这里已经做了“降低饱和度”的选择
% 最终还会在绘制时进一步“与白色混合”+“透明度渐变”，实现雾化高级感
POS_COLOR = [0.85 0.30 0.30];   % muted red
NEG_COLOR = [0.20 0.45 0.85];   % muted blue
NON_SIG_COLOR = [1 1 1]; 

% 雾化参数（越“高级/雾化”，就提高 whitenMix、降低 alpha 上限）
EDGE_WHITEN_MIX = 0.35;   % 0~1，与白色混合比例（建议 0.25~0.45）
EDGE_ALPHA_MIN  = 0.18;   % 弱边最小透明度
EDGE_ALPHA_MAX  = 0.75;   % 强边最大透明度

% 标题上移
TITLE_Y_NORM = 1.10;      % 标题归一化Y位置（建议 1.05~1.15）

png_dpi = 300;

% ROI MNI coords（保持原样；本脚本可视化不再用坐标，但保留不动）
roi_mni_map = struct();
roi_mni_map.NF   = [ 24  -11  52];
roi_mni_map.NPC1 = [ 22  -61  60];
roi_mni_map.NPC2 = [ 33  -40  52];
roi_mni_map.NPC3 = [ 45  -30  40];
roi_mni_map.NPO  = [ 25  -82  34];
roi_mni_map.NTO  = [ 44  -75  -4];

%% ===================== Step 1: Load subject-level FFX =====================
nSub = numel(subj_list);
beta_all = struct();
baseROI_list = [];

for ci = 1:numel(condLevels)
    beta_all.(condLevels{ci}) = [];
end

roi_list_by_subj = cell(nSub,1);

for s = 1:nSub
    subj = subj_list{s};
    ffx_file = fullfile(data_root, subj, sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral.mat', subj));
    assert(exist(ffx_file,'file')==2, 'Cannot find FFX file: %s', ffx_file);

    S = load(ffx_file, 'results_fix');
    results_fix = S.results_fix;

    base_list_this = normalize_roi_list(results_fix.baseROI_list);
    roi_list_by_subj{s} = base_list_this;

    if s == 1
        baseROI_list = base_list_this;
        nBase = numel(baseROI_list);

        for ci = 1:numel(condLevels)
            beta_all.(condLevels{ci}) = nan(nBase, nBase, nSub);
        end

        if exist('DEBUG_ROILIST_PRINT_EACH_SUBJ','var') && DEBUG_ROILIST_PRINT_EACH_SUBJ
            fprintf('\n================ [ROI LIST TEMPLATE] ================\n');
            fprintf('Template subject (subject-1) = %s\n', subj);
            fprintf('FFX file: %s\n', ffx_file);
            fprintf('nBase = %d\n', nBase);
            fprintf('baseROI_list (template): %s\n', strjoin(baseROI_list, ', '));
            for k = 1:nBase
                fprintf('  %2d: "%s"\n', k, baseROI_list{k});
            end
            fprintf('=====================================================\n\n');
        end
    else
        if ~isequal(baseROI_list, base_list_this)
            fprintf('\n================ [ROI LIST MISMATCH DETECTED] ================\n');
            fprintf('Mismatch vs subject-1 detected at subject = %s\n', subj);
            fprintf('FFX file: %s\n', ffx_file);
            report_roi_list_mismatch(baseROI_list, base_list_this, subj_list{1}, subj);
            fprintf('==============================================================\n\n');

            if exist('STOP_ON_ROILIST_MISMATCH','var') && STOP_ON_ROILIST_MISMATCH
                error('Subject %s baseROI_list mismatch vs subject-1 (see printed diagnostics above).', subj);
            else
                warning('Subject %s baseROI_list mismatch vs subject-1 (continuing).', subj);
            end
        else
            if exist('DEBUG_ROILIST_PRINT_EACH_SUBJ','var') && DEBUG_ROILIST_PRINT_EACH_SUBJ
                fprintf('[ROI LIST OK] %s matches template | nBase=%d | %s\n', ...
                    subj, numel(base_list_this), strjoin(base_list_this, ', '));
            end
        end
    end

    for ci = 1:numel(condLevels)
        condLabel = condLevels{ci};
        if isfield(results_fix, condLabel) && isfield(results_fix.(condLabel), 'beta_PPI_bilat')
            beta_all.(condLabel)(:,:,s) = results_fix.(condLabel).beta_PPI_bilat;
        else
            warning('Subject %s missing %s beta_PPI_bilat; keep NaN', subj, condLabel);
        end
    end
end

if exist('DEBUG_ROILIST_PRINT_EACH_SUBJ','var') && DEBUG_ROILIST_PRINT_EACH_SUBJ
    fprintf('\n================ [ROI LIST SUMMARY ALL SUBJECTS] ================\n');
    for s = 1:nSub
        fprintf('%s | n=%d | %s\n', subj_list{s}, numel(roi_list_by_subj{s}), strjoin(roi_list_by_subj{s}, ', '));
    end
    fprintf('==================================================================\n\n');
end

fprintf('Loaded %d subjects | nBase=%d\n', nSub, nBase);
fprintf('baseROI_list: %s\n', strjoin(baseROI_list, ', '));

%% ===================== Step 2: Group-level stats (UNDIRECTED) =====================
stats = struct();

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};
    beta_sub  = beta_all.(condLabel);

    beta_sym  = nan(nBase, nBase);
    p_sym     = nan(nBase, nBase);

    for i = 1:nBase
        for j = i+1:nBase
            vals_ij = squeeze(beta_sub(i,j,:));
            vals_ji = squeeze(beta_sub(j,i,:));

            m_sub = nan(nSub,1);
            for s = 1:nSub
                tmp = [vals_ij(s), vals_ji(s)];
                tmp = tmp(~isnan(tmp));
                if ~isempty(tmp)
                    m_sub(s) = mean(tmp);
                end
            end

            valid = ~isnan(m_sub);
            n_valid = sum(valid);

            if n_valid >= 2
                [~, p] = ttest(m_sub(valid), 0);
                bmean  = mean(m_sub(valid));
            elseif n_valid == 1
                p     = 1;
                bmean = m_sub(valid);
            else
                p     = NaN;
                bmean = NaN;
            end

            beta_sym(i,j) = bmean;
            beta_sym(j,i) = bmean;
            p_sym(i,j)    = p;
            p_sym(j,i)    = p;
        end
        beta_sym(i,i) = 0;
        p_sym(i,i)    = NaN;
    end

    mask_ut = triu(~isnan(p_sym), 1);
    p_vec   = p_sym(mask_ut);
    sig_fdr = false(size(p_sym));
    p_crit  = NaN;

    if ~isempty(p_vec)
        [h_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr(mask_ut) = h_vec;
        sig_fdr = sig_fdr | sig_fdr.';
        fprintf('[%s] FDR(q=%.3f): %d/%d edges pass (p_crit=%.3g)\n', ...
            condLabel, q_fdr, sum(h_vec), numel(h_vec), p_crit);
    else
        fprintf('[%s] No valid p values for FDR\n', condLabel);
    end

    stats.(condLabel).beta_sym = beta_sym;
    stats.(condLabel).p_sym    = p_sym;
    stats.(condLabel).sig_fdr  = sig_fdr;
    stats.(condLabel).p_crit   = p_crit;
end

%% ===================== Step 3: connect - unconnect (UNDIRECTED) =====================
if all(isfield(beta_all, {'connect','unconnect'}))
    beta_diff_sub = beta_all.connect - beta_all.unconnect;

    beta_sym  = nan(nBase, nBase);
    p_sym     = nan(nBase, nBase);

    for i = 1:nBase
        for j = i+1:nBase
            vals_ij = squeeze(beta_diff_sub(i,j,:));
            vals_ji = squeeze(beta_diff_sub(j,i,:));

            m_sub = nan(nSub,1);
            for s = 1:nSub
                tmp = [vals_ij(s), vals_ji(s)];
                tmp = tmp(~isnan(tmp));
                if ~isempty(tmp)
                    m_sub(s) = mean(tmp);
                end
            end

            valid = ~isnan(m_sub);
            n_valid = sum(valid);

            if n_valid >= 2
                [~, p] = ttest(m_sub(valid), 0);
                bmean  = mean(m_sub(valid));
            elseif n_valid == 1
                p     = 1;
                bmean = m_sub(valid);
            else
                p     = NaN;
                bmean = NaN;
            end

            beta_sym(i,j) = bmean;
            beta_sym(j,i) = bmean;
            p_sym(i,j)    = p;
            p_sym(j,i)    = p;
        end
        beta_sym(i,i) = 0;
        p_sym(i,i)    = NaN;
    end

    mask_ut = triu(~isnan(p_sym), 1);
    p_vec   = p_sym(mask_ut);
    sig_fdr = false(size(p_sym));
    p_crit  = NaN;

    if ~isempty(p_vec)
        [h_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr(mask_ut) = h_vec;
        sig_fdr = sig_fdr | sig_fdr.';
        fprintf('[diff] FDR(q=%.3f): %d/%d edges pass (p_crit=%.3g)\n', ...
            q_fdr, sum(h_vec), numel(h_vec), p_crit);
    else
        fprintf('[diff] No valid p values for FDR\n');
    end

    stats.diff.beta_sym = beta_sym;
    stats.diff.p_sym    = p_sym;
    stats.diff.sig_fdr  = sig_fdr;
    stats.diff.p_crit   = p_crit;
else
    warning('Missing connect/unconnect in beta_all; skip diff');
end

%% ===================== Step 4: Circle-graph visualization (STRAIGHT + FOGGY) =====================
vis_mode = ternary(USE_FDR, 'FDR', sprintf('UNC(p<%.3g)', alpha_unc));

out_dir = fullfile(data_root, '_GROUP_CircleGraph', sprintf('PPI_circle_sigEdges_%s', vis_mode));
out_dir = sanitize_dirname(out_dir);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fprintf('Output dir: %s\n', out_dir);

tags = {'connect','unconnect'};
if isfield(stats,'diff'), tags{end+1} = 'diff'; end

% 节点扇区颜色（柔和）
NODE_RING_COLORS = [
    0.56 0.70 0.60;  % pastel green
    0.86 0.80 0.55;  % pastel yellow
    0.60 0.70 0.84;  % pastel blue
    0.93 0.86 0.70;  % pastel beige
    0.70 0.82 0.88;  % light blue
    0.80 0.80 0.80;  % light gray
];

% 固定绘图顺序（顺时针，从下方开始）
DESIRED_ORDER = {'NTO','NPO','NPC1','NPC2','NPC3','NF'};
[tf, loc] = ismember(DESIRED_ORDER, baseROI_list);
assert(all(tf), 'DESIRED_ORDER contains ROI not in baseROI_list. Missing: %s', strjoin(DESIRED_ORDER(~tf), ', '));

new2old = loc(:)'; % newIndex -> oldIndex
old2new = nan(1, nBase);
for k = 1:nBase
    old2new(new2old(k)) = k;
end

labels_plot = DESIRED_ORDER;

for ti = 1:numel(tags)
    tag = tags{ti};
    assert(isfield(stats, tag), 'stats.%s not found', tag);

    st = stats.(tag);
    beta_sym = st.beta_sym;
    p_sym    = st.p_sym;
    sig_fdr  = st.sig_fdr;

    if all(isnan(beta_sym(:)))
        warning('[%s] beta_sym is all-NaN; skip plot', tag);
        continue;
    end

    if USE_FDR
        mask_sig = sig_fdr;
    else
        mask_sig = isfinite(p_sym) & (p_sym < alpha_unc);
    end

    % 收集显著边（原始索引）
    edgeList_old = struct('u',{},'v',{},'beta',{},'p',{});
    kout = 0;
    for u = 1:nBase-1
        for v = u+1:nBase
            b = beta_sym(u,v);
            if isnan(b), continue; end
            if mask_sig(u,v)
                kout = kout + 1;
                edgeList_old(kout).u = u;
                edgeList_old(kout).v = v;
                edgeList_old(kout).beta = b;
                edgeList_old(kout).p = p_sym(u,v);
            end
        end
    end

    % 映射到绘图顺序
    edgeList = struct('u',{},'v',{},'beta',{},'p',{});
    for k = 1:numel(edgeList_old)
        u_old = edgeList_old(k).u;
        v_old = edgeList_old(k).v;
        u_new = old2new(u_old);
        v_new = old2new(v_old);
        if any(isnan([u_new v_new]))
            error('ROI index mapping failed. Check DESIRED_ORDER vs baseROI_list.');
        end
        uu = min(u_new, v_new);
        vv = max(u_new, v_new);

        edgeList(k).u = uu; %#ok<AGROW>
        edgeList(k).v = vv;
        edgeList(k).beta = edgeList_old(k).beta;
        edgeList(k).p    = edgeList_old(k).p;
    end

    % ---- draw ----
    fig = figure('Color','w', 'Renderer','opengl'); % 透明度需要 opengl
    set(fig, 'Position', [100 100 900 900]);
    ax = axes(fig);
    hold(ax, 'on');

    plot_circle_sigEdges(ax, labels_plot, edgeList, ...
        'posColor', POS_COLOR, 'negColor', NEG_COLOR, ...
        'edgeSizeRatio', EDGE_SIZE_RATIO, ...
        'nodeRingColors', NODE_RING_COLORS, ...
        'startFromBottom', true, ...
        'edgeWhitenMix', EDGE_WHITEN_MIX, ...
        'edgeAlphaMin', EDGE_ALPHA_MIN, ...
        'edgeAlphaMax', EDGE_ALPHA_MAX);

    if SHOW_TEXTBOX_INFO
        add_textbox_edges(fig, tag, vis_mode, edgeList, labels_plot, MAX_EDGES_IN_TEXTBOX);
    end

    % 标题上移
    t = title(ax, sprintf('Group PPI | %s | %s', tag, vis_mode), 'Interpreter','none');
    t.Units = 'normalized';
    t.Position(2) = TITLE_Y_NORM;

    % 给标题留空间，避免导出裁剪
    ax.Position(2) = ax.Position(2) - 0.02;
    ax.Position(4) = ax.Position(4) - 0.02;

    % ---- save ----
    out_png = fullfile(out_dir, sprintf('CircleGraph_PPI_%s_%s.png', vis_mode, tag));
    out_pdf = fullfile(out_dir, sprintf('CircleGraph_PPI_%s_%s.pdf', vis_mode, tag));

    % PNG：保留透明度效果（opengl）
    print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi), '-opengl');

    % PDF：若用 painters 可能丢透明度；这里用 opengl 保证风格一致（但可能是栅格化）
    print(fig, out_pdf, '-dpdf', '-opengl');

    close(fig);
    fprintf('[%s] saved: %s\n', tag, out_png);
end

fprintf('[DONE] CircleGraph outputs: %s\n', out_dir);

%% ===================== Local FDR (BH) =====================
function [h, p_crit] = fdr_bh_local(p_vals, q)
p = p_vals(:);
h = false(size(p));
valid = ~isnan(p);
pv = p(valid);
if isempty(pv)
    p_crit = NaN;
    h = reshape(h, size(p_vals));
    return;
end
pv_sorted = sort(pv);
m = numel(pv_sorted);
crit = (1:m)' * (q / m);
pass = pv_sorted <= crit;
if any(pass)
    kmax = find(pass, 1, 'last');
    p_crit = pv_sorted(kmax);
    h(valid) = (pv <= p_crit);
else
    p_crit = 0;
end
h = reshape(h, size(p_vals));
end

%% ===================== Circle-graph plotter (STRAIGHT + FOGGY) =====================
function plot_circle_sigEdges(ax, roi_labels, edgeList, varargin)
% circle graph（节点圆周 + 直线连接）
% 雾化高级感：
%   - 颜色与白色混合（edgeWhitenMix）
%   - 透明度随强度变化（edgeAlphaMin~edgeAlphaMax）
% 说明：透明度依赖 OpenGL renderer

p = inputParser;
p.addParameter('posColor', [0.85 0.30 0.30]);
p.addParameter('negColor', [0.20 0.45 0.85]);
p.addParameter('edgeSizeRatio', 0.25);
p.addParameter('nodeRingColors', []);
p.addParameter('startFromBottom', true);
p.addParameter('edgeWhitenMix', 0.35);
p.addParameter('edgeAlphaMin', 0.18);
p.addParameter('edgeAlphaMax', 0.75);
p.parse(varargin{:});
opt = p.Results;

n = numel(roi_labels);

% 节点角度（顺时针）
theta0 = ternary(opt.startFromBottom, -pi/2, pi/2);
theta = linspace(theta0, theta0-2*pi, n+1);
theta(end) = [];

R_node_outer  = 1.18;
R_node_inner  = 1.02;
R_edge_anchor = 1.00;
labelR        = 1.28;

% 节点颜色
if isempty(opt.nodeRingColors)
    cc = lines(n);
else
    cc = opt.nodeRingColors;
    if size(cc,1) < n
        cc = repmat(cc, ceil(n/size(cc,1)), 1);
    end
    cc = cc(1:n,:);
end

% 画扇区
dtheta = 2*pi/n;
for i = 1:n
    t1 = theta(i) - dtheta/2;
    t2 = theta(i) + dtheta/2;
    draw_ring_sector(ax, t1, t2, R_node_inner, R_node_outer, cc(i,:));
end

% 标签
for i = 1:n
    ang = theta(i);
    x = labelR * cos(ang);
    y = labelR * sin(ang);

    rot = rad2deg(ang);
    ha = 'left';
    if rot < -90 || rot > 90
        rot = rot + 180;
        ha = 'right';
    end

    text(ax, x, y, roi_labels{i}, ...
        'Rotation', rot, 'HorizontalAlignment', ha, ...
        'VerticalAlignment', 'middle', 'FontSize', 12, ...
        'Interpreter', 'none');
end

% 线宽/透明度归一化
if isempty(edgeList)
    maxAbs = 1;
else
    maxAbs = max(abs([edgeList.beta]));
    if maxAbs <= 0 || isnan(maxAbs), maxAbs = 1; end
end

LW_MIN   = 0.8;
LW_SCALE = 8.0;

% 边颜色做淡化（与白混合）
posCol_soft = mix_with_white(opt.posColor, opt.edgeWhitenMix);
negCol_soft = mix_with_white(opt.negColor, opt.edgeWhitenMix);

% 画边（直线 + 透明度）
for k = 1:numel(edgeList)
    u = edgeList(k).u;
    v = edgeList(k).v;
    b = edgeList(k).beta;

    p0 = R_edge_anchor * [cos(theta(u)), sin(theta(u))];
    p3 = R_edge_anchor * [cos(theta(v)), sin(theta(v))];

    strength = abs(b) / maxAbs;  % 0~1
    lw = LW_MIN + strength * (LW_SCALE * opt.edgeSizeRatio);

    % 透明度：弱边更透明、强边更实
    a = opt.edgeAlphaMin + strength * (opt.edgeAlphaMax - opt.edgeAlphaMin);

    if b >= 0
        col = posCol_soft;
    else
        col = negCol_soft;
    end

    % 用 patch 画“只有边的线段”，这样可以使用 EdgeAlpha
    patch(ax, 'XData', [p0(1) p3(1)], 'YData', [p0(2) p3(2)], ...
        'FaceColor', 'none', 'EdgeColor', col, ...
        'LineWidth', lw, 'EdgeAlpha', a);
end

% 图例（图例不显示 alpha 属于 MATLAB 限制；含义仍准确）
h1 = plot(ax, nan, nan, '-', 'Color', posCol_soft, 'LineWidth', 2);
h2 = plot(ax, nan, nan, '-', 'Color', negCol_soft, 'LineWidth', 2);
legend(ax, [h1 h2], {'beta > 0', 'beta < 0'}, ...
    'Location', 'southoutside', 'Orientation','horizontal', 'Box','off');

axis(ax, 'equal'); axis(ax, 'off');
xlim(ax, [-1.45 1.45]); ylim(ax, [-1.45 1.45]);
end

function c = mix_with_white(c0, w)
% w=0: 原色；w=1: 白色
w = max(0, min(1, w));
c = (1-w)*c0 + w*[1 1 1];
end

function draw_ring_sector(ax, t1, t2, r_in, r_out, faceColor)
tt = linspace(t1, t2, 80);
xo = r_out*cos(tt); yo = r_out*sin(tt);
xi = r_in*cos(fliplr(tt)); yi = r_in*sin(fliplr(tt));
x = [xo xi];
y = [yo yi];
patch(ax, x, y, faceColor, 'EdgeColor','w', 'LineWidth', 1);
end

%% ===================== Helpers =====================
function add_textbox_edges(fig, tag, vis_mode, edgeList, roi_labels, max_lines)
lines = {};
lines{end+1} = sprintf('PPI Group | %s', tag);
lines{end+1} = sprintf('Sig: %s', vis_mode);
lines{end+1} = sprintf('Edges(sig): %d', numel(edgeList));
lines{end+1} = '---';

if ~isempty(edgeList)
    betas = [edgeList.beta];
    [~, ord] = sort(abs(betas), 'descend');
    edgeList = edgeList(ord);
end

n_show = min(numel(edgeList), max_lines);
for k = 1:n_show
    u = edgeList(k).u; v = edgeList(k).v;
    b = edgeList(k).beta;
    p = edgeList(k).p;
    lines{end+1} = sprintf('%s--%s  beta=%+.3f  p=%.3g', roi_labels{u}, roi_labels{v}, b, p);
end
if numel(edgeList) > n_show
    lines{end+1} = sprintf('... (%d more)', numel(edgeList) - n_show);
end

annotation(fig, 'textbox', [0.01 0.01 0.40 0.35], ...
    'String', lines, 'FitBoxToText','on', ...
    'BackgroundColor','w', 'EdgeColor',[0.2 0.2 0.2], ...
    'Interpreter','none', 'FontSize', 9);
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function out = sanitize_dirname(in)
out = strrep(in, ':', '_');
out = strrep(out, '(', '');
out = strrep(out, ')', '');
out = strrep(out, '<', '');
out = strrep(out, '>', '');
out = strrep(out, '=', '');
out = strrep(out, ' ', '_');
end

function L = normalize_roi_list(in)
if isempty(in)
    L = {};
    return;
end
if isstring(in)
    in = cellstr(in);
elseif ischar(in)
    in = cellstr(in);
elseif iscell(in)
else
    try
        in = cellstr(in);
    catch
        error('normalize_roi_list: unsupported type: %s', class(in));
    end
end
in = in(:);
L = cell(size(in));
for i = 1:numel(in)
    if isempty(in{i})
        L{i} = '';
    else
        L{i} = strtrim(in{i});
    end
end
end

function report_roi_list_mismatch(refList, curList, subj_ref, subj_cur)
n1 = numel(refList);
n2 = numel(curList);

fprintf('Template (subject-1) = %s | nBase=%d\n', subj_ref, n1);
fprintf('Current  (subject)   = %s | nBase=%d\n', subj_cur, n2);

fprintf('\n-- Template list (with indices) --\n');
for i = 1:n1
    fprintf('  %2d: "%s"\n', i, refList{i});
end

fprintf('\n-- Current list (with indices) --\n');
for i = 1:n2
    fprintf('  %2d: "%s"\n', i, curList{i});
end

m = min(n1, n2);
firstDiff = [];
if m > 0
    eqpos = strcmp(refList(1:m), curList(1:m));
    firstDiff = find(~eqpos, 1, 'first');
end

fprintf('\n-- Position-wise comparison --\n');
if isempty(firstDiff) && n1 == n2
    fprintf('Same length, but still mismatch: likely hidden differences.\n');
elseif isempty(firstDiff) && n1 ~= n2
    fprintf('Prefix identical for first %d items, but lengths differ.\n', m);
else
    fprintf('First mismatch index = %d\n', firstDiff);
    if ~isempty(firstDiff)
        fprintf('  template[%d] = "%s"\n', firstDiff, refList{firstDiff});
        fprintf('  current [%d] = "%s"\n', firstDiff, curList{firstDiff});
    end
end

missing_in_cur = setdiff(refList, curList, 'stable');
extra_in_cur   = setdiff(curList, refList, 'stable');

fprintf('\n-- Set difference (ignore order) --\n');
fprintf('Missing in current: %s\n', strjoin(missing_in_cur, ', '));
fprintf('Extra in current  : %s\n', strjoin(extra_in_cur, ', '));

sameSet = isempty(missing_in_cur) && isempty(extra_in_cur);
if sameSet
    fprintf('\n-- Same set but different order: suggested reindex --\n');
    [tf2, loc2] = ismember(refList, curList);
    if all(tf2)
        fprintf('loc = %s\n', mat2str(loc2(:)'));
        fprintf('M_reordered = M(loc, loc)\n');
    end
else
    fprintf('\n-- Not the same set: not only ordering issue. --\n');
end
fprintf('\n');
end