function ppi_group_brainnet_signColor_noArrow()
% =========================================================
% Visualization ONLY for group ROI–ROI PPI results (CONTRASTS ONLY).
% - Read saved group MAT (no recomputation)
% - Plot ONLY contrasts from pair_list (supports 3+ conditions, multiple contrasts)
% - Undirected display (merge directions), keep sign color:
%       beta > 0 : red
%       beta < 0 : blue
% - Skip tags with ZERO significant edges
% - Optional textbox listing significant edges
%
% NOTE:
%   Circle-graph style uses transparency (alpha). PDF may be rasterized.
%   If you need vector PDF: set EDGE_ALPHA_MIN=1; EDGE_ALPHA_MAX=1;
%   and use print(...,'-dpdf','-painters') (see below).
% =========================================================

clear; clc; close all;

%% ===================== USER SETTINGS =====================
data_root = '/home/xue/data/prf_result/results_ppi_roi_runwise_3level_scanMap';

% --- Significance control ---
USE_FDR   = false;    % true: use st.sig_fdr; false: use uncorrected p<thr
alpha_unc = 0.05;    % used only when USE_FDR=false

% --- Textbox control ---
SHOW_TEXTBOX_INFO = false
MAX_EDGES_IN_TEXTBOX = 50;

% --- Edge thickness scale (smaller -> thinner) ---
EDGE_SIZE_RATIO = 0.25;   % try 0.15 ~ 0.60

% --- Sign colors (circle scheme) ---
POS_COLOR = [0.85 0.30 0.30];   % muted red
NEG_COLOR = [0.20 0.45 0.85];   % muted blue

% --- Circle style controls (fog/soft look) ---
EDGE_WHITEN_MIX = 0.35;   % mix with white (0~1)
EDGE_ALPHA_MIN  = 0.18;   % weak edge alpha
EDGE_ALPHA_MAX  = 0.75;   % strong edge alpha

% --- Node ring (pastel blocks) ---
NODE_RING_COLORS = [
    0.56 0.70 0.60;
    0.86 0.80 0.55;
    0.60 0.70 0.84;
    0.93 0.86 0.70;
    0.70 0.82 0.88;
    0.80 0.80 0.80;
];

% --- Optional: enforce a plotting order (must match ROIs present) ---
% If you want to force a specific circular order, edit here.
% If not all names match baseROI_list, it falls back to baseROI_list order.
DESIRED_ORDER = {'NTO','NPO','NPC1','NPC2','NPC3','NF'};

png_dpi = 300;

%% ===================== Locate input MAT (auto) =====================
mat_u = fullfile(data_root, '_GROUP_3level_undirected', 'ppi_group_stats_3level_scanMap_undirected.mat');
mat_d = fullfile(data_root, '_GROUP_3level_directed',   'ppi_group_stats_3level_scanMap_directed.mat');

if exist(mat_u,'file')==2
    in_mat = mat_u;
elseif exist(mat_d,'file')==2
    in_mat = mat_d;
else
    error('Cannot find group MAT. Tried:\n  %s\n  %s', mat_u, mat_d);
end

S = load(in_mat, 'stats','baseROI_list','condLevels','pair_list','alpha_unc','q_fdr','edge_mode');
stats        = S.stats;
baseROI_list = S.baseROI_list;
condLevels   = S.condLevels; %#ok<NASGU>  % keep loaded (not used; contrasts only)
pair_list    = S.pair_list;

vis_mode = ternary(USE_FDR, 'FDR', sprintf('UNC(p<%.3g)', alpha_unc));

fprintf('Loaded MAT: %s\n', in_mat);
fprintf('MAT edge_mode=%s | vis_mode=%s | EDGE_SIZE_RATIO=%.3f | SHOW_TEXTBOX_INFO=%d\n', ...
    S.edge_mode, vis_mode, EDGE_SIZE_RATIO, SHOW_TEXTBOX_INFO);

%% ===================== Output directory =====================
[out_root, ~, ~] = fileparts(in_mat);
out_dir = fullfile(out_root, sprintf('Circle_PPI_signColor_CONTRASTONLY_%s', vis_mode));
out_dir = sanitize_dirname(out_dir);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
fprintf('Output dir: %s\n', out_dir);

%% ===================== Build tag list (CONTRAST ONLY; supports 3+ conditions) =====================
% 原脚本：tags = [condLevels, pair_list(:,3)]
% 你的需求：只画 pair_list 中所有 contrasts（多个）
assert(iscell(pair_list) && ~isempty(pair_list), 'pair_list is empty: MAT does not contain contrasts to plot.');

if iscell(pair_list{1})
    % pair_list = { {condA, condB, tag}, ... }
    tags = cellfun(@(x) x{3}, pair_list, 'UniformOutput', false);
else
    % pair_list = Nx3 cell
    tags = pair_list(:,3)';
end

fprintf('Contrast tags to plot (from pair_list):\n');
fprintf('  %s\n', tags{:});

%% ===================== Plot-order remapping =====================
nBase = numel(baseROI_list);

use_desired = (numel(DESIRED_ORDER)==nBase) && all(ismember(DESIRED_ORDER, baseROI_list));
if use_desired
    labels_plot = DESIRED_ORDER;
    [~, loc] = ismember(DESIRED_ORDER, baseROI_list);
    old2new = nan(1, nBase);
    for k = 1:nBase
        old2new(loc(k)) = k; % old index -> new index
    end
else
    labels_plot = baseROI_list;
    old2new = 1:nBase;
end

%% ===================== Plot loop =====================
for ti = 1:numel(tags)
    tag = tags{ti};
    assert(isfield(stats, tag), 'stats.%s not found in MAT', tag);
    st = stats.(tag);

    % --- significance mask from MAT ---
    if USE_FDR
        mask = st.sig_fdr;
    else
        mask = isfinite(st.p_mat) & (st.p_mat < alpha_unc);
    end

    n_sig_entries = nnz(mask & ~eye(nBase));
    fprintf('[%s] sig entries (mask) = %d\n', tag, n_sig_entries);

    if n_sig_entries == 0
        fprintf('  -> skip (no significant edges)\n');
        continue;
    end

    % -------- Convert to undirected edges (no arrows) --------
    W = zeros(nBase);
    edgeList = struct('u',{},'v',{},'beta',{},'p',{},'src',{},'dst',{});
    kout = 0;

    for u = 1:nBase-1
        for v = u+1:nBase
            b_uv = st.beta_mat(u,v);  p_uv = st.p_mat(u,v);  m_uv = mask(u,v);
            b_vu = st.beta_mat(v,u);  p_vu = st.p_mat(v,u);  m_vu = mask(v,u);

            if ~(m_uv || m_vu)
                continue;
            end

            cand_beta = [];
            cand_p    = [];
            cand_src  = [];
            cand_dst  = [];

            if m_uv && isfinite(b_uv)
                cand_beta(end+1,1) = b_uv; %#ok<AGROW>
                cand_p(end+1,1)    = p_uv; %#ok<AGROW>
                cand_src(end+1,1)  = v;    %#ok<AGROW>
                cand_dst(end+1,1)  = u;    %#ok<AGROW>
            end
            if m_vu && isfinite(b_vu)
                cand_beta(end+1,1) = b_vu; %#ok<AGROW>
                cand_p(end+1,1)    = p_vu; %#ok<AGROW>
                cand_src(end+1,1)  = u;    %#ok<AGROW>
                cand_dst(end+1,1)  = v;    %#ok<AGROW>
            end

            if isempty(cand_beta)
                continue;
            end

            [~, idx] = max(abs(cand_beta));
            beta = cand_beta(idx);
            pval = cand_p(idx);
            src  = cand_src(idx);
            dst  = cand_dst(idx);

            W(u,v) = beta; W(v,u) = beta;

            kout = kout + 1;
            edgeList(kout).u    = u;
            edgeList(kout).v    = v;
            edgeList(kout).beta = beta;
            edgeList(kout).p    = pval;
            edgeList(kout).src  = src;
            edgeList(kout).dst  = dst;
        end
    end

    n_edges = nnz(triu(W,1));
    fprintf('  undirected edges kept = %d | beta range=[%.4g, %.4g]\n', ...
        n_edges, min(W(:)), max(W(:)));

    if n_edges == 0
        fprintf('  -> skip (empty after merge)\n');
        continue;
    end

    % -------- Remap edges to plotting order --------
    edgeList_plot = struct('u',{},'v',{},'beta',{},'p',{});
    for k = 1:numel(edgeList)
        u_new = old2new(edgeList(k).u);
        v_new = old2new(edgeList(k).v);
        uu = min(u_new, v_new);
        vv = max(u_new, v_new);

        edgeList_plot(k).u    = uu; %#ok<AGROW>
        edgeList_plot(k).v    = vv;
        edgeList_plot(k).beta = edgeList(k).beta;
        edgeList_plot(k).p    = edgeList(k).p;
    end

    % -------- Circle visualization --------
    fig = figure('Color','w', 'Renderer','opengl');
    set(fig, 'Position', [100 100 900 900]);
    ax = axes(fig); hold(ax, 'on');

    plot_circle_sigEdges(ax, labels_plot, edgeList_plot, ...
        'posColor', POS_COLOR, 'negColor', NEG_COLOR, ...
        'edgeSizeRatio', EDGE_SIZE_RATIO, ...
        'nodeRingColors', NODE_RING_COLORS, ...
        'startFromBottom', true, ...
        'edgeWhitenMix', EDGE_WHITEN_MIX, ...
        'edgeAlphaMin', EDGE_ALPHA_MIN, ...
        'edgeAlphaMax', EDGE_ALPHA_MAX);

    % 修改标题字体：Times New Roman 粗体，字号18
    title(ax, sprintf('Group PPI | %s | %s', tag, vis_mode), ...
        'Interpreter','none', 'FontName','Times New Roman', ...
        'FontWeight','bold', 'FontSize', 18);

    if SHOW_TEXTBOX_INFO
        add_textbox_edges(fig, tag, vis_mode, edgeList_plot, labels_plot, MAX_EDGES_IN_TEXTBOX);
    end

    out_png = fullfile(out_dir, sprintf('Circle_%s_%s.png', vis_mode, tag));
    out_pdf = fullfile(out_dir, sprintf('Circle_%s_%s.pdf', vis_mode, tag));

    print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi), '-opengl');

    % 默认：保留透明效果（可能栅格化）
    print(fig, out_pdf, '-dpdf', '-opengl');

    % 如果必须“纯矢量PDF”：把 EDGE_ALPHA_MIN/MAX 设为 1，并改用下面这一行：
    % print(fig, out_pdf, '-dpdf', '-painters');

    close(fig);

    fprintf('  saved: %s\n', out_png);
end

fprintf('[DONE] outputs: %s\n', out_dir);

end

%% ===================== Helpers =====================

function plot_circle_sigEdges(ax, roi_labels, edgeList, varargin)
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

theta0 = ternary(opt.startFromBottom, -pi/2, pi/2);
theta = linspace(theta0, theta0-2*pi, n+1);
theta(end) = [];

R_node_outer  = 1.18;
R_node_inner  = 1.02;
R_edge_anchor = 1.00;
% 将标签半径增大，使其远离节点环（原1.28，现改为1.42）
labelR        = 1.42;

if isempty(opt.nodeRingColors)
    cc = lines(n);
else
    cc = opt.nodeRingColors;
    if size(cc,1) < n
        cc = repmat(cc, ceil(n/size(cc,1)), 1);
    end
    cc = cc(1:n,:);
end

dtheta = 2*pi/n;
for i = 1:n
    t1 = theta(i) - dtheta/2;
    t2 = theta(i) + dtheta/2;
    draw_ring_sector(ax, t1, t2, R_node_inner, R_node_outer, cc(i,:));
end

% 绘制 ROI 标签（水平、Times New Roman 粗体、字号14）
for i = 1:n
    ang = theta(i);
    x = labelR * cos(ang);
    y = labelR * sin(ang);

    % 水平放置，根据节点所在半圆选择对齐方式，使文本向外辐射
    if cos(ang) > 0   % 右半圆，文本向左展开
        ha = 'right';
    else               % 左半圆，文本向右展开
        ha = 'left';
    end

    text(ax, x, y, roi_labels{i}, ...
        'Rotation', 0, ...                         % 保持水平
        'HorizontalAlignment', ha, ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', ...
        'FontSize', 14, ...                         % 比原来12大一些
        'Interpreter', 'none');
end

if isempty(edgeList)
    maxAbs = 1;
else
    maxAbs = max(abs([edgeList.beta]));
    if maxAbs <= 0 || isnan(maxAbs), maxAbs = 1; end
end

LW_MIN   = 0.8;
LW_SCALE = 8.0;

posCol_soft = mix_with_white(opt.posColor, opt.edgeWhitenMix);
negCol_soft = mix_with_white(opt.negColor, opt.edgeWhitenMix);

for k = 1:numel(edgeList)
    u = edgeList(k).u;
    v = edgeList(k).v;
    b = edgeList(k).beta;

    p0 = R_edge_anchor * [cos(theta(u)), sin(theta(u))];
    p3 = R_edge_anchor * [cos(theta(v)), sin(theta(v))];

    strength = abs(b) / maxAbs;
    lw = LW_MIN + strength * (LW_SCALE * opt.edgeSizeRatio);
    a  = opt.edgeAlphaMin + strength * (opt.edgeAlphaMax - opt.edgeAlphaMin);

    col = posCol_soft;
    if b < 0, col = negCol_soft; end

    patch(ax, 'XData', [p0(1) p3(1)], 'YData', [p0(2) p3(2)], ...
        'FaceColor', 'none', 'EdgeColor', col, ...
        'LineWidth', lw, 'EdgeAlpha', a);
end

h1 = plot(ax, nan, nan, '-', 'Color', posCol_soft, 'LineWidth', 2);
h2 = plot(ax, nan, nan, '-', 'Color', negCol_soft, 'LineWidth', 2);
leg = legend(ax, [h1 h2], {'beta > 0', 'beta < 0'}, ...
    'Location', 'southoutside', 'Orientation','horizontal', 'Box','off');
% 设置图例字体：Times New Roman，正常粗细，字号12
set(leg, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 12);

axis(ax, 'equal'); axis(ax, 'off');
% 调整坐标轴范围以容纳更大的标签（原[-1.45 1.45]，现改为[-1.6 1.6]）
xlim(ax, [-1.6 1.6]); ylim(ax, [-1.6 1.6]);
end

function c = mix_with_white(c0, w)
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

function add_textbox_edges(fig, tag, vis_mode, edgeList, roi_labels, max_lines)
lines = {};
lines{end+1} = sprintf('PPI Group | %s', tag);
lines{end+1} = sprintf('Sig: %s', vis_mode);
lines{end+1} = sprintf('Edges: %d', numel(edgeList));
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

