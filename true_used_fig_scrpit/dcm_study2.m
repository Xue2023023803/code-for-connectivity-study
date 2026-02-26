function lastest_plot_PEB3_polyLQ_BrainNet_PP095_withValueAndDirection()
% ==============================================================
% BrainNet visualization for polyLQ PEB3 (mean/Linear/Quadratic)
% Controlled by ONE parameter: MATRIX_MODE = 'A' | 'B' | 'AB'
%
% [VIS MOD]
%   - muted red/blue edges + node dark gray + opacity
%   - optional uniform edge width
%   - remove BrainNet colorbar/legend
%   - ADD heatmap (square PNG + vector PDF) with bottom-left NTO, fixed layout
%   - [NEW] 手动Y轴文字标签：拉开与轴的距离
% ==============================================================

clear; clc; close all;

%% ===================== USER SETTINGS =====================
root_prf = '/home/xue/data/prf_result';
prf_ses  = 'ses-01';

MATRIX_MODE = 'AB';     % 'A' / 'B' / 'AB'
pp_thr      = 0.95;
effectsToPlot = {'mean','Linear','Quadratic'};

SHOW_TEXTBOX_INFO = false;
png_dpi = 300;

% Heatmap debug：打印 “target <- source : Ep”
DEBUG_HEATMAP_MAPPING = false;

%% --- 可视化方案(来自参考脚本) ---
VIS_SCHEME = 1;  % 1=directed arrow, 2=undirected（这里保持箭头方案）
EDGE_COLOR_MODE   = 'threshold';
EDGE_COLOR_THR    = 0;
EDGE_POS_COLOR    = [0.88 0.39 0.39];   % muted red
EDGE_NEG_COLOR    = [0.35 0.55 0.78];   % muted blue
EDGE_SIZE_ABS     = true;
EDGE_COLOR_ABS    = false;
EDGE_SIZE_SCALE   = 0.35;

UNIFORM_EDGE_WIDTH = true;
UNIFORM_EDGE_VALUE = 1;

EDGE_OPACITY_MODE = 1;
EDGE_OPACITY_SAME = 0.45;

NODE_RGB  = [0.18 0.20 0.24];
NODE_SIZE = 5;

DRAW_DIRECTED = (VIS_SCHEME == 1);

%% ===================== ROI order & coords =====================
roi_labels = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};
roi_mni = [ ...
     24  -11  52;
     22  -61  60;
     33  -40  52;
     45  -30  40;
     25  -82  34;
     44  -75  -4];
nROI = numel(roi_labels);

%% ===================== BrainNet detect =====================
bn_mapcfg = which('BrainNet_MapCfg');
assert(~isempty(bn_mapcfg), 'BrainNet_MapCfg not found. Please add BrainNetViewer to MATLAB path.');
BNV_DIR = fileparts(bn_mapcfg);
while ~exist(fullfile(BNV_DIR,'Data'),'dir')
    parent = fileparts(BNV_DIR);
    if strcmp(parent, BNV_DIR), error('Cannot locate BrainNet "Data" folder.'); end
    BNV_DIR = parent;
end
addpath(genpath(BNV_DIR));
fprintf('BrainNet root detected: %s\n', BNV_DIR);

meshFile = fullfile(BNV_DIR,'Data','SurfTemplate','BrainMesh_ICBM152Right_smoothed.nv');
if exist(meshFile,'file') ~= 2
    meshFile = fullfile(BNV_DIR,'Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv');
end
fprintf('Using mesh: %s\n', meshFile);

%% ===================== Resolve PEB mat =====================
MATRIX_MODE = upper(string(MATRIX_MODE));
assert(ismember(MATRIX_MODE, ["A","B","AB"]), 'MATRIX_MODE must be A/B/AB.');

peb_dir  = fullfile(root_prf, sprintf('PEB_group_AL3cond_poly%s', char(MATRIX_MODE)), prf_ses);
peb_mat1 = fullfile(peb_dir, sprintf('PEB3_polyLQ_%s.mat', char(MATRIX_MODE)));
peb_mat2 = fullfile(peb_dir, 'PEB3_polyLQ_core.mat');

if exist(peb_mat1,'file')==2
    peb_mat = peb_mat1;
elseif exist(peb_mat2,'file')==2
    peb_mat = peb_mat2;
else
    error('PEB mat not found. Tried:\n  %s\n  %s', peb_mat1, peb_mat2);
end
fprintf('Using PEB mat: %s\n', peb_mat);

out_dir = fullfile(peb_dir, sprintf('plot_BrainNet_polyLQ_PP%03d', round(pp_thr*100)));
if ~exist(out_dir,'dir'), mkdir(out_dir); end

node_file = fullfile(out_dir, 'nodes_6ROIs.node');
write_bnv_node(node_file, roi_mni, roi_labels, NODE_SIZE);

cfg_file = fullfile(out_dir, sprintf('cfg_edgeOpacity%.2f_nodeRGB_%s_scheme%d.mat', ...
    EDGE_OPACITY_SAME, rgb_tag(NODE_RGB), VIS_SCHEME));

make_bnv_cfg_edgeSign(cfg_file, EDGE_COLOR_MODE, EDGE_COLOR_THR, EDGE_POS_COLOR, EDGE_NEG_COLOR, ...
                      EDGE_SIZE_ABS, EDGE_COLOR_ABS, DRAW_DIRECTED, EDGE_SIZE_SCALE, ...
                      EDGE_OPACITY_MODE, EDGE_OPACITY_SAME, NODE_RGB);

fprintf('Using cfg: %s\n', cfg_file);

%% ===================== Load BMA3 =====================
BMA = load_BMA3_from_PEBmat(peb_mat);

Pnames = BMA.Pnames;
if isstring(Pnames), Pnames = cellstr(Pnames); end
if ischar(Pnames),   Pnames = cellstr(Pnames); end

Pp_all = full(BMA.Pp(:));
fprintf('BMA3.Pp range = [%.4f, %.4f]\n', min(Pp_all), max(Pp_all));

hasA = any(contains(Pnames,'A(') | contains(Pnames,'A{'));
hasB = any(contains(Pnames,'B(') | contains(Pnames,'B{'));
fprintf('[Info] Pnames contains A=%d | B=%d\n', hasA, hasB);

%% ===================== Decide fields to plot from MATRIX_MODE =====================
switch char(MATRIX_MODE)
    case 'A',  fieldsToPlot = {'A'};
    case 'B',  fieldsToPlot = {'B'};
    case 'AB', fieldsToPlot = {'A','B'};
end

%% ===================== Plot =====================
for pf = 1:numel(fieldsToPlot)
    fieldChar = upper(fieldsToPlot{pf});
    if strcmp(fieldChar,'A') && ~hasA
        warning('This mat has no A parameters. Skip A.');
        continue;
    end
    if strcmp(fieldChar,'B') && ~hasB
        warning('This mat has no B parameters. Skip B.');
        continue;
    end

    field_dir = fullfile(out_dir, sprintf('Field_%s', fieldChar));
    if ~exist(field_dir,'dir'), mkdir(field_dir); end

    for ei = 1:numel(effectsToPlot)
        eff = effectsToPlot{ei};

        edges = extract_edges_polyLQ(BMA, fieldChar, eff, pp_thr);
        fprintf('>>> [%s | %s] edges=%d at PP>=%.2f\n', eff, fieldChar, numel(edges), pp_thr);

        % 无边也仍然输出 heatmap（全0），用于明确展示“无显著连接”
        W_val = zeros(nROI);

        if ~isempty(edges)
            npos = sum([edges.Ep] > 0);
            nneg = sum([edges.Ep] < 0);
            fprintf('    Ep sign: +%d  -%d  (0=%d)\n', npos, nneg, numel(edges)-npos-nneg);

            txt_file = fullfile(field_dir, sprintf('edges_%s_%s_PP%03d.txt', eff, fieldChar, round(pp_thr*100)));
            write_edge_list_txt(txt_file, edges, roi_labels, pp_thr, fieldChar, eff);

            for k = 1:numel(edges)
                i = edges(k).i; j = edges(k).j;
                w = edges(k).Ep; % signed
                if ~isfinite(w), continue; end

                if DRAW_DIRECTED
                    if abs(w) > abs(W_val(i,j)), W_val(i,j) = w; end
                else
                    if abs(w) > abs(W_val(i,j)), W_val(i,j) = w; end
                    if abs(w) > abs(W_val(j,i)), W_val(j,i) = w; end
                end
            end
        end

        % 绘图用矩阵：可选线宽统一
        if UNIFORM_EDGE_WIDTH
            W_plot = zeros(nROI);
            nz = (W_val ~= 0);
            W_plot(nz) = sign(W_val(nz)) * UNIFORM_EDGE_VALUE;
        else
            W_plot = W_val;
        end

        fprintf('    W_val range: [%.4f, %.4f]\n', min(W_val(:)), max(W_val(:)));

        % ---- BrainNet（有边才画，避免生成“空脑图”）----
        if ~isempty(edges)
            edge_file = fullfile(field_dir, sprintf('edges_%s_%s_PP%03d.edge', eff, fieldChar, round(pp_thr*100)));
            dlmwrite(edge_file, W_plot, 'delimiter','\t');

            BrainNet_MapCfg(meshFile, node_file, edge_file, cfg_file);

            fig = gcf;
            set(fig, 'Color','w', 'Renderer','opengl'); drawnow;
            remove_brainnet_colorbar_legend(fig);

            if SHOW_TEXTBOX_INFO
                titleStr = sprintf('PEB3 polyLQ | %s | %s | PP>=%.2f', eff, fieldChar, pp_thr);
                add_textbox(fig, edges, roi_labels, titleStr, fieldChar, eff, pp_thr);
            end

            out_png = fullfile(field_dir, sprintf('BrainNet_%s_%s_PP%03d.png', eff, fieldChar, round(pp_thr*100)));
            out_pdf = fullfile(field_dir, sprintf('BrainNet_%s_%s_PP%03d.pdf', eff, fieldChar, round(pp_thr*100)));
            print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi), '-opengl');
            print(fig, out_pdf, '-dpdf', '-painters');
            close(fig);
        end

        % ---- Heatmap（始终输出：即使无边也输出全0）----
        hm_png = fullfile(field_dir, sprintf('Heatmap_%s_%s_PP%03d.png', eff, fieldChar, round(pp_thr*100)));
        hm_pdf = fullfile(field_dir, sprintf('Heatmap_%s_%s_PP%03d.pdf', eff, fieldChar, round(pp_thr*100)));

        titleMain = sprintf('PEB3 polyLQ | Field %s', fieldChar);
        plot_and_save_heatmap_squarePNG(W_val, roi_labels, titleMain, eff, hm_png, hm_pdf, ...
            png_dpi, pp_thr, DEBUG_HEATMAP_MAPPING, EDGE_NEG_COLOR, EDGE_POS_COLOR);
    end
end

fprintf('[DONE] out_dir = %s\n', out_dir);
end

%% ===================== Helpers =====================

function make_bnv_cfg_edgeSign(cfg_file, mode, thr0, posColor, negColor, size_abs, color_abs, directed, edge_size_scale, ...
                              opacity_mode, opacity_same, node_rgb)
if exist(cfg_file,'file')==2, return; end

global EC
EC = [];
try
    BrainNet;  % init EC
    if ~isempty(gcf), close(gcf); end
catch
    EC = struct();
end
if ~isfield(EC,'edg') || isempty(EC.edg), EC.edg = struct(); end
if ~isfield(EC,'nod') || isempty(EC.nod), EC.nod = struct(); end

EC.edg.size_abs  = double(logical(size_abs));
EC.edg.color_abs = double(logical(color_abs));
EC.edg.directed  = double(logical(directed));
edge_size_scale = max(0.1, min(10, edge_size_scale));
EC.edg.size_ratio = edge_size_scale;

if ~isfield(EC.edg,'CM') || isempty(EC.edg.CM) || size(EC.edg.CM,2)~=3 || size(EC.edg.CM,1)<64
    EC.edg.CM = ones(64,3)*0.5;
end

mode = lower(string(mode));
switch mode
    case "threshold"
        EC.edg.color = 3;
        EC.edg.color_threshold = thr0;
        EC.edg.CM(1,:)  = posColor;
        EC.edg.CM(64,:) = negColor;
    otherwise
        EC.edg.color = 3;
        EC.edg.color_threshold = thr0;
        EC.edg.CM(1,:)  = posColor;
        EC.edg.CM(64,:) = negColor;
end

opacity_same = max(0, min(1, opacity_same));
EC.edg.opacity = double(opacity_mode);
EC.edg.opacity_same = opacity_same;
if ~isfield(EC.edg,'opacity_max'), EC.edg.opacity_max = min(0.95, opacity_same+0.25); end
if ~isfield(EC.edg,'opacity_min'), EC.edg.opacity_min = max(0.05, opacity_same-0.25); end
if ~isfield(EC.edg,'opacity_abs'), EC.edg.opacity_abs = 1; end

if ~isfield(EC.nod,'CM') || isempty(EC.nod.CM) || size(EC.nod.CM,2)~=3 || size(EC.nod.CM,1)<64
    EC.nod.CM = zeros(64,3);
end
EC.nod.CM = repmat(node_rgb(:)', 64, 1);
EC.nod.color = 2;

save(cfg_file,'EC');
end

function BMA = load_BMA3_from_PEBmat(matfile)
S = load(matfile);
if isfield(S,'res_poly') && isfield(S.res_poly,'BMA3')
    BMA = S.res_poly.BMA3; return;
end
if isfield(S,'BMA3')
    BMA = S.BMA3; return;
end
error('No BMA3 found in: %s', matfile);
end

function edges = extract_edges_polyLQ(BMA, fieldChar, eff, pp_thr)
fieldChar = upper(char(fieldChar));
eff = char(eff);

Pnames = BMA.Pnames;
if isstring(Pnames), Pnames = cellstr(Pnames); end
if ischar(Pnames),   Pnames = cellstr(Pnames); end

try
    Ep = full(spm_vec(BMA.Ep));
catch
    Ep = full(BMA.Ep(:));
end
Ep = Ep(:);
Pp = full(BMA.Pp(:)); Pp = Pp(:);

n = min([numel(Pnames), numel(Ep), numel(Pp)]);
Pnames = Pnames(1:n); Ep = Ep(1:n); Pp = Pp(1:n);

pat = sprintf('^(?<eff>[^:]+)\\s*:\\s*(?<field>%s)\\s*(?:\\{\\s*(?<u>\\d+)\\s*\\})?\\s*\\(\\s*(?<i>\\d+)\\s*,\\s*(?<j>\\d+)\\s*\\)\\s*$', fieldChar);

edges = struct('i',{},'j',{},'Ep',{},'Pp',{},'pname',{});
kout = 0;
for k = 1:n
    pn = strtrim(Pnames{k});
    tok = regexp(pn, pat, 'names', 'once');
    if isempty(tok), continue; end
    if ~strcmpi(strtrim(tok.eff), strtrim(eff)), continue; end
    if ~isfinite(Pp(k)) || Pp(k) < pp_thr, continue; end
    i = str2double(tok.i); j = str2double(tok.j);
    if ~isfinite(i) || ~isfinite(j) || i==j, continue; end
    kout = kout + 1;
    edges(kout).i = i;
    edges(kout).j = j;
    edges(kout).Ep = double(Ep(k));
    edges(kout).Pp = double(Pp(k));
    edges(kout).pname = pn;
end
end

function write_bnv_node(nodeFile, coords, labels, nodeSize)
n = size(coords,1);
if nargin < 4 || isempty(nodeSize), nodeSize = 5; end
nodeValue = (1:n)';

fid = fopen(nodeFile,'w');
assert(fid>0, 'Cannot write node file: %s', nodeFile);
for i = 1:n
    fprintf(fid, '%g %g %g %g %g %s\n', ...
        coords(i,1), coords(i,2), coords(i,3), nodeValue(i), nodeSize, labels{i});
end
fclose(fid);
end

function write_edge_list_txt(txt_file, edges, roi_labels, pp_thr, fieldChar, eff)
fid = fopen(txt_file,'w');
assert(fid>0, 'Cannot write: %s', txt_file);
fprintf(fid, 'Effect = %s\nField = %s\nPPthr = %.4f\n', eff, fieldChar, pp_thr);
fprintf(fid, 'Direction: %s(i,j) => j -> i (source=j, target=i)\n\n', fieldChar);
for k = 1:numel(edges)
    i = edges(k).i; j = edges(k).j;
    fprintf(fid, '%s -> %s | Ep=%+.4f | PP=%.4f | %s\n', ...
        roi_labels{j}, roi_labels{i}, edges(k).Ep, edges(k).Pp, edges(k).pname);
end
fclose(fid);
end

function add_textbox(fig, edges, roi_labels, titleStr, fieldChar, eff, pp_thr)
lines = {titleStr; sprintf('Effect=%s | Field=%s | PP>=%.2f', eff, fieldChar, pp_thr); ...
         sprintf('Dir: %s(i,j)=>j->i', fieldChar); ' '};
for k = 1:numel(edges)
    i = edges(k).i; j = edges(k).j;
    lines{end+1} = sprintf('%s->%s Ep=%+.3f PP=%.3f', roi_labels{j}, roi_labels{i}, edges(k).Ep, edges(k).Pp);
end
annotation(fig, 'textbox', [0.01 0.01 0.35 0.98], ...
    'String', lines, 'FitBoxToText','on', 'BackgroundColor','w', ...
    'EdgeColor',[0.2 0.2 0.2], 'Interpreter','none', 'FontSize', 9);
end

function remove_brainnet_colorbar_legend(fig)
if nargin<1 || isempty(fig), fig = gcf; end
cbs = findall(fig, 'Type', 'ColorBar');
if ~isempty(cbs), delete(cbs); end
lgs = findall(fig, 'Type', 'Legend');
if ~isempty(lgs), delete(lgs); end
end

function s = rgb_tag(rgb)
rgb = max(0,min(1,rgb(:)'));
s = sprintf('%03d%03d%03d', round(255*rgb(1)), round(255*rgb(2)), round(255*rgb(3)));
end

%% ===================== Heatmap (square, consistent) =====================
function plot_and_save_heatmap_squarePNG(W, labels, titleMain, effTag, out_png, out_pdf, png_dpi, pp_thr, debug_map, blue, red)
% row = target, col = source；左下角从 NTO 开始；像素锁定保证正方形 PNG
% [NEW] Y轴标签使用手动 text() 绘制，拉开与轴的距离

desired = {'NTO','NPO','NPC3','NPC2','NPC1','NF'};
[tf, ord] = ismember(desired, labels);
assert(all(tf), 'Heatmap reorder failed: labels must contain %s', strjoin(desired, ', '));

W2 = W(ord, ord);
lab2 = cellstr(upper(string(desired)));
n = numel(lab2);

if debug_map
    fprintf('\n[HEATMAP MAPPING] row=target, col=source\n');
    for r = 1:n
        for c = 1:n
            if W2(r,c) ~= 0
                fprintf('%s <- %s : Ep=%+.4f\n', lab2{r}, lab2{c}, W2(r,c));
            end
        end
    end
    fprintf('[HEATMAP MAPPING END]\n\n');
end

cmap = muted_blue_white_red(256, blue, [1 1 1], red);

mx = max(abs(W2(:)));
if ~isfinite(mx) || mx == 0, mx = 1; end

safeTitle = char(string(titleMain));
safeTitle = strrep(safeTitle, '%', '%%');

% --------- 像素级布局：axes 正方形，不被 colorbar 影响 ---------
mL = 260;   % left margin (Y labels)
mB = 260;   % bottom margin (X labels)
mT = 200;   % top margin (title)
gap = 60;   % gap between heatmap and colorbar
cbW = 110;  % colorbar width
mR = 160;   % right margin
side = 1800; % heatmap square side in pixels

figW = mL + side + gap + cbW + mR;
figH = mB + side + mT;

fig = figure('Color','w','Units','pixels','Position',[100 80 figW figH]);
set(fig,'InvertHardcopy','off');

ax = axes('Parent',fig,'Units','pixels','Position',[mL mB side side]);
ax.ActivePositionProperty = 'position';

imagesc(ax, W2);
colormap(ax, cmap);
caxis(ax, [-mx mx]);
set(ax,'YDir','normal');   % 左下角为起点
axis(ax,'image');
xlim(ax,[0.5 n+0.5]); ylim(ax,[0.5 n+0.5]);

% ===================== [NEW] 手动Y轴标签 =====================
set(ax, 'XTick',1:n, 'XTickLabel',lab2, ...
        'YTick',1:n, 'YTickLabel',repmat({''},1,n), ...
        'TickLabelInterpreter','none');
xtickangle(ax, 0);

ax.FontName = 'Arial';
ax.FontSize = 38;
ax.LineWidth = 2.6;
ax.TickDir = 'out';

% 手动Y标签：x位置放在 0.5 左侧；Clipping off 保证显示
x_text = 0.40; % 越小越远（更靠左）
for yy = 1:n
    text(ax, x_text, yy, lab2{yy}, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontName','Arial', 'FontSize', ax.FontSize, ...
        'Interpreter','none', 'Clipping','off');
end

% （保留原 try 块，不依赖它）
try
    ax.XAxis.TickLabelGapOffset = 12;
    ax.YAxis.TickLabelGapOffset = 12;
catch
end

box(ax,'on');
rectangle(ax, 'Position',[0.5 0.5 n n], 'EdgeColor',[0.05 0.05 0.05], 'LineWidth',3.0);

cb = colorbar(ax);
cb.Units = 'pixels';
cb.Position = [mL + side + gap, mB, cbW, side];
cb.Box = 'off';
cb.Label.String = 'Ep';
cb.FontName = 'Arial';
cb.FontSize = 28;
cb.Label.FontSize = 28;

ax.Position = [mL mB side side];

tstr = sprintf('Ep heatmap | %s | %s | PP>=%.2f', safeTitle, effTag, pp_thr);
annotation(fig,'textbox', [0.00 (mB+side+30)/figH 1.00 0.06], ...
    'String', tstr, 'EdgeColor','none', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'Interpreter','none', 'FontName','Arial', 'FontSize', 44, 'Color',[0 0 0]);

try
    exportgraphics(fig, out_png, 'Resolution', png_dpi, 'BackgroundColor','white');
catch
    set(fig,'PaperPositionMode','auto');
    print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi));
end
close(fig);

% ===================== PDF（矢量优先，同步手动Y标签） =====================
fig = figure('Color','w');
ax = axes('Parent',fig);
imagesc(ax, W2); colormap(ax,cmap); caxis(ax,[-mx mx]);
set(ax,'YDir','normal'); axis(ax,'image');
xlim(ax,[0.5 n+0.5]); ylim(ax,[0.5 n+0.5]);

set(ax,'XTick',1:n,'XTickLabel',lab2, ...
       'YTick',1:n,'YTickLabel',repmat({''},1,n), ...
       'TickLabelInterpreter','none');
xtickangle(ax,0);

ax.FontName = 'Arial';
ax.FontSize = 16;

x_text = 0.20;
for yy = 1:n
    text(ax, x_text, yy, lab2{yy}, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontName','Arial', 'FontSize', ax.FontSize, ...
        'Interpreter','none', 'Clipping','off');
end

title(ax,tstr,'Interpreter','none');
cb = colorbar(ax); cb.Label.String='Ep';

try
    exportgraphics(fig, out_pdf, 'ContentType','vector', 'BackgroundColor','white');
catch
    print(fig, out_pdf, '-dpdf', '-painters');
end
close(fig);
end

function cm = muted_blue_white_red(n, blue, white, red)
if nargin<1, n = 256; end
n = max(3, round(n));
half = floor(n/2);
t1 = linspace(0,1,half)';
t2 = linspace(0,1,n-half)';

cm1 = blue  .* (1-t1) + white .* t1;
cm2 = white .* (1-t2) + red   .* t2;
cm  = [cm1; cm2];
cm  = cm(1:n,:);
end

