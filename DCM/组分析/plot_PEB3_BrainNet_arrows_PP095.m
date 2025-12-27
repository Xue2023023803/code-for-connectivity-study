%% plot_PEB3_BrainNet_PP095_withValueAndDirection.m
% BrainNet 3D brain + nodes + single-line edges (undirected),
% and annotate DCM direction + Ep + PP in a textbox on the figure.
%
% ROI order (index -> ROI):
%  1 NF, 2 NPC1, 3 NPC2, 4 NPC3, 5 NPO, 6 NTO
%
% Fix "blank saved image":
%  - BrainNet_MapCfg may create/switch figure. We save gcf AFTER plotting.
%  - Use drawnow + print -opengl.

clear; clc; close all;

%% ------------------ User settings ------------------
% Directory where BMA3_*.mat live
result_dir_in = '/home/xue/data/prf_result/PEB_3level_QCskip_PP095_3sets/ses-01';  % <<< 改这里

% Output directory
out_dir = fullfile(result_dir_in, 'plot_BrainNV');
paramField = 'A';   % 'A' or 'B'   （跑A就设'A'，跑B就设'B'）

pp_thr  = 0.95;

% Optional saved BrainNet cfg (*.mat). Leave empty if your BrainNet works without cfg.
cfgFile = '';

% For condEffect interpretation text (optional)
% +1 means Ep>0 corresponds to "connect > unconnect" if your contrast is connect-unconnect.
condEffect_sign = +1;

% Save resolution
png_dpi = 300;

%% ------------------ Auto-detect BrainNet Viewer root from MATLAB path ------------------
bn_mapcfg = which('BrainNet_MapCfg');
assert(~isempty(bn_mapcfg), 'BrainNet_MapCfg not found. BrainNet Viewer is not on MATLAB path.');

BNV_DIR = fileparts(bn_mapcfg);
while ~exist(fullfile(BNV_DIR,'Data'),'dir')
    parent = fileparts(BNV_DIR);
    if strcmp(parent, BNV_DIR)
        error('Cannot find BrainNet "Data" folder by walking up from: %s', bn_mapcfg);
    end
    BNV_DIR = parent;
end
addpath(genpath(BNV_DIR));
fprintf('BrainNet root detected: %s\n', BNV_DIR);

%% ------------------ Locate BrainNet mesh (your install: Data/SurfTemplate) ------------------
meshFile = fullfile(BNV_DIR,'Data','SurfTemplate','BrainMesh_ICBM152Right_smoothed.nv');
if exist(meshFile,'file') ~= 2
    meshFile = fullfile(BNV_DIR,'Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv');
end
if exist(meshFile,'file') ~= 2
    mesh_list = dir(fullfile(BNV_DIR,'Data','**','BrainMesh*_smoothed*.nv'));
    assert(~isempty(mesh_list), 'No BrainMesh*_smoothed*.nv found under %s', fullfile(BNV_DIR,'Data'));
    meshFile = fullfile(mesh_list(1).folder, mesh_list(1).name);
end
fprintf('Using mesh: %s\n', meshFile);

%% ------------------ ROI order & MNI coords (RH) ------------------
roi_labels = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};

roi_mni = [ ...
     24  -11  52;   % NF
     22  -61  60;   % NPC1
     33  -40  52;   % NPC2
     45  -30  40;   % NPC3
     25  -82  34;   % NPO
     44  -75  -4];  % NTO

nROI = numel(roi_labels);

%% ------------------ Output dir ------------------
fig_dir = fullfile(out_dir, sprintf('figs_PP%03d_BrainNet', round(pp_thr*100)));
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% ------------------ Write node file ------------------
node_file = fullfile(fig_dir, 'nodes_numerosity_6ROIs_RH.node');
write_bnv_node(node_file, roi_mni, roi_labels);

%% ------------------ Resolve BMA3 files robustly ------------------
file_connect    = find_bma_file(result_dir_in, 'BMA3_connect.mat');
file_unconnect  = find_bma_file(result_dir_in, 'BMA3_unconnect.mat');
file_condEffect = find_bma_file(result_dir_in, 'BMA3_condEffect.mat');

fprintf('BMA3 file resolved:\n');
fprintf('  connect   : %s\n', file_connect);
fprintf('  unconnect : %s\n', file_unconnect);
fprintf('  condEffect: %s\n', file_condEffect);

sets = { ...
    struct('tag','connect',    'mat', file_connect,    'title','PEB3: connect (mean B)',    'is_condEffect',false), ...
    struct('tag','unconnect',  'mat', file_unconnect,  'title','PEB3: unconnect (mean B)',  'is_condEffect',false), ...
    struct('tag','condEffect', 'mat', file_condEffect, 'title','PEB3: condEffect (connect vs unconnect)', 'is_condEffect',true) ...
};

for si = 1:numel(sets)
    S = sets{si};

    if isempty(S.mat) || exist(S.mat,'file') ~= 2
        warning('Missing file (skip %s)', S.tag);
        continue;
    end

    BMA   = load_BMA3(S.mat);
    edges = extract_B_edges_from_BMA(BMA, pp_thr);

    % Skip if no edges
    if isempty(edges)
        fprintf('>>> [%s] edges=0 at PP>%.2f -> SKIP plotting.\n', S.tag, pp_thr);
        continue;
    end

    % Save edge list txt (direction + value)
    out_txt = fullfile(fig_dir, sprintf('edges_%s_PP%03d.txt', S.tag, round(pp_thr*100)));
    write_edge_list(out_txt, edges, roi_labels, S.is_condEffect, condEffect_sign, pp_thr);

    % Build symmetric edge matrix for BrainNet (single line)
    W = zeros(nROI);
    for e = 1:numel(edges)
        i = edges(e).i;  % target
        j = edges(e).j;  % source
        w = abs(edges(e).Ep);
        W(i,j) = max(W(i,j), w);
        W(j,i) = max(W(j,i), w);
    end
    edge_file = fullfile(fig_dir, sprintf('edges_%s_PP%03d.edge', S.tag, round(pp_thr*100)));
    dlmwrite(edge_file, W, 'delimiter','\t');

    % ---- Plot with BrainNet (do NOT create a figure yourself) ----
    call_brainnet(meshFile, node_file, edge_file, cfgFile);

    % The BrainNet figure is now current:
    fig = gcf;
    set(fig, 'Color', 'w', 'Renderer', 'opengl');
    drawnow;

    % Add a textbox with direction + Ep + PP on the figure
    add_edge_textbox(fig, edges, roi_labels, S.is_condEffect, condEffect_sign, pp_thr, S.title);
    drawnow;

    % Save (avoid blank by saving the BrainNet figure handle after drawnow)
    out_png = fullfile(fig_dir, sprintf('DCM_%s_PP%03d_BrainNet.png', S.tag, round(pp_thr*100)));
    out_pdf = fullfile(fig_dir, sprintf('DCM_%s_PP%03d_BrainNet.pdf', S.tag, round(pp_thr*100)));

    print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi), '-opengl');
    try
        print(fig, out_pdf, '-dpdf', sprintf('-r%d', png_dpi), '-opengl');
    catch
    end

    close(fig);

    fprintf('>>> [%s] edges=%d, saved:\n    %s\n', S.tag, numel(edges), out_png);
end

fprintf('Done. Output dir:\n  %s\n', fig_dir);

%% =======================================================================
% Local functions
% =======================================================================

function f = find_bma_file(rootdir, fname)
f = fullfile(rootdir, fname);
if exist(f,'file')==2, return; end

d = dir(fullfile(rootdir, '**', fname));
if ~isempty(d)
    f = fullfile(d(1).folder, d(1).name);
    return;
end

parent = fileparts(rootdir);
if ~isempty(parent) && ~strcmp(parent, rootdir)
    d2 = dir(fullfile(parent, '**', fname));
    if ~isempty(d2)
        f = fullfile(d2(1).folder, d2(1).name);
        return;
    end
end
f = '';
end

function call_brainnet(meshFile, nodeFile, edgeFile, cfgFile)
if exist('BrainNet_MapCfg','file') ~= 2
    error('BrainNet_MapCfg not found. Check BrainNet Viewer path.');
end

if ~isempty(cfgFile) && exist(cfgFile,'file')==2
    try
        BrainNet_MapCfg(meshFile, nodeFile, edgeFile, cfgFile);
        return;
    catch
    end
end

try
    BrainNet_MapCfg(meshFile, nodeFile, edgeFile);
catch
    error(['BrainNet_MapCfg failed without cfg. ' ...
           'Open BrainNetViewer GUI once, save an option cfg (*.mat), ' ...
           'then set cfgFile to that path.']);
end
end

function write_bnv_node(nodeFile, coords, labels)
% .node format: x y z color size label
n = size(coords,1);
colorIdx = ones(n,1);
nodeSize = 4*ones(n,1);

fid = fopen(nodeFile,'w');
if fid < 0; error('Cannot write node file: %s', nodeFile); end
for i = 1:n
    fprintf(fid, '%g %g %g %g %g %s\n', ...
        coords(i,1), coords(i,2), coords(i,3), ...
        colorIdx(i), nodeSize(i), labels{i});
end
fclose(fid);
end

function add_edge_textbox(fig, edges, roi_labels, is_condEffect, condEffect_sign, pp_thr, titleStr)
lines = {};
lines{end+1} = sprintf('%s', titleStr);
lines{end+1} = sprintf('PP threshold = %.2f', pp_thr);
lines{end+1} = 'Direction: B(target,source) => source → target';
lines{end+1} = ' ';

for e = 1:numel(edges)
    i  = edges(e).i;  % target
    j  = edges(e).j;  % source
    Ep = edges(e).Ep;
    Pp = edges(e).Pp;

    src = roi_labels{j};
    tgt = roi_labels{i};

    if is_condEffect
        Ep2 = condEffect_sign * Ep;
        if Ep2 > 0
            tag = 'connect > unconnect';
        elseif Ep2 < 0
            tag = 'unconnect > connect';
        else
            tag = 'no effect';
        end
        lines{end+1} = sprintf('%s → %s : Ep=%.4f, PP=%.3f (%s)', src, tgt, Ep, Pp, tag);
    else
        lines{end+1} = sprintf('%s → %s : Ep=%.4f, PP=%.3f', src, tgt, Ep, Pp);
    end
end

txt = strjoin(lines, newline);

% Put textbox at upper-left; adjust if you want
annotation(fig, 'textbox', [0.02 0.65 0.36 0.33], ...
    'String', txt, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', [0 0 0], ...
    'BackgroundColor', [1 1 1], ...
    'Interpreter', 'none', ...
    'FontSize', 10);
end

function BMA = load_BMA3(matfile)
S = load(matfile);
if isfield(S, 'BMA3')
    BMA = S.BMA3;
    return;
end
fn = fieldnames(S);
for i = 1:numel(fn)
    v = S.(fn{i});
    if isstruct(v) && isfield(v, 'Pnames') && isfield(v, 'Ep') && isfield(v, 'Pp')
        BMA = v;
        return;
    end
end
error('No BMA3-like variable found in: %s', matfile);
end

function edges = extract_B_edges_from_BMA(BMA, pp_thr)
Pnames = BMA.Pnames;
if isstring(Pnames); Pnames = cellstr(Pnames); end
if ischar(Pnames);   Pnames = cellstr(Pnames); end

try
    Ep = full(spm_vec(BMA.Ep));
catch
    Ep = full(BMA.Ep(:));
end
Ep = Ep(:);
Pp = full(BMA.Pp(:)); Pp = Pp(:);

n = min([numel(Pnames), numel(Ep), numel(Pp)]);

edges = struct('i', {}, 'j', {}, 'Ep', {}, 'Pp', {}, 'pname', {});
kout = 0;

for k = 1:n
    pn = Pnames{k};
    if isstring(pn); pn = char(pn); end
    if isempty(pn) || ~contains(pn, 'B('); continue; end

    tok = regexp(pn, 'B\((\d+),(\d+)\)', 'tokens', 'once');
    if isempty(tok); continue; end

    i = str2double(tok{1}); % target
    j = str2double(tok{2}); % source
    if ~(isfinite(i) && isfinite(j)); continue; end
    if i == j; continue; end

    pp = double(full(Pp(k)));
    if ~isfinite(pp) || pp < pp_thr; continue; end

    kout = kout + 1;
    edges(kout).i = i;
    edges(kout).j = j;
    edges(kout).Ep = double(full(Ep(k)));
    edges(kout).Pp = pp;
    edges(kout).pname = pn;
end
end

function write_edge_list(out_txt, edges, roi_list, is_condEffect, condEffect_sign, pp_thr)
fid = fopen(out_txt, 'w');
if fid < 0; warning('Cannot write: %s', out_txt); return; end

fprintf(fid, 'PP threshold = %.4f\n', pp_thr);
fprintf(fid, 'Direction: B(target, source) => source -> target\n\n');

if isempty(edges)
    fprintf(fid, '(No edges pass threshold)\n');
    fclose(fid);
    return;
end

for e = 1:numel(edges)
    i  = edges(e).i;  % target
    j  = edges(e).j;  % source
    Ep = edges(e).Ep;
    Pp = edges(e).Pp;

    src = roi_list{j};
    tgt = roi_list{i};

    if is_condEffect
        Ep2 = condEffect_sign * Ep;
        if Ep2 > 0
            dirtxt = 'connect > unconnect';
        elseif Ep2 < 0
            dirtxt = 'unconnect > connect';
        else
            dirtxt = 'no effect';
        end
        fprintf(fid, '%s -> %s | Ep=%.4f | PP=%.4f | %s\n', src, tgt, Ep, Pp, dirtxt);
    else
        fprintf(fid, '%s -> %s | Ep=%.4f | PP=%.4f\n', src, tgt, Ep, Pp);
    end
end

fclose(fid);
end

