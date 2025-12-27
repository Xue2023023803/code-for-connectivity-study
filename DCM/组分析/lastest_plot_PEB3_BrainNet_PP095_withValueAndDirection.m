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

pp_thr  = 0.95;

SHOW_LEGEND_INFO = false;   % true=输出图包含textbox说明；false=不包含


% Choose which DCM parameter field(s) to visualize:
%   {'B'}        : modulatory connections (default)
%   {'A'}        : intrinsic/baseline connections
%   {'A','B'}    : run both (creates separate output folders)
paramFields = {'B'};

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

%% ------------------ Resolve BMA3 files robustly ------------------
file_connect    = find_bma_file(result_dir_in, 'BMA3_connect.mat');
file_unconnect  = find_bma_file(result_dir_in, 'BMA3_unconnect.mat');
file_condEffect = find_bma_file(result_dir_in, 'BMA3_condEffect.mat');

fprintf('BMA3 file resolved:\n');
fprintf('  connect   : %s\n', file_connect);
fprintf('  unconnect : %s\n', file_unconnect);
fprintf('  condEffect: %s\n', file_condEffect);

%% ================== Main plotting (loop over A/B) ==================
for pf = 1:numel(paramFields)
    paramField = paramFields{pf};

    %% ------------------ Output dir ------------------
    fig_dir = fullfile(out_dir, sprintf('figs_PP%03d_BrainNet_%s', round(pp_thr*100), paramField));
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    %% ------------------ Write node file ------------------
    node_file = fullfile(fig_dir, sprintf('nodes_numerosity_6ROIs_RH_%s.node', paramField));
    write_bnv_node(node_file, roi_mni, roi_labels);

    %% ------------------ Which sets to plot ------------------
    sets = { ...
        struct('tag','connect',    'mat', file_connect,    'title',sprintf('PEB3: connect (mean %s)', paramField),    'is_condEffect',false), ...
        struct('tag','unconnect',  'mat', file_unconnect,  'title',sprintf('PEB3: unconnect (mean %s)', paramField),  'is_condEffect',false), ...
        struct('tag','condEffect', 'mat', file_condEffect, 'title',sprintf('PEB3: condEffect (connect vs unconnect) [%s]', paramField), 'is_condEffect',true) ...
    };

    for si = 1:numel(sets)
        S = sets{si};

        if isempty(S.mat) || exist(S.mat,'file') ~= 2
            warning('Missing file (skip %s)', S.tag);
            continue;
        end

        BMA   = load_BMA3(S.mat);
        edges = extract_edges_from_BMA(BMA, pp_thr, paramField);

        % Skip if no edges
        if isempty(edges)
            fprintf('>>> [%s | %s] edges=0 at PP>%.2f -> SKIP plotting.\n', S.tag, paramField, pp_thr);
            continue;
        end

        % Save edge list txt (direction + value)
        out_txt = fullfile(fig_dir, sprintf('edges_%s_%s_PP%03d.txt', S.tag, paramField, round(pp_thr*100)));
        write_edge_list(out_txt, edges, roi_labels, S.is_condEffect, condEffect_sign, pp_thr, paramField);

        % Build symmetric edge matrix for BrainNet (single line)
        W = zeros(nROI);
        for e = 1:numel(edges)
            i = edges(e).i;  % target
            j = edges(e).j;  % source
            w = abs(edges(e).Ep);
            W(i,j) = max(W(i,j), w);
            W(j,i) = max(W(j,i), w);
        end
        edge_file = fullfile(fig_dir, sprintf('edges_%s_%s_PP%03d.edge', S.tag, paramField, round(pp_thr*100)));
        dlmwrite(edge_file, W, 'delimiter','\t');

        % ---- Plot with BrainNet (do NOT create a figure yourself) ----
        call_brainnet(meshFile, node_file, edge_file, cfgFile);

        % The BrainNet figure is now current:
        fig = gcf;
        set(fig, 'Color', 'w', 'Renderer', 'opengl');

        % Add textbox annotation with direction & values
        if SHOW_LEGEND_INFO
            add_edge_textbox(fig, edges, roi_labels, S.is_condEffect, condEffect_sign, pp_thr, S.title, paramField);
        end

        drawnow;

        % Save
        out_png = fullfile(fig_dir, sprintf('BrainNet_%s_%s_PP%03d.png', S.tag, paramField, round(pp_thr*100)));
        out_pdf = fullfile(fig_dir, sprintf('BrainNet_%s_%s_PP%03d.pdf', S.tag, paramField, round(pp_thr*100)));

        % Use print to avoid blank saved
        try
            print(fig, out_png, '-dpng', sprintf('-r%d', png_dpi), '-opengl');
            print(fig, out_pdf, '-dpdf', '-painters');
        catch
            saveas(fig, out_png);
            saveas(fig, out_pdf);
        end

        fprintf('Saved: %s\n', out_png);
        fprintf('Saved: %s\n', out_pdf);

        close(fig);
    end
end  % for pf (paramFields)

%% ========================= Helper functions =========================
function call_brainnet(meshFile, nodeFile, edgeFile, cfgFile)
% Call BrainNet_MapCfg with optional cfg
if ~isempty(cfgFile) && exist(cfgFile,'file') == 2
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

function add_edge_textbox(fig, edges, roi_labels, is_condEffect, condEffect_sign, pp_thr, titleStr, paramField)
if nargin < 8 || isempty(paramField)
    paramField = 'B';
end
lines = {};
lines{end+1} = sprintf('%s', titleStr);
lines{end+1} = sprintf('PP threshold = %.2f', pp_thr);
lines{end+1} = sprintf('Direction: %s(target,source) => source -> target', paramField);
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
            dirtxt = 'connect > unconnect';
        elseif Ep2 < 0
            dirtxt = 'unconnect > connect';
        else
            dirtxt = 'no effect';
        end
        lines{end+1} = sprintf('%s -> %s | Ep=%.3f | PP=%.3f | %s', src, tgt, Ep, Pp, dirtxt);
    else
        lines{end+1} = sprintf('%s -> %s | Ep=%.3f | PP=%.3f', src, tgt, Ep, Pp);
    end
end

annotation(fig, 'textbox', [0.01 0.01 0.35 0.98], ...
    'String', lines, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.2 0.2 0.2], ...
    'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', 9);
end

function f = find_bma_file(result_dir_in, fname)
% Find a file possibly located in result_dir_in or its subfolders
cand = fullfile(result_dir_in, fname);
if exist(cand,'file') == 2
    f = cand;
    return;
end
d = dir(fullfile(result_dir_in, '**', fname));
if ~isempty(d)
    f = fullfile(d(1).folder, d(1).name);
else
    f = '';
end
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

function edges = extract_edges_from_BMA(BMA, pp_thr, paramField)
% Extract directed edges (i=target, j=source) for a given parameter field
% from a BMA/PEB result structure (supports 'A' or 'B', including B{1}(i,j)).
%
% Inputs:
%   BMA        : struct with fields Pnames, Ep, Pp (SPM PEB/BMA output)
%   pp_thr     : posterior probability threshold (e.g., 0.95)
%   paramField : 'A' or 'B'
%
% Output:
%   edges(e).i (target), edges(e).j (source), edges(e).Ep, edges(e).Pp, edges(e).pname

if nargin < 3 || isempty(paramField)
    paramField = 'B';
end
if isstring(paramField); paramField = char(paramField); end
paramField = upper(paramField);

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

% Match: A(i,j) or B(i,j) or B{1}(i,j)
pat = sprintf('%s(?:\\{\\d+\\})?\\((\\d+),(\\d+)\\)', paramField);

for k = 1:n
    pn = Pnames{k};
    if isstring(pn); pn = char(pn); end
    if isempty(pn); continue; end

    tok = regexp(pn, pat, 'tokens', 'once');
    if isempty(tok); continue; end

    i = str2double(tok{1}); % target
    j = str2double(tok{2}); % source
    if ~(isfinite(i) && isfinite(j)); continue; end
    if i == j; continue; end

    pp = double(full(Pp(k)));
    if pp < pp_thr; continue; end

    kout = kout + 1;
    edges(kout).i = i;
    edges(kout).j = j;
    edges(kout).Ep = double(full(Ep(k)));
    edges(kout).Pp = pp;
    edges(kout).pname = pn;
end
end

function write_edge_list(out_txt, edges, roi_list, is_condEffect, condEffect_sign, pp_thr, paramField)
if nargin < 7 || isempty(paramField)
    paramField = 'B';
end
fid = fopen(out_txt, 'w');
if fid < 0; warning('Cannot write: %s', out_txt); return; end

fprintf(fid, 'PP threshold = %.4f\n', pp_thr);
fprintf(fid, 'Direction: %s(target, source) => source -> target\n\n', paramField);

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
