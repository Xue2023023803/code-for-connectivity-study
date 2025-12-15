%% plot_PEB3_directed_graphs_PP095_MNIbg.m
% Visualize directed DCM connections (B-parameters) with PP > 0.95 on an
% anatomical background (SPM canonical avg152T1, sagittal slice).
%
% Inputs expected in result_dir:
%   - BMA3_connect.mat       (variable: BMA3)
%   - BMA3_unconnect.mat     (variable: BMA3)
%   - BMA3_condEffect.mat    (variable: BMA3)   % connect vs unconnect effect
%   - metadata_DCM_index.mat (REQUIRED; used to find DCMidx and/or roi_labels)
%
% Output:
%   - PNG/PDF figures and edge-list txt files under result_dir/figs_PPxxx_MNIbg/
%
% Notes:
%   - An arrow is drawn from SOURCE -> TARGET for parameter B(TARGET, SOURCE).
%   - For condEffect: Ep > 0 means "connect > unconnect" ONLY if your design
%     coded condEffect = connect - unconnect. If you coded the opposite, flip.
%
% Requirements:
%   - SPM on MATLAB path (spm, spm_vol, spm_read_vols).
% -------------------------------------------------------------------------

clear; clc;

%% ------------------ User settings ------------------
result_dir = '/home/xue/data/prf_result/PEB_3level_QCskip_PP095_3sets/ses-01';

pp_thr      = 0.95;   % exploratory threshold
node_radius = 6;      % shorten arrows so they don't end inside nodes

% If you want a specific sagittal slice in MNI x(mm), set slice_x_mni = 0.
slice_x_mni = 0;

% For condEffect interpretation:
%   +1 => connect > unconnect  (default)
%   -1 => unconnect > connect
condEffect_sign = +1;

%% ------------------ Locate ROI order ------------------
% IMPORTANT:
%  - Do not use default/fallback ROI order.
%  - Prefer reading roi_labels/baseROI_list/roi_list from metadata.
%  - If missing, derive ROI order from DCM.Y.name in DCM files listed by DCMidx,
%    then append-save roi_labels back to metadata.
%  - If still not found, ERROR.

meta_file = fullfile(result_dir, 'metadata_DCM_index.mat');
if exist(meta_file, 'file') ~= 2
    error('Cannot find metadata_DCM_index.mat in result_dir: %s', result_dir);
end

M = load(meta_file);

baseROI_list = {};

% (A) Try direct fields first (if your pipeline already wrote it)
if isfield(M, 'roi_labels') && ~isempty(M.roi_labels)
    baseROI_list = M.roi_labels;
elseif isfield(M, 'baseROI_list') && ~isempty(M.baseROI_list)
    baseROI_list = M.baseROI_list;
elseif isfield(M, 'roi_list') && ~isempty(M.roi_list)
    baseROI_list = M.roi_list;
end

% Normalize if found
if ~isempty(baseROI_list)
    baseROI_list = normalize_roi_list(baseROI_list);
end

% (B) If missing, derive from DCM.Y.name via DCMidx.file and save back
if isempty(baseROI_list)
    baseROI_list = derive_roi_labels_from_DCMidx_and_save(meta_file, M);
end

nBase = numel(baseROI_list);

fprintf('>>> ROI order used for plotting (index -> ROI):\n');
for i = 1:nBase
    fprintf('    %2d -> %s\n', i, baseROI_list{i});
end

%% ------------------ Load SPM template as background ------------------
if exist('spm','file') ~= 2
    error('SPM not found on MATLAB path. Please add SPM and retry.');
end

spm_dir = spm('Dir');
template_nii = fullfile(spm_dir, 'canonical', 'avg152T1.nii');
assert(exist(template_nii,'file')==2, 'Cannot find template: %s', template_nii);

V = spm_vol(template_nii);
Y = spm_read_vols(V);
dim = V.dim;

% Choose sagittal slice index based on MNI x(mm)
ix = xindex_from_mni_x(V, slice_x_mni);
ix = max(1, min(dim(1), ix));

% Node positions: (j,k) in template voxel space computed from Harvey y,z
[node_xy, roi_yz_world] = compute_roi_pixel_positions_and_worldYZ(baseROI_list, V);

% Ensure all ROI coordinates are known
if any(any(isnan(node_xy)))
    error('Some ROI names do not have predefined coordinates. Update compute_roi_pixel_positions_and_worldYZ().');
end

%% ------------------ Load BMA files & plot ------------------
fig_dir = fullfile(result_dir, sprintf('figs_PP%03d_MNIbg', round(pp_thr*100)));
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

sets = { ...
    struct('tag','connect',     'mat', fullfile(result_dir,'BMA3_connect.mat'),     'title','PEB3: connect (mean B)',     'is_condEffect',false), ...
    struct('tag','unconnect',   'mat', fullfile(result_dir,'BMA3_unconnect.mat'),   'title','PEB3: unconnect (mean B)',   'is_condEffect',false), ...
    struct('tag','condEffect',  'mat', fullfile(result_dir,'BMA3_condEffect.mat'),  'title','PEB3: condEffect (connect vs unconnect)', 'is_condEffect',true) ...
};

for si = 1:numel(sets)
    S = sets{si};

    if exist(S.mat, 'file') ~= 2
        warning('Missing file: %s (skip %s)', S.mat, S.tag);
        continue;
    end

    BMA = load_BMA3(S.mat);
    edges = extract_B_edges_from_BMA(BMA, pp_thr);

    % Write edge list
    out_txt = fullfile(fig_dir, sprintf('edges_%s_PP%03d.txt', S.tag, round(pp_thr*100)));
    write_edge_list(out_txt, edges, baseROI_list, S.is_condEffect, condEffect_sign, pp_thr);

    % Plot
    fig_title = sprintf('%s | PP>%.2f | sagittal x=%gmm', S.title, pp_thr, slice_x_mni);
    out_png = fullfile(fig_dir, sprintf('DCM_%s_PP%03d_MNIbg.png', S.tag, round(pp_thr*100)));
    out_pdf = fullfile(fig_dir, sprintf('DCM_%s_PP%03d_MNIbg.pdf', S.tag, round(pp_thr*100)));

    plot_dcm_mni_network(V, Y, ix, baseROI_list, node_xy, roi_yz_world, edges, fig_title, node_radius);

    saveas(gcf, out_png);
    try
        exportgraphics(gcf, out_pdf, 'ContentType','vector');
    catch
        % older MATLAB: ignore PDF export failure
    end
    close(gcf);

    fprintf('>>> [%s] edges=%d, figure saved to:\n    %s\n', S.tag, numel(edges), out_png);
    fprintf('    edge list: %s\n', out_txt);
end

fprintf('\nDone.\n');

%% =======================================================================
% Local functions
% =======================================================================

function roi_list = normalize_roi_list(x)
% Normalize ROI list to cell array of char, column vector.
if isstring(x)
    roi_list = cellstr(x(:));
elseif ischar(x)
    roi_list = cellstr(x);
elseif iscell(x)
    roi_list = x(:);
else
    error('ROI list type not supported. Expect cell/string/char.');
end

for i = 1:numel(roi_list)
    v = roi_list{i};
    if isstring(v); v = char(v); end
    if iscell(v) && numel(v)==1; v = v{1}; end
    if issparse(v); v = full(v); end
    if ~ischar(v); v = char(string(v)); end
    roi_list{i} = v;
end

if isempty(roi_list)
    error('ROI list is empty after normalization.');
end
end

function roi_list = derive_roi_labels_from_DCMidx_and_save(meta_file, M)
% Derive ROI order from DCM.Y.name using DCM files referenced by metadata.DCMidx.file.
% Also enforces that ROI order is identical across all included DCMs.
% Finally appends roi_labels back to metadata_DCM_index.mat.

if ~isfield(M, 'DCMidx') || isempty(M.DCMidx)
    error('metadata_DCM_index.mat lacks roi_labels/baseROI_list/roi_list and also has no DCMidx to derive ROI names.');
end

% Helper to resolve possibly-relative DCM file paths
resolve_dcm_path = @(f, entry) local_resolve_path(f, M, entry);

% 1) Pick a reference DCM file (first non-skip that exists)
ref_file = '';
ref_entry = [];
for k = 1:numel(M.DCMidx)
    entry = M.DCMidx(k);
    if isfield(entry,'skip') && entry.skip
        continue;
    end
    if ~isfield(entry,'file') || isempty(entry.file)
        continue;
    end
    f = entry.file;
    if isstring(f); f = char(f); end
    f2 = resolve_dcm_path(f, entry);
    if exist(f2, 'file') == 2
        ref_file = f2;
        ref_entry = entry; %#ok<NASGU>
        break;
    end
end

if isempty(ref_file)
    error('Cannot find any existing non-skip DCM file from metadata.DCMidx.file.');
end

S = load(ref_file, 'DCM');
if ~isfield(S,'DCM')
    error('DCM variable not found in: %s', ref_file);
end
DCM = S.DCM;

if ~isfield(DCM,'Y') || ~isfield(DCM.Y,'name') || isempty(DCM.Y.name)
    error('Cannot find ROI names in DCM.Y.name in: %s', ref_file);
end

roi_list = normalize_roi_list(DCM.Y.name);

% 2) Consistency check across ALL included DCMs
for k = 1:numel(M.DCMidx)
    entry = M.DCMidx(k);

    if isfield(entry,'skip') && entry.skip
        continue;
    end
    if ~isfield(entry,'file') || isempty(entry.file)
        continue;
    end

    f = entry.file;
    if isstring(f); f = char(f); end
    f2 = resolve_dcm_path(f, entry);

    if exist(f2, 'file') ~= 2
        continue;
    end

    T = load(f2, 'DCM');
    if ~isfield(T,'DCM') || ~isfield(T.DCM,'Y') || ~isfield(T.DCM.Y,'name') || isempty(T.DCM.Y.name)
        error('Cannot read DCM.Y.name from: %s', f2);
    end

    roi2 = normalize_roi_list(T.DCM.Y.name);

    if ~isequal(roi2, roi_list)
        error('ROI order mismatch detected!\n  Reference: %s\n  Mismatch : %s', ref_file, f2);
    end
end

% 3) Append-save to metadata for future runs
roi_labels = roi_list; %#ok<NASGU>
save(meta_file, 'roi_labels', '-append');
fprintf('>>> [ROI] Derived roi_labels from DCM.Y.name and appended to metadata:\n    %s\n', meta_file);
end

function f2 = local_resolve_path(f, M, entry)
% If f exists, use it.
% If not, try constructing full path using root_prf/session_label/dcm_subdir/subj when possible.
f2 = f;

if exist(f2,'file')==2
    return;
end

% If f is just a filename, try build:
% root_prf/sub-XX/ses-YY/<dcm_subdir>/<filename>
try
    if isfield(M,'root_prf') && isfield(M,'session_label') && isfield(M,'dcm_subdir') ...
            && isfield(entry,'subj') && ~isempty(M.root_prf) && ~isempty(M.session_label) && ~isempty(M.dcm_subdir)
        root_prf = M.root_prf; if isstring(root_prf); root_prf = char(root_prf); end
        ses      = M.session_label; if isstring(ses); ses = char(ses); end
        dcm_dir  = M.dcm_subdir; if isstring(dcm_dir); dcm_dir = char(dcm_dir); end
        subj     = entry.subj; if isstring(subj); subj = char(subj); end

        cand = fullfile(root_prf, subj, ses, dcm_dir, f);
        if exist(cand,'file')==2
            f2 = cand;
            return;
        end
    end
catch
    % ignore and keep original
end
end

function ix = xindex_from_mni_x(V, x_mni)
% Convert MNI x(mm) to template voxel x-index.
M = V.mat;
invM = inv(M);
v = invM * [x_mni; 0; 0; 1];
ix = round(v(1));
end

function BMA = load_BMA3(matfile)
% Load variable named BMA3 (preferred) or fallback to any struct with Pnames/Ep/Pp.
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
% Extract significant directed B edges from BMA based on PP threshold.
%
% Output edges is struct array with fields:
%   i (target), j (source), Ep, Pp, pname

Pnames = BMA.Pnames;
if isstring(Pnames); Pnames = cellstr(Pnames); end
if ischar(Pnames);   Pnames = cellstr(Pnames); end

% Robust Ep extraction (handles sparse)
try
    Ep = full(spm_vec(BMA.Ep));
catch
    Ep = full(BMA.Ep(:));
end
Ep = Ep(:);

Pp = full(BMA.Pp(:));
Pp = Pp(:);

n = min([numel(Pnames), numel(Ep), numel(Pp)]);

edges = struct('i', {}, 'j', {}, 'Ep', {}, 'Pp', {}, 'pname', {});
kout = 0;

for k = 1:n
    pn = Pnames{k};
    if isstring(pn); pn = char(pn); end
    if isempty(pn) || ~contains(pn, 'B(')
        continue;
    end

    tok = regexp(pn, 'B\((\d+),(\d+)\)', 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    i = str2double(tok{1});
    j = str2double(tok{2});
    if ~(isfinite(i) && isfinite(j))
        continue;
    end
    if i == j
        continue;
    end

    pp = double(full(Pp(k)));
    if ~isfinite(pp) || pp < pp_thr
        continue;
    end

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
if fid < 0
    warning('Cannot write: %s', out_txt);
    return;
end

fprintf(fid, 'PP threshold = %.4f (exploratory)\n', pp_thr);
fprintf(fid, 'Direction convention: B(target, source) => source -> target\n\n');

if isempty(edges)
    fprintf(fid, '(No edges pass threshold)\n');
    fclose(fid);
    return;
end

fprintf(fid, 'Edges (source -> target):\n');
fprintf(fid, '------------------------------------------------------------\n');

for e = 1:numel(edges)
    i = edges(e).i;  % target
    j = edges(e).j;  % source
    Ep = edges(e).Ep;
    Pp = edges(e).Pp;

    src = safe_roi_name(roi_list, j);
    tgt = safe_roi_name(roi_list, i);

    if is_condEffect
        if condEffect_sign > 0
            if Ep > 0
                dirtxt = 'connect > unconnect';
            elseif Ep < 0
                dirtxt = 'unconnect > connect';
            else
                dirtxt = 'no effect';
            end
        else
            if Ep > 0
                dirtxt = 'unconnect > connect';
            elseif Ep < 0
                dirtxt = 'connect > unconnect';
            else
                dirtxt = 'no effect';
            end
        end
        fprintf(fid, '%s -> %s | Ep=%.4f | PP=%.4f | %s\n', src, tgt, Ep, Pp, dirtxt);
    else
        fprintf(fid, '%s -> %s | Ep=%.4f | PP=%.4f\n', src, tgt, Ep, Pp);
    end
end

fclose(fid);
end

function name = safe_roi_name(roi_list, idx)
% No fallback naming: if mismatch, error (to avoid silent label errors).
if ~(idx >= 1 && idx <= numel(roi_list))
    error('ROI index %d out of range (1..%d). Check ROI order / metadata.', idx, numel(roi_list));
end

name = roi_list{idx};
if isstring(name); name = char(name); end
if iscell(name) && numel(name)==1; name = name{1}; end
if issparse(name); name = full(name); end
if ~ischar(name); name = char(string(name)); end
end

function [node_xy, roi_yz_world] = compute_roi_pixel_positions_and_worldYZ(baseROI_list, V)
% Compute node position in template voxel space using Harvey & Dumoulin (2017)
% y,z (mm) averaged across hemispheres.
%
% node_xy: [nROI x 2] in voxel indices (j,k) corresponding to (y,z)
% roi_yz_world: [nROI x 2] with [y_mm, z_mm] (NaN if unknown)

harvey_yz = struct();

% Mean (y,z) across hemispheres, from Harvey & Dumoulin 2017:
harvey_yz.NTO  = [mean([-75, -77]), mean([-4, -3])];
harvey_yz.NPO  = [mean([-82, -80]), mean([34, 32])];
harvey_yz.NPC1 = [mean([-61, -59]), mean([60, 61])];
harvey_yz.NPC2 = [mean([-40, -43]), mean([52, 48])];
harvey_yz.NPC3 = [mean([-30, -29]), mean([40, 34])];
harvey_yz.NF   = [mean([-11, -11]), mean([52, 50])];

nBase = numel(baseROI_list);
node_xy = nan(nBase, 2);
roi_yz_world = nan(nBase, 2);

M = V.mat;
a_y = M(2, 2); b_y = M(2, 4);     % y_world = a_y * j + b_y
a_z = M(3, 3); b_z = M(3, 4);     % z_world = a_z * k + b_z

for b = 1:nBase
    roi = baseROI_list{b};
    if isstring(roi); roi = char(roi); end

    if isfield(harvey_yz, roi)
        yz = harvey_yz.(roi);  % [y_mm, z_mm]
        y_w = yz(1);
        z_w = yz(2);

        j = (y_w - b_y) / a_y;
        k = (z_w - b_z) / a_z;

        node_xy(b,:) = [j, k];
        roi_yz_world(b,:) = [y_w, z_w];
    else
        error('ROI "%s" has no predefined coordinate in compute_roi_pixel_positions_and_worldYZ(). Add its (y,z) MNI coords.', roi);
    end
end
end

function plot_dcm_mni_network(V, Y, ix, baseROI_list, node_xy, roi_yz_world, edges, fig_title, node_radius)
% Plot background + directed edges + node labels.

dim = V.dim;
slice = squeeze(Y(ix, :, :));   % [Ny x Nz]

% Define colors: positive Ep (red) / negative Ep (blue)
col_pos = [1.0, 0.3, 0.3];
col_neg = [0.3, 0.8, 1.0];
col_txt = [1 1 1];

% Scale line width by |Ep|
absEp = arrayfun(@(e) abs(e.Ep), edges);
if isempty(absEp)
    max_abs = 1;
else
    max_abs = max(absEp);
    if max_abs == 0; max_abs = 1; end
end

fig = figure('Name', fig_title, 'Color', [0 0 0]);
ax  = axes('Parent', fig);
hold(ax, 'on');
axis(ax, 'image');
set(ax, 'YDir', 'normal');
set(ax, 'Color', [0 0 0]);

lo = prctile(slice(:), 5);
hi = prctile(slice(:), 95);
imagesc(slice', 'Parent', ax, [lo, hi]);
colormap(ax, gray);

xlim(ax, [1 dim(2)]);
ylim(ax, [1 dim(3)]);

% Edges
for e = 1:numel(edges)
    i = edges(e).i;  % target
    j = edges(e).j;  % source

    if i < 1 || i > numel(baseROI_list) || j < 1 || j > numel(baseROI_list)
        continue;
    end

    x1 = node_xy(j, 1); y1 = node_xy(j, 2);
    x2 = node_xy(i, 1); y2 = node_xy(i, 2);
    if any(isnan([x1 y1 x2 y2]))
        continue;
    end

    % Shorten for node radius
    [xs, ys, xt, yt] = shorten_segment(x1, y1, x2, y2, node_radius);

    dx = xt - xs;
    dy = yt - ys;

    % Style
    w = 1.2 + 3.5 * (abs(edges(e).Ep) / max_abs);
    if edges(e).Ep >= 0
        c = col_pos;
    else
        c = col_neg;
    end

    % Arrow
    quiver(ax, xs, ys, dx, dy, 0, ...
        'Color', c, ...
        'LineWidth', w, ...
        'MaxHeadSize', 0.8);
end

% Nodes + labels (with y,z in mm if available)
node_size = 140;
for b = 1:numel(baseROI_list)
    x = node_xy(b, 1);
    y = node_xy(b, 2);
    if any(isnan([x y])); continue; end

    scatter(ax, x, y, node_size, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 1.5);

    yz = roi_yz_world(b,:);
    if all(isfinite(yz))
        lab = sprintf('%s\\n(y=%d,z=%d)', baseROI_list{b}, round(yz(1)), round(yz(2)));
    else
        lab = baseROI_list{b};
    end

    text(ax, x, y - 6, lab, ...
        'Color', col_txt, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');
end

title(ax, fig_title, 'Color', [1 1 1], 'FontSize', 13, 'FontWeight', 'bold');

note = 'Arrow: source -> target for B(target,source). Color: Ep>0 (red), Ep<0 (blue).';
text(ax, 5, 8, note, 'Color', [1 1 1], 'FontSize', 9, ...
    'FontWeight', 'bold', 'Interpreter', 'none');
end

function [xs, ys, xt, yt] = shorten_segment(x1, y1, x2, y2, r)
% Shorten both ends of a segment by r (in data units), if possible.
vx = x2 - x1;
vy = y2 - y1;
d  = sqrt(vx^2 + vy^2);
if d <= 2*r || d == 0
    xs = x1; ys = y1; xt = x2; yt = y2;
    return;
end
ux = vx / d;
uy = vy / d;
xs = x1 + ux * r;
ys = y1 + uy * r;
xt = x2 - ux * r;
yt = y2 - uy * r;
end

