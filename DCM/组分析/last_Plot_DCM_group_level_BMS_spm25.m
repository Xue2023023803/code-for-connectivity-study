%% DCM_group_level_visualization_spm25_Yname.m
%  用 PEB-BMA 的 B 参数画 ROI-ROI 网络图
%  ROI 名从个体 DCM 的 DCM.Y.name 读取（因为你是手动 build 的 DCM，没有 xY）

clear; clc; close all;

%% -------------------- 结果文件路径 --------------------
res_file = '/home/xue/data/prf_result/DCM_group_BMS/DCM_BMS_simple_ses-01_spm25.mat';

fprintf('>>> 读取结果文件: %s\n', res_file);
if ~exist(res_file, 'file')
    error('结果文件不存在：%s', res_file);
end

S = load(res_file);

PEB_results        = S.PEB_results;
conds              = S.conds;
subjects           = S.subjects;
model_names        = S.model_names;
best_model_idx_all = S.best_model_idx_all;

% 这些在你的计算脚本里是一起 save 出来的
if isfield(S, 'root_dir')
    root_dir = S.root_dir;
else
    root_dir = '/home/xue/data/prf_result';
end
if isfield(S, 'analysis_dirname')
    analysis_dirname = S.analysis_dirname;
else
    analysis_dirname = 'DCM_manual_spm25_multiModel';
end
if isfield(S, 'session_label')
    session_label = S.session_label;
else
    session_label = 'ses-01';
end

fprintf('\n================= Visualization: ROI-to-ROI connections (PEB-BMA, PP>0.95) =================\n');

if exist('digraph','class') ~= 8
    fprintf('[Plot] 当前 MATLAB 版本不支持 digraph，退出。\n');
    return;
end

%% -------------------- 全局 ROI label，占用一次 --------------------
roi_labels_global = {};   % 先空着，等知道 nROI 后再从 DCM.Y.name 里读

thrPP = 0.85;   % PP 阈值

for ic = 1:numel(conds)
    cond = conds{ic};
    fprintf('\n[Plot] Condition: %s (PEB-BMA B, PP > %.2f)\n', cond, thrPP);

    res = PEB_results{ic};
    if isempty(res) || ~isfield(res,'BMA') || isempty(res.BMA) ...
            || ~isfield(res,'PP')  || isempty(res.PP) ...
            || ~isfield(res,'Pnames') || isempty(res.Pnames) ...
            || ~isfield(res,'Ep_vec') || isempty(res.Ep_vec)
        fprintf('  [Plot] cond=%s: PEB/BMA 结果不完整，跳过。\n', cond);
        continue;
    end

    Pnames = res.Pnames(:);
    Pp_vec = res.PP(:);
    Ep_vec = res.Ep_vec(:);

    if numel(Pnames) ~= numel(Pp_vec) || numel(Pnames) ~= numel(Ep_vec)
        fprintf('  [Plot] cond=%s: Pnames / PP / Ep_vec 长度不一致，跳过。\n', cond);
        continue;
    end

    %% ---- 1) 从 B(i,j) 的名字推断 ROI 数 ----
    max_idx = 0;
    for k = 1:numel(Pnames)
        name_k = Pnames{k};
        if isempty(name_k) || name_k(1) ~= 'B'
            continue;
        end
        tok = regexp(name_k, '^B\((\d+),(\d+)(?:,(\d+))?\)$', 'tokens', 'once');
        if isempty(tok), continue; end
        ii = str2double(tok{1});   % 行 target
        jj = str2double(tok{2});   % 列 source
        if ~isnan(ii), max_idx = max(max_idx, ii); end
        if ~isnan(jj), max_idx = max(max_idx, jj); end
    end

    if max_idx == 0
        fprintf('  [Plot] cond=%s: Pnames 中没有任何 B(i,j) 参数，跳过。\n', cond);
        continue;
    end
    nROI = max_idx;

    %% ---- 2) 如果还没 ROI label，就从某个 DCM.Y.name 里读一次 ----
    if isempty(roi_labels_global)
        roi_labels_global = local_get_roi_labels_from_DCM_Yname( ...
            root_dir, analysis_dirname, session_label, subjects, nROI);
    end

    %% ---- 3) 针对当前条件，构造本地 roi_labels（与 nROI 对齐） ----
    if ~isempty(roi_labels_global)
        roi_labels = roi_labels_global;
        if numel(roi_labels) < nROI
            roi_labels(numel(roi_labels)+1:nROI) = ...
                arrayfun(@(ii) sprintf('ROI%d',ii), ...
                         numel(roi_labels)+1:nROI, 'UniformOutput', false);
        elseif numel(roi_labels) > nROI
            roi_labels = roi_labels(1:nROI);
        end
    else
        fprintf('  [Plot] 无法从 DCM 读取 ROI 名，使用默认 ROI1..ROI%d。\n', nROI);
        roi_labels = arrayfun(@(ii) sprintf('ROI%d',ii), 1:nROI, 'UniformOutput', false);
    end

    %% ---- 4) 把 Ep / PP 映射回 ROI×ROI 矩阵 ----
    Ep_B = zeros(nROI, nROI);
    Pp_B = zeros(nROI, nROI);

    for k = 1:numel(Pnames)
        name_k = Pnames{k};
        if isempty(name_k) || name_k(1) ~= 'B'
            continue;
        end
        tok = regexp(name_k, '^B\((\d+),(\d+)(?:,(\d+))?\)$', 'tokens', 'once');
        if isempty(tok), continue; end

        i = str2double(tok{1});   % 行 target
        j = str2double(tok{2});   % 列 source
        if isnan(i) || isnan(j) || i<1 || j<1 || i>nROI || j>nROI
            continue;
        end

        if abs(Ep_vec(k)) > abs(Ep_B(i,j))
            Ep_B(i,j) = Ep_vec(k);
            Pp_B(i,j) = Pp_vec(k);
        end
    end

    sig_mask = (Pp_B > thrPP);
    sig_mask(1:nROI+1:end) = false;  % 去掉对角线

    if ~any(sig_mask(:))
        fprintf('  [Plot] cond=%s: 没有 PP > %.2f 的 B 连接。\n', cond, thrPP);
        continue;
    end

    [from_idx, to_idx] = find(sig_mask);
    fprintf('  [Plot] cond=%s, PP > %.2f 的 B 连接（行=target, 列=source，source -> target）：\n', ...
        cond, thrPP);
    for e = 1:numel(from_idx)
        i = from_idx(e);
        j = to_idx(e);
        fprintf('     %s <- %s: Ep=%.4f, PP=%.3f\n', ...
            roi_labels{i}, roi_labels{j}, Ep_B(i,j), Pp_B(i,j));
    end

    %% ---- 5) 构造 digraph 并画图 ----
    % digraph(A) 里 A(u,v) ~= 0 表示 u -> v
    % B(i,j) 是 source=j -> target=i，所以 A(j,i) = sig_mask(i,j)
    A = sig_mask.';   % A(source,target)
    G = digraph(A);

    figName = sprintf('PEB-BMA B (PP>%.2f) - %s', thrPP, cond);
    figure('Name', figName, 'Color', 'w');
    h = plot(G, 'Layout','layered', 'NodeLabel',roi_labels);
    title(figName, 'Interpreter','none');

    % 线宽按 |Ep| 缩放
    if numedges(G) > 0
        w = zeros(numedges(G),1);
        for e = 1:numedges(G)
            s = G.Edges.EndNodes(e,1);  % source
            t = G.Edges.EndNodes(e,2);  % target
            w(e) = abs(Ep_B(t,s));
        end
        if any(w > 0)
            h.LineWidth = 1 + 4*w./max(w);
        end
    end
end

fprintf('\n可视化脚本结束。\n');

%% ======================================================================
%  辅助函数：从任意一个 DCM 的 DCM.Y.name 中读取 ROI 名
% ======================================================================
function roi_labels = local_get_roi_labels_from_DCM_Yname(root_dir, analysis_dirname, session_label, subjects, nROI)
roi_labels = {};
fprintf('[Plot/ROI] 尝试从 DCM.Y.name 中读取 ROI label...\n');

for is = 1:numel(subjects)
    subj    = subjects{is};
    dcm_dir = fullfile(root_dir, subj, session_label, analysis_dirname);
    if ~exist(dcm_dir, 'dir')
        fprintf('  [Plot/ROI] subj=%s: 目录不存在: %s\n', subj, dcm_dir);
        continue;
    end

    files = dir(fullfile(dcm_dir, '*.mat'));
    fprintf('  [Plot/ROI] subj=%s: 在 %s 下找到 %d 个 .mat 文件。\n', ...
        subj, dcm_dir, numel(files));
    if isempty(files), continue; end

    for k = 1:numel(files)
        fname    = files(k).name;
        fullpath = fullfile(dcm_dir, fname);
        try
            S = load(fullpath, 'DCM');
        catch
            fprintf('    [Plot/ROI] 载入失败: %s\n', fullpath);
            continue;
        end
        if ~isfield(S, 'DCM')
            continue;
        end
        DCM = S.DCM;

        % 关键：你是手动 build 的 DCM，真正的 ROI 名在 DCM.Y.name 里
        if ~isfield(DCM, 'Y') || ~isfield(DCM.Y, 'name') || isempty(DCM.Y.name)
            continue;
        end

        tmp = DCM.Y.name;
        if ischar(tmp)
            tmp = cellstr(tmp);
        end
        tmp = tmp(:);

        if numel(tmp) ~= nROI
            fprintf('    [Plot/ROI] %s: DCM.Y.name 数量(%d) != nROI(%d)，跳过。\n', ...
                fullpath, numel(tmp), nROI);
            continue;
        end

        roi_labels = tmp;
        fprintf('  [Plot/ROI] 成功从 %s 读取到 ROI 名:\n', fullpath);
        disp(roi_labels(:)');
        return;
    end
end

fprintf('  [Plot/ROI] 未能从任何 DCM.Y.name 中获得 ROI 名，将由主程序退回 ROI1..ROIx。\n');
end


