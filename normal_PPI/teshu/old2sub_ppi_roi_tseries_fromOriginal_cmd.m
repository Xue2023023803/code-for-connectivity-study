%% ppi_roi_tseries_fromOriginal_cmd_old2subs.m
% 只针对旧实验的 sub-02（保持 sub-02）和 sub-03（在输出中作为 sub-08）
% 从 Gray/Original/TSeries 中提取多 run 的 ROI 原始时序，用于后续 PPI
%
% 主要改动：
%   - 手动指定每个被试需要的 Scan（connect / unconnect），不再使用“前半/后半+奇偶”规则；
%   - pRF dataType 使用旧工程命名：
%         Averages_single_layer_01      -> unconnect
%         Averages_connected_layer_01   -> connect
%   - ROI 只保留以 L/R/l/r 开头的文件，大小写不敏感；
%   - sub-02 的输出仍为 sub-02，sub-03 在输出和 meta 中记为 sub-08。

clear; clc;

%% ------------ 全局路径与被试配置 ------------
% root_old：sub-02 / sub-03 顶层所在目录（请按实际情况修改）
% 例如当前路径类似：
%   /home/xue/daifenpei/sub-02/results1/sub-02/ses-01/...
% 则 root_old = '/home/xue/daifenpei';
root_old = '/home/xue/daifenpei';
ses      = 'ses-01';

% 两个被试的配置：
% - dirID   : 实际文件夹名（用于拼路径）
% - labelID : 输出和 meta 中使用的被试名
% - connectScans / unconnectScans : 需要使用的 Scan 编号
subjects = struct([]);

% 被试 3：在输出中作为 sub-08
subjects(1).dirID          = 'sub-03';
subjects(1).labelID        = 'sub-08';
subjects(1).connectScans   = [2 4 9 11 15 16 20 22 27];
subjects(1).unconnectScans = [3 5 7 10 14 18 21 23 25];


%% ------------ 循环处理每个被试 ------------
for sIdx = 1:numel(subjects)
    dirID          = subjects(sIdx).dirID;      % 实际目录名
    subj_label     = subjects(sIdx).labelID;    % 输出用的被试名
    connectScans   = subjects(sIdx).connectScans(:)';   % row
    unconnectScans = subjects(sIdx).unconnectScans(:)'; % row

    fprintf('\n==================== Subject %s (label=%s) ====================\n', ...
        dirID, subj_label);

    %% 0. 基本配置
    % 会话目录：root_old/sub-XX/results1/sub-XX/ses-01
    sess_dir  = fullfile(root_old, dirID, 'results1', dirID, ses);
    gray_dir  = fullfile(sess_dir, 'Gray');

    % 两套 pRF dataType：旧工程命名
    % Averages_single_layer_01      -> 视为 unconnect
    % Averages_connected_layer_01   -> 视为 connect
    dt_prf_list = { ...
        'Averages_single_layer_01', ...     % unconnect pRF
        'Averages_connected_layer_01' ...   % connect pRF
    };
    scanNum_prf  = 1;      % 两套 pRF 的 scan 号通常都是 Scan1

    % 原始 BOLD dataType
    dt_raw    = 'Original';

    % TR 和每个 run 的 TR 数（已切掉 NSS）
    TR          = 1.950;   % s
    NTR_expect  = 176;     % 每个 run 的帧数（已切掉 dummy NSS 之后的长度）

    % LN* 数量 ROI
    roi_glob   = fullfile(gray_dir, 'ROIs', '*N*.mat');

    % 输出目录
    out_dir = fullfile(sess_dir, 'results_ppi_roi_tseries_fromOriginal');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end;

    fprintf('>>> ppi_roi_tseries_fromOriginal_cmd_old2subs：\n');
    fprintf('    被试目录 : %s\n', sess_dir);
    fprintf('    ROI 模式 : %s\n', roi_glob);
    fprintf('    pRF dt   : %s\n', strjoin(dt_prf_list, ', '));
    fprintf('    原始 dt  : %s\n', dt_raw);
    fprintf('    期望 TR  : %d, TR = %.3f s\n\n', NTR_expect, TR);

    cd(sess_dir);

    %% 1. 读取 Gray coords
    coords_file = fullfile(gray_dir, 'coords.mat');
    assert(exist(coords_file, 'file')==2, '找不到 Gray coords: %s', coords_file);
    Scoords = load(coords_file);
    assert(isfield(Scoords,'coords'), 'coords.mat 中没有变量 coords');
    coords_gray = Scoords.coords;   % [3 × nGray] 或更多

    if size(coords_gray,1) > 3
        coords_gray = coords_gray(1:3, :);
    end
    nGray = size(coords_gray, 2);
    fprintf('    [Info] Gray voxel 数 = %d\n\n', nGray);

    %% 2. 启动 mrVista，并获取 dataTYPES
    mrVista 3;
    global VOLUME dataTYPES;
    assert(~isempty(VOLUME), 'mrVista 没有正确启动（VOLUME 为空）');

    all_dt_names = {dataTYPES.name};
    fprintf('    [Info] 当前会话中的 dataTYPES:\n');
    for ii = 1:numel(dataTYPES)
        fprintf('        %2d: %s\n', ii, dataTYPES(ii).name);
    end
    fprintf('\n');



    %% 2.1 确认原始 dataType，并根据手动指定的 Scan 生成 scanList_raw
    dt_idx_raw = find(strcmp(all_dt_names, dt_raw), 1);
    assert(~isempty(dt_idx_raw), 'dataType "%s" 不在 dataTYPES 中', dt_raw);
    nScans_raw = numel(dataTYPES(dt_idx_raw).scanParams);

    allScans_needed = sort(unique([connectScans, unconnectScans]));
    assert(all(allScans_needed >= 1 & allScans_needed <= nScans_raw), ...
        '有手动指定的 Scan 号超出 dataType "%s" 的范围 [1,%d]。', dt_raw, nScans_raw);

    scanList_raw = allScans_needed;   % 只使用这些 Scan

    fprintf('    [Info] 原始 dataType=%s 的 Scan 总数 = %d\n', dt_raw, nScans_raw);
    fprintf('           本脚本实际使用的 Scan 列表 = %s\n\n', mat2str(scanList_raw));

    %% 3. 收集 ROI 列表，并解析左右半球 / base 名
    roi_files_all = dir(roi_glob);
    assert(~isempty(roi_files_all), '在 %s 下找不到 ROI (*.mat)', roi_glob);

    nROI_all = numel(roi_files_all);
    roiNames_all = cell(nROI_all,1);

    for r = 1:nROI_all
        [~, roiNames_all{r}] = fileparts(roi_files_all(r).name);
    end

    fprintf('    [Info] 初始共发现 %d 个 ROI 文件:\n', nROI_all);
    for r = 1:nROI_all
        fprintf('       %-12s (%s)\n', roiNames_all{r}, roi_files_all(r).name);
    end
    fprintf('\n');

    % 仅保留以 L/R/l/r 开头的 ROI（大小写不敏感）
    use_mask    = false(nROI_all,1);
    roiHemi_all = cell(nROI_all,1);
    roiBase_all = cell(nROI_all,1);

    for r = 1:nROI_all
        rn = roiNames_all{r};
        if isempty(rn), continue; end
        firstChar = upper(rn(1));

        if firstChar == 'L'
            use_mask(r)    = true;
            roiHemi_all{r} = 'L';
            roiBase_all{r} = rn(2:end);
        elseif firstChar == 'R'
            use_mask(r)    = true;
            roiHemi_all{r} = 'R';
            roiBase_all{r} = rn(2:end);
        else
            fprintf('    [Info] 跳过非 L/R 开头 ROI: %s\n', rn);
        end
    end

    roi_files = roi_files_all(use_mask);
    roiNames  = roiNames_all(use_mask);
    roiHemi   = roiHemi_all(use_mask);
    roiBase   = roiBase_all(use_mask);

    nROI = numel(roiNames);
    assert(nROI > 0, '筛选 L/R 开头后，没有剩余 ROI。');

    fprintf('    [Info] 筛选后保留 %d 个 L/R ROI:\n', nROI);
    for r = 1:nROI
        fprintf('       %-12s (hemi=%s, base=%s)\n', roiNames{r}, roiHemi{r}, roiBase{r});
    end
    fprintf('\n');

    % 所有出现的 base ROI
    baseROI_all = unique(roiBase(~cellfun(@isempty, roiBase)));
    nBase = numel(baseROI_all);

    idxL = nan(nBase,1);
    idxR = nan(nBase,1);
    for b = 1:nBase
        bn = baseROI_all{b};
        idxL(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'L'), 1);
        idxR(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'R'), 1);
        if isnan(idxL(b)) || isnan(idxR(b))
            fprintf('    [Warn] base ROI %s 在某个半球缺失 (L idx=%g, R idx=%g)。\n', ...
                bn, idxL(b), idxR(b));
        end
    end

    %% 4. 构建 run→条件(connect/unconnect) 信息（用手动 Scan 列表）
    nScan = numel(scanList_raw);

    scanInfo = struct();
    scanInfo.scanNums   = scanList_raw(:);        % 实际 Scan 编号
    scanInfo.sessionIdx = ones(nScan,1);          % 旧数据已 merge 到 ses-01，这里统一记为 1
    scanInfo.cond       = cell(nScan,1);          % 'connect' or 'unconnect'

    for i = 1:nScan
        sNum = scanList_raw(i);
        if ismember(sNum, connectScans)
            scanInfo.cond{i} = 'connect';
        elseif ismember(sNum, unconnectScans)
            scanInfo.cond{i} = 'unconnect';
        else
            error('Scan%d 既不在 connectScans 也不在 unconnectScans 中，请检查配置。', sNum);
        end
    end

    fprintf('    [Info] run 条件标记:\n');
    for i = 1:nScan
        fprintf('        Scan%2d -> ses%d, %s\n', ...
            scanInfo.scanNums(i), scanInfo.sessionIdx(i), scanInfo.cond{i});
    end
    fprintf('\n');

    run_idx_connect   = find(strcmp(scanInfo.cond, 'connect'));
    run_idx_unconnect = find(strcmp(scanInfo.cond, 'unconnect'));

    %% 5. 在两套 pRF dataType 下，分别计算 ROI 的阈值掩膜和 Gray 索引

    roiPref        = struct();
    roiAmp         = struct();
    roiVe          = struct();
    roiMask_prf    = struct();
    roiNvox_total  = struct();
    roiNvox_keep   = struct();
    roiGrayIdx_all = struct();
    roiGrayIdx_keep= struct();
    model_lin      = struct();

    for c = 1:numel(dt_prf_list)
        dt_prf = dt_prf_list{c};

        % 根据 dataType 名推断 condLabel
        % Averages_single_layer_01 -> 'unconnect'
        % Averages_connected_layer_01 (含 connect) -> 'connect'
        if contains(dt_prf, 'single')
            condLabel = 'unconnect';
        elseif contains(dt_prf, 'connect')
            condLabel = 'connect';
        else
            error('无法从 dataType "%s" 推断 condLabel（connect/unconnect/single）', dt_prf);
        end

        fprintf('>>> [pRF] dataType = %s (cond = %s)\n', dt_prf, condLabel);

        dt_idx_prf = find(strcmp(all_dt_names, dt_prf), 1);
        assert(~isempty(dt_idx_prf), 'dataType "%s" 不在 dataTYPES 中', dt_prf);

        VOLUME{1} = viewSet(VOLUME{1}, 'curdt',   dt_idx_prf);
        VOLUME{1} = viewSet(VOLUME{1}, 'curscan', scanNum_prf);

        % 检查该 pRF dt 的 TSeries 帧数
        tseries_prf_file = fullfile(gray_dir, dt_prf, 'TSeries', ...
            sprintf('Scan%d', scanNum_prf), 'tSeries1.mat');
        assert(exist(tseries_prf_file,'file')==2, '找不到 %s', tseries_prf_file);
        Sprf = load(tseries_prf_file);
        assert(isfield(Sprf,'tSeries'), 'pRF tSeries 文件中没有变量 tSeries');
        [nFrames_prf, nVox_prf] = size(Sprf.tSeries);
        fprintf('    [Info] pRF dt=%s, scan=%d：帧数=%d, voxel=%d\n', ...
            dt_prf, scanNum_prf, nFrames_prf, nVox_prf);
        assert(nFrames_prf == NTR_expect, ...
            'pRF dt=%s 帧数(%d) != 期望(%d)，请检查 NTR_expect 或该 pRF 会话', ...
            dt_prf, nFrames_prf, NTR_expect);

        % 找该 dt 下的 linearized 模型
        model_dir = fullfile(gray_dir, dt_prf);
        lin_list  = dir(fullfile(model_dir, '*gFit-gFit_linearized.mat'));
        assert(~isempty(lin_list),  '在 %s 下找不到 *gFit-gFit_linearized.mat', model_dir);
        lin_file  = fullfile(model_dir, lin_list(1).name);
        model_lin.(condLabel) = lin_list(1).name;

        fprintf('    [Info] pRF linearized 模型文件：%s\n', lin_list(1).name);

        % 对每个 ROI，在该 dt 下获取 pref/amp/ve，并映射到 Gray coords
        for r = 1:nROI
            rn       = roiNames{r};
            roi_path = fullfile(roi_files(r).folder, roi_files(r).name);

            fprintf('    [pRF ROI] (%s) %-8s : ', condLabel, rn);

            [VOLUME{1}, ok] = loadROI(VOLUME{1}, roi_path, 1, 'b', 1); %#ok<ASGLU>
            roiNum = viewGet(VOLUME{1}, 'numROIs');

            VOLUME{1} = rmSelect(VOLUME{1}, 1, lin_file);
            VOLUME{1} = rmLoadDefault(VOLUME{1}, 0);

            pref = getCurDataROI(VOLUME{1}, 'map');
            amp  = getCurDataROI(VOLUME{1}, 'amp');
            ve   = getCurDataROI(VOLUME{1}, 'co');

            pref_r = pref(:);
            amp_r  = amp(:);
            ve_r   = ve(:);

            n_vox_roi = numel(pref_r);

            coords_roi = viewGet(VOLUME{1}, 'ROICoords', roiNum);  % [3 × nVoxROI]
            if size(coords_roi,1) > 3
                coords_roi = coords_roi(1:3, :);
            end
            n_vox_coord = size(coords_roi, 2);

            if n_vox_coord ~= n_vox_roi
                warning('ROI %s (%s): pref 长度(%d) 与 coords voxel 数(%d) 不一致，截取较小值。', ...
                    rn, condLabel, n_vox_roi, n_vox_coord);
                n_min      = min(n_vox_roi, n_vox_coord);
                pref_r     = pref_r(1:n_min);
                amp_r      = amp_r(1:n_min);
                ve_r       = ve_r(1:n_min);
                coords_roi = coords_roi(:,1:n_min);
                n_vox_roi  = n_min;
            end

            [~, idx_gray] = ismember(coords_roi', coords_gray', 'rows');
            if any(idx_gray == 0)
                warning('ROI %s (%s): 有 %d 个 voxel 没有在 Gray coords 中找到，将丢弃。', ...
                    rn, condLabel, sum(idx_gray==0));
                keep_ok        = idx_gray > 0;
                idx_gray       = idx_gray(keep_ok);
                pref_r         = pref_r(keep_ok);
                amp_r          = amp_r(keep_ok);
                ve_r           = ve_r(keep_ok);
                coords_roi     = coords_roi(:, keep_ok);
                n_vox_roi      = numel(pref_r);
            end

            mask_roi = (ve_r > 0.3) & (pref_r >= 1.05) & (pref_r <= 64);
            n_keep   = sum(mask_roi);

            if n_keep == 0
                warning('ROI %s (%s): 无 voxel 通过阈值，后续将用全部 voxel 做平均。', rn, condLabel);
                idx_keep_gray = idx_gray;
            else
                idx_keep_gray = idx_gray(mask_roi);
            end

            % 按 condLabel 和 ROI 保存
            roiPref.(condLabel).(rn)        = pref_r;
            roiAmp.(condLabel).(rn)         = amp_r;
            roiVe.(condLabel).(rn)          = ve_r;

            roiMask_prf.(condLabel).(rn)    = mask_roi;
            roiNvox_total.(condLabel).(rn)  = n_vox_roi;
            roiNvox_keep.(condLabel).(rn)   = n_keep;

            roiGrayIdx_all.(condLabel).(rn) = idx_gray(:);
            roiGrayIdx_keep.(condLabel).(rn)= idx_keep_gray(:);

            fprintf('voxel 总数=%d, 通过阈值=%d\n', n_vox_roi, n_keep);
        end

        fprintf('\n');
    end

    %% 6. 在 Gray/Original 中，对每个条件分别提取 ROI 原始时序（平均 + 体素级）

    condList = {'connect', 'unconnect'};

    roiSeed   = struct();   % roiSeed.(condLabel).(ROI): [NTR × nRuns_cond]
    roiTS_vox = struct();   % roiTS_vox.(condLabel).(ROI){k}: [NTR × nKeep]

    for c = 1:numel(condList)
        condLabel = condList{c};

        if strcmp(condLabel, 'connect')
            run_idx = run_idx_connect;
        else
            run_idx = run_idx_unconnect;
        end

        nRun_cond = numel(run_idx);
        if nRun_cond == 0
            warning('没有任何 run 标注为 %s，跳过该条件的时序提取。', condLabel);
            continue;
        end

        fprintf('>>> [原始时序] 条件 = %s，nRuns = %d\n', condLabel, nRun_cond);

        % 为该条件下的每个 ROI 分配存储
        for r = 1:nROI
            rn = roiNames{r};
            roiSeed.(condLabel).(rn)   = nan(NTR_expect, nRun_cond);
            roiTS_vox.(condLabel).(rn) = cell(1, nRun_cond);
        end

        % 遍历该条件下的所有 run（注意这里的 k 是该条件内的索引）
        for k = 1:nRun_cond
            i       = run_idx(k);             % 在本脚本使用的所有 run 中的索引
            scanNum = scanList_raw(i);

            fprintf('    条件=%s, run %d/%d: Scan%d (ses%d)\n', ...
                condLabel, k, nRun_cond, scanNum, scanInfo.sessionIdx(i));

            tseries_raw_file = fullfile(gray_dir, dt_raw, 'TSeries', ...
                sprintf('Scan%d', scanNum), 'tSeries1.mat');
            assert(exist(tseries_raw_file,'file')==2, '找不到 %s', tseries_raw_file);

            Sraw = load(tseries_raw_file);
            assert(isfield(Sraw,'tSeries'), '原始 tSeries 文件中没有变量 tSeries');
            ts_all = Sraw.tSeries;   % [NTR × nGray]

            [nT, nVox] = size(ts_all);
            fprintf('        [Info] Scan%d: 帧数=%d, voxel=%d\n', scanNum, nT, nVox);

            assert(nT == NTR_expect, ...
                'Scan%d: 帧数(%d) != 期望(%d)，请检查 NSS TR 是否一致。', ...
                scanNum, nT, NTR_expect);
            assert(nVox == nGray, ...
                'Scan%d: voxel 数(%d) != Gray coords voxel 数(%d)，请检查数据对齐。', ...
                scanNum, nVox, nGray);

            for r = 1:nROI
                rn = roiNames{r};

                % 对于当前条件 condLabel，使用该条件下 pRF 阈值得到的 voxel 索引
                idx_keep_gray = roiGrayIdx_keep.(condLabel).(rn);   % [nKeep × 1]
                ts_roi_vox    = ts_all(:, idx_keep_gray);           % [NTR × nKeep]

                seed = mean(ts_roi_vox, 2, 'omitnan');              % [NTR × 1]

                roiSeed.(condLabel).(rn)(:, k) = seed;
                roiTS_vox.(condLabel).(rn){k}  = ts_roi_vox;
            end
        end

        fprintf('\n');
    end

    %% 7. 保存结果

    meta = struct();
    meta.subj        = subj_label;       % 使用逻辑被试名（sub-02 / sub-08）
    meta.ses         = ses;
    meta.dirID       = dirID;            % 实际目录名
    meta.TR          = TR;
    meta.NTR         = NTR_expect;
    meta.dt_prf_list = dt_prf_list;
    meta.dt_raw      = dt_raw;
    meta.scanList    = scanList_raw;
    meta.scanInfo    = scanInfo;
    meta.runIdx_connect   = run_idx_connect;
    meta.runIdx_unconnect = run_idx_unconnect;

    meta.roiNames        = roiNames;     % 带 L/R 的 ROI 名
    meta.roiHemi         = roiHemi;      % 'L'/'R'
    meta.roiBaseNames    = roiBase;      % base ROI 名
    meta.baseROI_all     = baseROI_all;  % 所有出现过的 base ROI
    meta.idxL            = idxL;         % 每个 base 的 L 索引（在 roiNames 中，可为 NaN）
    meta.idxR            = idxR;         % 每个 base 的 R 索引（在 roiNames 中，可为 NaN）

    meta.roiNvox     = roiNvox_total;     % 按 condLabel 再分 ROI
    meta.roiNkeep    = roiNvox_keep;      % 按 condLabel 再分 ROI
    meta.model_lin   = model_lin;         % meta.model_lin.connect / .unconnect
    meta.note        = [...
        'roiSeed.(cond).(ROI): [NTR × nRuns_cond]，每列为该条件下某个 run 的平均时序（PPI 生理量）； ', ...
        'roiTS_vox.(cond).(ROI){k}: [NTR × nKeep]，该条件第 k 个 run 的体素级原始时序； ', ...
        'roiGrayIdx_keep.(cond).(ROI): 在 Gray coords 空间中，通过 pRF 阈值的 voxel 索引； ', ...
        'cond="connect"/"unconnect" 是根据手动指定的 Scan 列表（connectScans/unconnectScans）赋值的； ', ...
        'ROI 只保留以 L/R/l/r 开头者，并在 meta.roiHemi / roiBaseNames / baseROI_all / idxL / idxR 中给出左右半球信息。' ...
    ];

    out_name = fullfile(out_dir, sprintf('%s_%s_Original_roiTS_multiRun_dualPRF.mat', subj_label, ses));
    save(out_name, 'roiSeed', 'roiTS_vox', ...
        'roiPref', 'roiAmp', 'roiVe', ...
        'roiMask_prf', 'roiGrayIdx_all', 'roiGrayIdx_keep', ...
        'meta', '-v7.3');

    fprintf('>>> 已保存 ROI 多 run 原始时序 (connect/unconnect 并行) 到: %s\n', out_name);
    fprintf('=== ppi_roi_tseries_fromOriginal_cmd_old2subs 完成: %s ===\n\n', subj_label);
end
