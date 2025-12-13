%% prf_roi_tseries_cmd.m
% 纯命令行 pRF ROI 时序提取
% 依赖：vistasoft/mrVista 在 path 中，且该被试听过 pRF 并在 Gray 里有
%       Averages_connect_layer_01 / Averages_unconnect_layer_01 等 dataType。
%
% 功能：
%   - 对指定会话（subj+ses）：
%       * 读取 Gray/ROIs 中的 LN* 数量 ROI；
%       * 对每个指定 dataType（连线/不连线）：
%           + 用 *gFit-gFit_linearized.mat 读取 map/amp/co（pref/width/ve）；
%           + 用 *gFit-gFit.mat + rmPlotGUI 拿 ROI 内所有体素的 tSeries；
%           + 用 ve>0.3 且 1.05<=pref<=64 做阈值；
%           + 输出两个版本：
%               - roiSeed.(ROI)：阈值后体素平均，得到 1×NTR 的 ROI seed
%               - roiTS_vox.(ROI)：阈值后每个体素的原始 time series，NTR×nKeep
%       * 把结果保存到 results_prf_roi_tseries_cmd 下的 mat 文件。

clear; clc;

%% 0. 基本配置（根据需要修改）
subj      = 'sub-03';
ses       = 'ses-01';
root_prf  = '/home/xue/data/prf_result';

% 会话目录（pRF 会话）
sess_dir  = fullfile(root_prf, subj, ses);

% LN* 数量 ROI
roi_glob  = fullfile(sess_dir, 'Gray', 'ROIs', '*N*.mat');

% 想要提取的 dataType，可以只留一个，也可以两个都提
dataTypes_target = { ...
    'Averages_unconnect_layer_01', ...  % 不连线
    'Averages_connect_layer_01'        % 连线
};

scanNum    = 1;      % pRF 的 scan 号（你现在就是 Scan1）
NTR_expect = 176;    % pRF run 的帧数

% 输出目录
out_dir = fullfile(sess_dir, 'results_prf_roi_tseries_cmd');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fprintf('>>> prf_roi_tseries_cmd：\n');
fprintf('    会话目录 : %s\n', sess_dir);
fprintf('    ROI 模式 : %s\n', roi_glob);
fprintf('    使用 scan: %d\n', scanNum);
fprintf('    期望 NTR : %d\n\n', NTR_expect);

cd(sess_dir);

%% 1. 启动 mrVista，并查看 dataTYPES
mrVista 3;   % 会弹 GUI 窗口，但你不用动它

global VOLUME dataTYPES;
assert(~isempty(VOLUME), 'mrVista 没有正确启动（VOLUME 为空）');

fprintf('    [Info] 会话中的 dataTYPES:\n');
for ii = 1:numel(dataTYPES)
    fprintf('        %d: %s\n', ii, dataTYPES(ii).name);
end
fprintf('\n');

%% 2. 收集 ROI 列表（LN* 数量 ROI）
roi_files = dir(roi_glob);
assert(~isempty(roi_files), '在 %s 下找不到 ROI (*.mat)', roi_glob);
nROI = numel(roi_files);

roiNames = cell(nROI,1);
for r = 1:nROI
    [~, roiNames{r}] = fileparts(roi_files(r).name);
end

fprintf('    [Info] 共发现 %d 个 ROI:\n', nROI);
for r = 1:nROI
    fprintf('       %-8s (%s)\n', roiNames{r}, roi_files(r).name);
end
fprintf('\n');

%% 3. 对每一个 dataType 做 ROI 种子提取
for dt_i = 1:numel(dataTypes_target)
    dtName = dataTypes_target{dt_i};

    fprintf('>>> DataType: %s\n', dtName);

    % 每一轮循环都从当前 dataTYPES 中定位 dataType index（更稳）
    global dataTYPES VOLUME;
    all_dt_names = {dataTYPES.name};
    dt_idx = find(strcmp(all_dt_names, dtName), 1);
    assert(~isempty(dt_idx), 'dataType "%s" 不在 dataTYPES 中', dtName);

    %% 3.1 切换到对应 dataType & scan，并检查 TSeries 帧数
    VOLUME{1} = viewSet(VOLUME{1}, 'curdt',   dt_idx);
    VOLUME{1} = viewSet(VOLUME{1}, 'curscan', scanNum);

    tseries_file = fullfile(sess_dir, 'Gray', dataTYPES(dt_idx).name, 'TSeries', ...
                            sprintf('Scan%d', scanNum), 'tSeries1.mat');
    assert(exist(tseries_file,'file')==2, '找不到 %s', tseries_file);
    S = load(tseries_file);
    assert(isfield(S,'tSeries'), 'tSeries 文件中没有变量 tSeries');
    [nFrames, ~] = size(S.tSeries);
    fprintf('    [Info] dataType=%s, scan=%d 的帧数 = %d\n', dataTYPES(dt_idx).name, scanNum, nFrames);
    assert(nFrames == NTR_expect, ...
           '帧数(%d) != 期望(%d)，请检查 pRF 会话与 NTR_expect 设置', ...
           nFrames, NTR_expect);

    %% 3.2 找到模型文件：linearized（取 map/amp/co）和 gFit-gFit（取 tSeries/pred）
    model_dir = fullfile(sess_dir, 'Gray', dataTYPES(dt_idx).name);

    lin_list  = dir(fullfile(model_dir, '*gFit-gFit_linearized.mat'));
    gfit_list = dir(fullfile(model_dir, '*gFit-gFit.mat'));

    assert(~isempty(lin_list),  '在 %s 下找不到 *gFit-gFit_linearized.mat', model_dir);
    assert(~isempty(gfit_list), '在 %s 下找不到 *gFit-gFit.mat',           model_dir);

    lin_file  = fullfile(model_dir, lin_list(1).name);
    gfit_file = fullfile(model_dir, gfit_list(1).name);

    fprintf('    [Info] 参数图模型 (map/amp/co) ：%s\n', lin_list(1).name);
    fprintf('    [Info] 时序模型   (tSeries)    ：%s\n', gfit_list(1).name);

    %% 3.3 逐个 ROI 提取：先用 linearized 取 map/amp/co，再用 gFit-gFit 取 tSeries
    roiSeed    = struct();   % 每个 ROI 的平均 time series（NTR×1）
    roiTS_vox  = struct();   % 每个 ROI 的 voxel-wise time series（NTR×nKeep）
    roiMask    = struct();   % 每个 ROI 的逻辑 mask（nVox×1），哪些 voxel 被保留
    roiNvox    = struct();   % 每个 ROI 的总体素数
    roiNkeep   = struct();   % 每个 ROI 的通过阈值体素数
    roiPref    = struct();   % 每个 ROI 的 pref map（nVox×1）
    roiAmp     = struct();   % 每个 ROI 的 amp map（nVox×1）
    roiVe      = struct();   % 每个 ROI 的 ve map （nVox×1）

    for r = 1:nROI
        rn       = roiNames{r};
        roi_path = fullfile(roi_files(r).folder, roi_files(r).name);

        fprintf('    [ROI] %-8s : ', rn);

        %--- 3.3.1 加载 ROI（只需加载一次，后面切模型即可重用）
        [VOLUME{1}, ok] = loadROI(VOLUME{1}, roi_path, 1, 'b', 1); %#ok<ASGLU>

        %--- 3.3.2 选 linearized 模型，读取 map/amp/co（pref/width/ve）
        VOLUME{1} = rmSelect(VOLUME{1}, 1, lin_file);
        VOLUME{1} = rmLoadDefault(VOLUME{1}, 0);

        pref = getCurDataROI(VOLUME{1}, 'map');  % 偏好 numerosity
        amp  = getCurDataROI(VOLUME{1}, 'amp');  % tuning width
        ve   = getCurDataROI(VOLUME{1}, 'co');   % variance explained

        pref_r = pref(:);
        amp_r  = amp(:);
        ve_r   = ve(:);

        %--- 3.3.3 换成 gFit-gFit 模型，用 rmPlotGUI 拿 tSeries
        VOLUME{1} = rmSelect(VOLUME{1}, 1, gfit_file);
        VOLUME{1} = rmLoadDefault(VOLUME{1}, 0);

        M      = rmPlotGUI(VOLUME{1});   % 命令行调用，只用 M.tSeries
        ts_all = M.tSeries;              % [NTR x nVox]

        [nT, nVox] = size(ts_all);
        if nT ~= NTR_expect
            warning('tSeries 帧数 %d != 期望 %d，跳过该 ROI', nT, NTR_expect);
            continue;
        end

        %--- 3.3.4 体素筛选：ve>0.3 & 1.05<=pref<=64（完全按你原来的逻辑）
        mask = (ve_r > 0.3) & (pref_r >= 1.05) & (pref_r <= 64);

        n_vox_total = nVox;
        n_vox_keep  = sum(mask);

        if n_vox_keep == 0
            warning('无体素通过阈值，退化为对全部体素做平均');
            seed      = mean(ts_all, 2, 'omitnan');  % NTR×1
            ts_keep   = ts_all;                      % NTR×nVox（全部保留）
            mask_keep = true(nVox,1);
        else
            seed      = mean(ts_all(:, mask), 2, 'omitnan');  % NTR×1
            ts_keep   = ts_all(:, mask);                      % NTR×nKeep
            mask_keep = mask;                                 % nVox×1
        end

        %--- 3.3.5 保存到结构体
        roiSeed.(rn)   = seed(:);      % NTR×1
        roiTS_vox.(rn) = ts_keep;      % NTR×nKeep
        roiMask.(rn)   = mask_keep;    % nVox×1

        roiNvox.(rn)   = n_vox_total;
        roiNkeep.(rn)  = n_vox_keep;

        roiPref.(rn)   = pref_r;       % nVox×1
        roiAmp.(rn)    = amp_r;        % nVox×1
        roiVe.(rn)     = ve_r;         % nVox×1

        fprintf('体素总数=%d, 通过阈值=%d, seed 长度=%d, voxel TS 尺寸=[%d × %d]\n', ...
                n_vox_total, n_vox_keep, numel(seed), size(ts_keep,1), size(ts_keep,2));
    end

    %% 3.4 保存结果
    out_name = fullfile(out_dir, sprintf('%s_%s_roiSeed.mat', subj, dataTYPES(dt_idx).name));
    meta = struct();
    meta.subj      = subj;
    meta.ses       = ses;
    meta.dtName    = dataTYPES(dt_idx).name;
    meta.scan      = scanNum;
    meta.NTR       = NTR_expect;
    meta.roiNames  = roiNames;
    meta.roiNvox   = roiNvox;
    meta.roiNkeep  = roiNkeep;
    meta.model_lin = lin_list(1).name;
    meta.model_gfit= gfit_list(1).name;
    meta.note      = ['roiSeed: NTR×1 平均时序; ' ...
                      'roiTS_vox: NTR×nKeep 体素级时序(按 roiMask 筛选)'];

    save(out_name, 'roiSeed', 'roiTS_vox', 'roiMask', ...
                    'roiPref', 'roiAmp', 'roiVe', 'meta');

    fprintf('>>> 已保存 ROI 种子 + voxel 时序到: %s\n\n', out_name);
end

fprintf('=== prf_roi_tseries_cmd 完成 ===\n');


