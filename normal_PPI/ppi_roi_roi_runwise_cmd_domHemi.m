%% ppi_roi_roi_runwise_cmd_domHemi.m
% 从 Gray/Original 的 run-wise ROI 时序 + fMRIPrep confounds 做 ROI-ROI PPI
%   - 每个 run 单独建模（run-wise）
%   - confounds 使用 fMRIPrep 输出，按自定义规则重算（motion, FD, aCompCor WM/CSF）
%   - 本版本：增加左右半球固定效应合并 + 秩感知 QR 求解 + confounds 对 Z 的正交化：
%       * 在 ROI(含 L/R) 空间中做 run-wise PPI
%       * 设计矩阵 X = [Z, N_ortho]，其中 Z = [1, psych, phys, ppi]，
%         N_ortho 为 confounds 对 Z 正交化后的部分（不改变 PPI 本身）
%       * 用 safe_glm_rankaware 在秩 r 子空间求解 GLM，自动处理秩亏
%       * 在被试内部先对每个半球单独跨 run 平均，得到 β_L、β_R
%       * 再对 β_L、β_R 简单平均，得到 base ROI–base ROI 级别的 bilateral β

clear; clc;
fprintf('\n================ ROI-ROI PPI (run-wise, from Original, bilateral FFX) ================\n');

%% 0. 基本配置（按需要修改）
cfg = struct();

% 被试与路径
cfg.subj       = 'sub-09';
cfg.prf_ses    = 'ses-01';  % prf_result 下面的会话名（Gray/Original 在这里）
cfg.root_prf   = '/home/xue/data/prf_result';      % 有 Gray/Original 的地方
cfg.root_bids  = '/home/xue/data/BIDS_prep';       % fMRIPrep 输出的 BIDS_prep 根
cfg.t1_ses     = 'ses-03';  % anat / WM/CSF 概率图所在的 session

% ROI multi-run 原始时序 .mat（前一步脚本生成）
cfg.roi_ts_dir = 'results_ppi_roi_tseries_fromOriginal';
cfg.roi_ts_file= sprintf('%s_%s_Original_roiTS_multiRun_dualPRF.mat', ...
                         cfg.subj, cfg.prf_ses);

% events（所有 run 共用）
cfg.events_tsv = '/media/xue/new_B1/glm/events_output/events_bids.tsv';

% 时间信息
cfg.TR         = 1.95;
cfg.n_dummy    = 6;         % 原始 BOLD / confounds 前面的 NSS 数量
cfg.NTR_target = 176;       % 最终每个 run 的长度

% 心理条件（trial_type）名：1 = numerosity, 2 = baseline
cfg.cond_names = {'numerosity', 'baseline'};

% 是否在 GLM 中加入混杂项（true = 加入；false = 暂不加入，只构建）
cfg.use_confounds = true;

% fMRIPrep confounds 相关参数
cfg.fd_thr        = 0.3;    % FD 阈值（mm）
cfg.acompcor_thr  = 0.9;    % （旧版 PCA 用的阈值，现在主要是记录）
cfg.nComp_WM      = 5;      % WM aCompCor 组件数
cfg.nComp_CSF     = 5;      % CSF aCompCor 组件数

% 输出目录
cfg.out_dir = fullfile(cfg.root_bids, 'derivatives', 'ppi_roi_runwise', cfg.subj);
if ~exist(cfg.out_dir, 'dir')
    fprintf('  [Info] 创建输出目录: %s\n', cfg.out_dir);
    mkdir(cfg.out_dir);
end

fprintf('  [Config] subj=%s\n', cfg.subj);
fprintf('  [Config] TR=%.3f, NTR_target=%d, n_dummy=%d\n', ...
    cfg.TR, cfg.NTR_target, cfg.n_dummy);
fprintf('  [Config] cond_names: %s\n', strjoin(cfg.cond_names, ', '));
fprintf('  [Config] use_confounds = %d\n', cfg.use_confounds);

%% 1. 载入 ROI 时序 (run-wise, Original) 并解析 L/R + base ROI 信息
fprintf('\n[Step 1] 载入 run-wise ROI 时序 (Original)...\n');

sess_dir  = fullfile(cfg.root_prf, cfg.subj, cfg.prf_ses);
roi_ts_mat= fullfile(sess_dir, cfg.roi_ts_dir, cfg.roi_ts_file);
assert(exist(roi_ts_mat, 'file')==2, '找不到 ROI ts 文件: %s', roi_ts_mat);

S = load(roi_ts_mat);
req_fields = {'roiSeed','meta'};
missing = setdiff(req_fields, fieldnames(S));
assert(isempty(missing), 'ROI ts .mat 中缺少字段: %s', strjoin(missing, ', '));

roiSeed = S.roiSeed;
meta    = S.meta;

fprintf('    [Info] meta.NTR = %d, cfg.NTR_target = %d\n', meta.NTR, cfg.NTR_target);
assert(meta.NTR == cfg.NTR_target, ...
    'meta.NTR (%d) 与 cfg.NTR_target (%d) 不一致', meta.NTR, cfg.NTR_target);

roiNames   = meta.roiNames(:);
nROI       = numel(roiNames);
scanInfo   = meta.scanInfo;
nScan      = numel(scanInfo.scanNums);
nRuns_ses1 = meta.nRuns_ses1;
nRuns_ses2 = meta.nRuns_ses2;

runIdx_connect   = meta.runIdx_connect;
runIdx_unconnect = meta.runIdx_unconnect;

fprintf('    [Info] ROI 数 = %d\n', nROI);
fprintf('    [Info] nScan 总共 = %d, ses1=%d, ses2=%d\n', ...
        nScan, nRuns_ses1, nRuns_ses2);

% 解析 ROI 的 hemis/base 名
roiHemi = cell(nROI,1);
roiBase = cell(nROI,1);
for i = 1:nROI
    rn = roiNames{i};
    if startsWith(rn,'L')
        roiHemi{i} = 'L';
        roiBase{i} = rn(2:end);
    elseif startsWith(rn,'R')
        roiHemi{i} = 'R';
        roiBase{i} = rn(2:end);
    else
        roiHemi{i} = '';
        roiBase{i} = rn;
    end
end

baseROI_all = unique(roiBase(~cellfun(@isempty, roiBase)));
% 仅保留同时存在 L/R 的 base ROI，作为 bilateral 网络的节点
hasL = false(size(baseROI_all));
hasR = false(size(baseROI_all));
for b = 1:numel(baseROI_all)
    bn   = baseROI_all{b};
    hasL(b) = any(strcmp(roiBase,bn) & strcmp(roiHemi,'L'));
    hasR(b) = any(strcmp(roiBase,bn) & strcmp(roiHemi,'R'));
end
baseROI_list = baseROI_all(hasL & hasR);
nBase = numel(baseROI_list);

idxL = nan(nBase,1);
idxR = nan(nBase,1);
for b = 1:nBase
    bn = baseROI_list{b};
    idxL(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'L'), 1);
    idxR(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'R'), 1);
end

fprintf('    [Info] bilateral base ROI 数 = %d\n', nBase);
for b = 1:nBase
    fprintf('       base=%-8s : L-idx=%d (%s), R-idx=%d (%s)\n', ...
        baseROI_list{b}, ...
        idxL(b), roiNames{idxL(b)}, ...
        idxR(b), roiNames{idxR(b)});
end

%% 1.1 构建每个 scan 的 session / 条件 / runInSes 信息（为 BIDS 对齐做准备）
runMeta = struct('scanIdx', cell(1,nScan), ...
                 'scanNum', [], 'sessionIdx', [], 'cond', [], ...
                 'runInSes', [], 'sesLabel', [], 'runLabel', []);

for i = 1:nScan
    sNum = scanInfo.scanNums(i);          % 通常等于 i
    sesIdx = scanInfo.sessionIdx(i);      % 1 or 2
    cond_i = scanInfo.cond{i};            % 'connect' / 'unconnect'

    if sesIdx == 1
        runInSes = sNum;                  % ses1 内部的 run 号
    else
        runInSes = sNum - nRuns_ses1;     % ses2 内部的 run 号
    end

    runMeta(i).scanIdx   = i;
    runMeta(i).scanNum   = sNum;
    runMeta(i).sessionIdx= sesIdx;
    runMeta(i).cond      = cond_i;
    runMeta(i).runInSes  = runInSes;
    runMeta(i).sesLabel  = sprintf('ses-%02d', sesIdx);
    runMeta(i).runLabel  = sprintf('run-%03d', runInSes);

    fprintf('    Scan%2d -> %s, %s, %s\n', ...
        sNum, runMeta(i).sesLabel, runMeta(i).runLabel, cond_i);
end

%% 2. 构造心理回归量 Psych_diff（numerosity vs baseline），每个 run 共用
fprintf('\n[Step 2] 从 events_bids.tsv 构造 Psych_diff...\n');
[psych_mat, psych_names] = build_psych(cfg);
NTR_psych = size(psych_mat,1);
assert(NTR_psych == cfg.NTR_target, ...
    'Psych_diff 长度(%d) != cfg.NTR_target(%d)', NTR_psych, cfg.NTR_target);

%% 3. 为每个 run 构建 fMRIPrep confounds（motion + FD + aCompCor）
fprintf('\n[Step 3] 为每个 run 构建 fMRIPrep 混杂项...\n');

% WM / CSF 概率图路径（T1 space）
wm_nii  = fullfile(cfg.root_bids, cfg.subj, cfg.t1_ses, 'anat', ...
                   sprintf('%s_%s_label-WM_probseg.nii',  cfg.subj, cfg.t1_ses));
csf_nii = fullfile(cfg.root_bids, cfg.subj, cfg.t1_ses, 'anat', ...
                   sprintf('%s_%s_label-CSF_probseg.nii', cfg.subj, cfg.t1_ses));

if ~exist(wm_nii, 'file')
    wm_nii_gz = [wm_nii '.gz'];
    assert(exist(wm_nii_gz,'file')==2, 'WM_probseg 文件不存在: %s / %s', wm_nii, wm_nii_gz);
    wm_nii = wm_nii_gz;
end

if ~exist(csf_nii, 'file')
    csf_nii_gz = [csf_nii '.gz'];
    assert(exist(csf_nii_gz,'file')==2, 'CSF_probseg 文件不存在: %s / %s', csf_nii, csf_nii_gz);
    csf_nii = csf_nii_gz;
end

% 针对 connect / unconnect 分别存 confounds（顺序要和 roiSeed 一致）
condLevels = {'connect','unconnect'};
confounds = struct();

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    if strcmp(condLabel,'connect')
        run_idx_cond = runIdx_connect;
    else
        run_idx_cond = runIdx_unconnect;
    end

    nRun_cond = numel(run_idx_cond);
    fprintf('\n    [Cond=%s] nRun_cond = %d\n', condLabel, nRun_cond);

    confounds.(condLabel).X      = cell(1, nRun_cond);  % 每个 run 一个 [NTR×P] 矩阵
    confounds.(condLabel).names  = cell(1, nRun_cond);  % 每列名
    confounds.(condLabel).meta   = cell(1, nRun_cond);  % 文件名等信息

    for k = 1:nRun_cond
        i = run_idx_cond(k);      % 全部 run 中的索引
        rm = runMeta(i);

        sesDir  = fullfile(cfg.root_bids, cfg.subj, rm.sesLabel, 'func');
        baseStub= sprintf('%s_%s_task-social_%s', ...
                          cfg.subj, rm.sesLabel, rm.runLabel);

        % confounds tsv
        conf_tsv = fullfile(sesDir, [baseStub '_desc-confounds_timeseries.tsv']);
        assert(exist(conf_tsv,'file')==2, '找不到 confounds: %s', conf_tsv);

        % BOLD (T1w space) —— 不再用于 aCompCor，仅记录
        bold_nii1 = fullfile(sesDir, [baseStub '_space-T1w_desc-preproc_bold.nii']);
        bold_nii2 = [bold_nii1 '.gz'];
        if exist(bold_nii1,'file')
            bold_nii = bold_nii1;
        elseif exist(bold_nii2,'file')
            bold_nii = bold_nii2;
        else
            error('找不到 BOLD 文件: %s / %s', bold_nii1, bold_nii2);
        end

        fprintf('        Run %d/%d: %s, %s\n', k, nRun_cond, rm.sesLabel, rm.runLabel);
        fprintf('            confounds: %s\n', conf_tsv);
        fprintf('            BOLD     : %s\n', bold_nii);

        % 调用子函数构建一个 run 的 confounds
        [Xc, cnames, cmeta] = build_confounds_one_run(conf_tsv, bold_nii, ...
            wm_nii, csf_nii, cfg);

        confounds.(condLabel).X{k}     = Xc;
        confounds.(condLabel).names{k} = cnames;
        confounds.(condLabel).meta{k}  = cmeta;
    end
end

%% 4. run-wise ROI-ROI PPI：每个 run 单独建模（ROI 仍是带 L/R 的 ROI 列表）
fprintf('\n[Step 4] run-wise ROI-ROI PPI GLM...\n');

% ---------- DCT 高通滤波矩阵（对所有 run 通用） ----------
N  = cfg.NTR_target;
TR = cfg.TR;

% 高通截止周期（秒），可按需要修改；这里 256 s
hp_cutoff = 256;

% 参考 SPM 计算 DCT 基函数个数：n = fix(2 * N * TR / hp_cutoff) + 1
n_dct = fix(2 * N * TR / hp_cutoff) + 1;
n_dct = max(n_dct, 1);  % 至少 1 个基函数

if exist('spm_dctmtx','file') == 2
    X0 = spm_dctmtx(N, n_dct);      % [N × n_dct]
else
    t = (0:N-1)';                   
    k = 0:(n_dct-1);                
    X0 = cos( (t + 0.5) * pi / N * k );   % [N × n_dct]
end

% 构造高通滤波矩阵：把 DCT 低频分量投影出去
H = eye(N) - X0 * pinv(X0);

results = struct();
results.cfg         = cfg;
results.roiNames    = roiNames;
results.roiHemi     = roiHemi;
results.roiBase     = roiBase;
results.baseROI_list= baseROI_list;
results.idxL        = idxL;
results.idxR        = idxR;
results.condLevels  = condLevels;
results.psych_names = psych_names;
results.psych_mat   = psych_mat;
results.runMeta     = runMeta;
results.confounds   = confounds;   % 如果太大，也可以之后不保存

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    if strcmp(condLabel,'connect')
        run_idx_cond = runIdx_connect;
    else
        run_idx_cond = runIdx_unconnect;
    end

    if ~isfield(roiSeed, condLabel)
        warning('roiSeed 中没有条件 %s，跳过', condLabel);
        continue;
    end

    nRun_cond = numel(run_idx_cond);

    % 预分配结果：beta_PPI, beta_phys, beta_psych, R2
    beta_PPI  = nan(nROI, nROI, nRun_cond);
    beta_phys = nan(nROI, nROI, nRun_cond);
    beta_psych= nan(nROI, nROI, nRun_cond);
    R2_mat    = nan(nROI, nROI, nRun_cond);

    for k = 1:nRun_cond
        % 对应的 confounds
        X_conf = confounds.(condLabel).X{k};  % [NTR × P] 或 []

        % 每个 ROI 的 seed 时序 (这一 run)
        Y_all = struct();
        for r = 1:nROI
            rn = roiNames{r};
            Y_all.(rn) = roiSeed.(condLabel).(rn)(:, k);  % [NTR × 1]
        end

        % 所有 ROI 两两 PPI
        for ii = 1:nROI
            rn_seed = roiNames{ii};
            phys    = Y_all.(rn_seed);   % 生理量（未滤波原始）

            % 预处理：z-score
            phys_z = zscore(phys);

            for jj = 1:nROI
                if ii == jj
                    continue;    % 自身→自身可以跳过
                end

                rn_targ = roiNames{jj};
                y       = Y_all.(rn_targ);  % 目标 ROI（未滤波原始）

                % 心理量：所有 run 共用同一条 Psych_diff（已 z-score）
                psych = psych_mat;

                % PPI 交互项（注意此时还未高通）
                ppi = phys_z .* psych;

                % ---- 构造设计矩阵：先分块 Z / N，再对 N 做正交化 ----
                % Z: [常数, psych, phys_z, ppi]
                Z = [ones(cfg.NTR_target,1), psych, phys_z, ppi];

                if cfg.use_confounds && ~isempty(X_conf)
                    % N0: 原始 confounds
                    N0 = X_conf;

                    % 将 confounds 在 Z 的子空间上回归并去除，使其对 Z 正交
                    % 注意：这里使用 \ 做最小二乘，可能本身有秩亏，但只影响
                    % confounds 的内部重参数化，不改变 Z 的解释空间
                    N_ortho = N0 - Z * (Z \ N0);

                    X = [Z, N_ortho];
                else
                    X = Z;
                end

                % ---------- 高通滤波：对 y 和所有回归量做相同的 DCT 高通 ----------
                y_filt = H * y;   % [NTR × 1]
                X_filt = H * X;   % [NTR × P]

                % ---------- 秩感知 OLS 拟合（QR + 列主元） ----------
                beta = safe_glm_rankaware(y_filt, X_filt);

                yhat = X_filt * beta;
                res  = y_filt - yhat;
                % R2 也在高通后的空间里计算
                R2   = 1 - var(res,0,1)/var(y_filt,0,1);

                % 记录：beta(4) = PPI 项的回归系数
                beta_PPI(ii,jj,k)   = beta(4);
                beta_phys(ii,jj,k)  = beta(3);
                beta_psych(ii,jj,k) = beta(2);
                R2_mat(ii,jj,k)     = R2;
            end
        end
    end

    results.(condLabel).beta_PPI   = beta_PPI;
    results.(condLabel).beta_phys  = beta_phys;
    results.(condLabel).beta_psych = beta_psych;
    results.(condLabel).R2         = R2_mat;
end

%% 5. 保存 run-wise 结果
out_name = fullfile(cfg.out_dir, ...
    sprintf('%s_ppi_roi_roi_runwise_fromOriginal.mat', cfg.subj));

save(out_name, 'results', '-v7.3');

fprintf('\n>>> 已保存 run-wise ROI-ROI PPI 结果到: %s\n', out_name);

%% 6. 被试级固定效应（跨 run + 只保留“优势半球”的 PPI 网络）
% 思路：
%   1) 在 ROI(L/R) 空间中，先对每个 ROI-pair (i,j) 在 run 维度上做固定效应平均
%   2) 用数量 pRF 的 voxel 数（meta.roiNkeep）来判定该被试整体的“优势半球”（L vs R）：
%        - 对所有 base ROI × 所有条件（connect / unconnect）的 voxel 数求和
%        - L_total > R_total -> 用 L 半球；R_total > L_total -> 用 R 半球
%        - 若两侧完全相等或 roiNkeep 不存在 -> 回退为双侧平均（保持原逻辑）
%   3) 在 base ROI 空间中，对每个 base-pair (B1,B2)：
%        - 根据“优势半球”只取 L(B1)->L(B2) 或 R(B1)->R(B2) 一条 ROI 边
%        - 该 ROI 边在 run 维度上的 FFX β 作为最终 β_B1→B2
%
% 注意：为兼容后续脚本，仍然使用字段名 beta_PPI_bilat，但现在表示“优势半球网络”。

fprintf('\n[Step 6] 计算被试级固定效应（跨 run + 优势半球）...\n');

results_fix = struct();
results_fix.cfg         = cfg;
results_fix.roiNames    = roiNames;
results_fix.roiHemi     = roiHemi;
results_fix.roiBase     = roiBase;
results_fix.baseROI_list= baseROI_list;
results_fix.idxL        = idxL;
results_fix.idxR        = idxR;
results_fix.condLevels  = condLevels;
results_fix.psych_names = psych_names;
results_fix.runMeta     = runMeta;
results_fix.source_runwise_file = out_name;  % 记录 run-wise 结果文件名

%% 6.0 基于 meta.roiNkeep 判定整体优势半球（L vs R）
domHemi_global = 'LR';   % 默认：LR -> 如果无法判定，就回退为“双侧平均”
domHemi_ratio  = NaN;
L_total_all    = NaN;
R_total_all    = NaN;

if exist('meta','var') && isfield(meta, 'roiNkeep')
    roiNkeep_struct = meta.roiNkeep;

    L_total_all = 0;
    R_total_all = 0;

    fprintf('\n[Step 6.0] 基于数量 pRF 体素数（roiNkeep）判定优势半球...\n');
    fprintf('  条件用于加总：%s\n', strjoin(condLevels, ', '));

    % 按 base ROI 打印各 ROI 左右 voxel 数
    for b = 1:numel(baseROI_list)
        baseName = baseROI_list{b};
        roiL = ['L' baseName];
        roiR = ['R' baseName];

        nL = 0;
        nR = 0;

        for ci = 1:numel(condLevels)
            condLabel = condLevels{ci};
            if ~isfield(roiNkeep_struct, condLabel)
                continue;
            end
            Skeep = roiNkeep_struct.(condLabel);

            if isfield(Skeep, roiL)
                nL = nL + double(Skeep.(roiL));
            end
            if isfield(Skeep, roiR)
                nR = nR + double(Skeep.(roiR));
            end
        end

        L_total_all = L_total_all + nL;
        R_total_all = R_total_all + nR;

        if nL + nR == 0
            fprintf('  ROI %-6s: L = %4d, R = %4d  (两侧均无有效 voxel)\n', baseName, nL, nR);
        else
            if nL > nR
                diff_ratio = (nL - nR) / (nL + nR) * 100;
                fprintf('  ROI %-6s: L = %4d, R = %4d  -> L 优势 (多 %.1f%% voxel)\n', ...
                    baseName, nL, nR, diff_ratio);
            elseif nR > nL
                diff_ratio = (nR - nL) / (nL + nR) * 100;
                fprintf('  ROI %-6s: L = %4d, R = %4d  -> R 优势 (多 %.1f%% voxel)\n', ...
                    baseName, nL, nR, diff_ratio);
            else
                fprintf('  ROI %-6s: L = %4d, R = %4d  -> L/R 相等\n', baseName, nL, nR);
            end
        end
    end

    if (L_total_all + R_total_all) > 0
        if L_total_all > R_total_all
            domHemi_global = 'L';
            domHemi_ratio  = (L_total_all - R_total_all) / (L_total_all + R_total_all);
        elseif R_total_all > L_total_all
            domHemi_global = 'R';
            domHemi_ratio  = (R_total_all - L_total_all) / (L_total_all + R_total_all);
        else
            domHemi_global = 'LR';
            domHemi_ratio  = 0;
        end

        fprintf('\n  [整体汇总] L_total = %d, R_total = %d\n', L_total_all, R_total_all);
        switch domHemi_global
            case 'L'
                fprintf('  ==> 优势半球：LEFT  (L 比 R 多 %.1f%% voxel)\n', domHemi_ratio*100);
            case 'R'
                fprintf('  ==> 优势半球：RIGHT (R 比 L 多 %.1f%% voxel)\n', domHemi_ratio*100);
            case 'LR'
                fprintf('  ==> L_total 与 R_total 完全相等，回退为 L/R 双侧平均模式。\n');
        end
    else
        fprintf('  [整体汇总] L_total 和 R_total 都为 0，无法判定优势半球，回退为双侧平均。\n');
        domHemi_global = 'LR';
    end
else
    fprintf('\n[Step 6.0] meta.roiNkeep 不存在，无法用 voxel 数判定优势半球，回退为双侧平均模式。\n');
    domHemi_global = 'LR';
end

results_fix.domHemi_global = domHemi_global;
results_fix.domHemi_ratio  = domHemi_ratio;
results_fix.L_total_vox    = L_total_all;
results_fix.R_total_vox    = R_total_all;
if exist('meta','var') && isfield(meta, 'roiNkeep')
    results_fix.roiNkeep = meta.roiNkeep;  % 方便之后检查
end

%% 6.1 按每个条件做 ROI 空间 FFX，然后投到 base ROI（只保留优势半球）
for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    if ~isfield(results, condLabel)
        warning('[FFX] results 中没有条件 %s，跳过', condLabel);
        continue;
    end

    beta_PPI  = results.(condLabel).beta_PPI;    % [nROI × nROI × nRun_cond]
    beta_phys = results.(condLabel).beta_phys;
    beta_psych= results.(condLabel).beta_psych;
    R2_mat    = results.(condLabel).R2;

    if isempty(beta_PPI)
        warning('[FFX] 条件 %s 的 beta_PPI 为空，跳过', condLabel);
        continue;
    end

    nRun_cond = size(beta_PPI,3);

    % ---- (1) ROI(L/R) 空间：跨 run 固定效应 ----
    beta_PPI_fix_roi   = mean(beta_PPI,  3, 'omitnan');   % [nROI × nROI]
    beta_phys_fix_roi  = mean(beta_phys, 3, 'omitnan');
    beta_psych_fix_roi = mean(beta_psych,3, 'omitnan');
    R2_fix_roi         = mean(R2_mat,    3, 'omitnan');

    % ---- 各 ROI-pair 实际参与的 run 数（非 NaN 数量）----
    Nruns_PPI_roi   = sum(~isnan(beta_PPI),  3);   % [nROI × nROI]
    Nruns_phys_roi  = sum(~isnan(beta_phys), 3);
    Nruns_psych_roi = sum(~isnan(beta_psych),3);
    Nruns_R2_roi    = sum(~isnan(R2_mat),    3);

    % ---- (2) base ROI 空间：只保留“优势半球”的边 ----
    nBase = numel(baseROI_list);
    beta_PPI_dom   = NaN(nBase, nBase);
    beta_phys_dom  = NaN(nBase, nBase);
    beta_psych_dom = NaN(nBase, nBase);
    R2_dom         = NaN(nBase, nBase);

    Nruns_PPI_dom   = NaN(nBase, nBase);
    Nruns_phys_dom  = NaN(nBase, nBase);
    Nruns_psych_dom = NaN(nBase, nBase);
    Nruns_R2_dom    = NaN(nBase, nBase);

    for b1 = 1:nBase
        for b2 = 1:nBase
            iL = idxL(b1);  jL = idxL(b2);
            iR = idxR(b1);  jR = idxR(b2);

            switch domHemi_global
                case 'L'
                    % 只用 L 半球，若 L 不存在则尝试 R（极端情况下）
                    if iL > 0 && jL > 0
                        beta_PPI_dom(b1,b2)   = beta_PPI_fix_roi(iL,jL);
                        beta_phys_dom(b1,b2)  = beta_phys_fix_roi(iL,jL);
                        beta_psych_dom(b1,b2) = beta_psych_fix_roi(iL,jL);
                        R2_dom(b1,b2)         = R2_fix_roi(iL,jL);

                        Nruns_PPI_dom(b1,b2)   = Nruns_PPI_roi(iL,jL);
                        Nruns_phys_dom(b1,b2)  = Nruns_phys_roi(iL,jL);
                        Nruns_psych_dom(b1,b2) = Nruns_psych_roi(iL,jL);
                        Nruns_R2_dom(b1,b2)    = Nruns_R2_roi(iL,jL);
                    elseif iR > 0 && jR > 0
                        % fallback：如果 L 上完全没有该 ROI，就退而求其次用 R
                        beta_PPI_dom(b1,b2)   = beta_PPI_fix_roi(iR,jR);
                        beta_phys_dom(b1,b2)  = beta_phys_fix_roi(iR,jR);
                        beta_psych_dom(b1,b2) = beta_psych_fix_roi(iR,jR);
                        R2_dom(b1,b2)         = R2_fix_roi(iR,jR);

                        Nruns_PPI_dom(b1,b2)   = Nruns_PPI_roi(iR,jR);
                        Nruns_phys_dom(b1,b2)  = Nruns_phys_roi(iR,jR);
                        Nruns_psych_dom(b1,b2) = Nruns_psych_roi(iR,jR);
                        Nruns_R2_dom(b1,b2)    = Nruns_R2_roi(iR,jR);
                    end

                case 'R'
                    % 只用 R 半球，若 R 不存在则尝试 L
                    if iR > 0 && jR > 0
                        beta_PPI_dom(b1,b2)   = beta_PPI_fix_roi(iR,jR);
                        beta_phys_dom(b1,b2)  = beta_phys_fix_roi(iR,jR);
                        beta_psych_dom(b1,b2) = beta_psych_fix_roi(iR,jR);
                        R2_dom(b1,b2)         = R2_fix_roi(iR,jR);

                        Nruns_PPI_dom(b1,b2)   = Nruns_PPI_roi(iR,jR);
                        Nruns_phys_dom(b1,b2)  = Nruns_phys_roi(iR,jR);
                        Nruns_psych_dom(b1,b2) = Nruns_psych_roi(iR,jR);
                        Nruns_R2_dom(b1,b2)    = Nruns_R2_roi(iR,jR);
                    elseif iL > 0 && jL > 0
                        % fallback：如果 R 上完全没有该 ROI，就退而求其次用 L
                        beta_PPI_dom(b1,b2)   = beta_PPI_fix_roi(iL,jL);
                        beta_phys_dom(b1,b2)  = beta_phys_fix_roi(iL,jL);
                        beta_psych_dom(b1,b2) = beta_psych_fix_roi(iL,jL);
                        R2_dom(b1,b2)         = R2_fix_roi(iL,jL);

                        Nruns_PPI_dom(b1,b2)   = Nruns_PPI_roi(iL,jL);
                        Nruns_phys_dom(b1,b2)  = Nruns_phys_roi(iL,jL);
                        Nruns_psych_dom(b1,b2) = Nruns_psych_roi(iL,jL);
                        Nruns_R2_dom(b1,b2)    = Nruns_R2_roi(iL,jL);
                    end

                case 'LR'
                    % 回退为原始“双侧平均”逻辑
                    vals_PPI  = [];
                    vals_phys = [];
                    vals_psych= [];
                    vals_R2   = [];

                    n_PPI  = [];
                    n_phys = [];
                    n_psych= [];
                    n_R2   = [];

                    if iL > 0 && jL > 0
                        vals_PPI(end+1)   = beta_PPI_fix_roi(iL,jL);
                        vals_phys(end+1)  = beta_phys_fix_roi(iL,jL);
                        vals_psych(end+1) = beta_psych_fix_roi(iL,jL);
                        vals_R2(end+1)    = R2_fix_roi(iL,jL);

                        n_PPI(end+1)   = Nruns_PPI_roi(iL,jL);
                        n_phys(end+1)  = Nruns_phys_roi(iL,jL);
                        n_psych(end+1) = Nruns_psych_roi(iL,jL);
                        n_R2(end+1)    = Nruns_R2_roi(iL,jL);
                    end
                    if iR > 0 && jR > 0
                        vals_PPI(end+1)   = beta_PPI_fix_roi(iR,jR);
                        vals_phys(end+1)  = beta_phys_fix_roi(iR,jR);
                        vals_psych(end+1) = beta_psych_fix_roi(iR,jR);
                        vals_R2(end+1)    = R2_fix_roi(iR,jR);

                        n_PPI(end+1)   = Nruns_PPI_roi(iR,jR);
                        n_phys(end+1)  = Nruns_phys_roi(iR,jR);
                        n_psych(end+1) = Nruns_psych_roi(iR,jR);
                        n_R2(end+1)    = Nruns_R2_roi(iR,jR);
                    end

                    if ~isempty(vals_PPI)
                        beta_PPI_dom(b1,b2)   = mean(vals_PPI,   'omitnan');
                        beta_phys_dom(b1,b2)  = mean(vals_phys,  'omitnan');
                        beta_psych_dom(b1,b2) = mean(vals_psych, 'omitnan');
                        R2_dom(b1,b2)         = mean(vals_R2,    'omitnan');

                        Nruns_PPI_dom(b1,b2)   = sum(n_PPI);
                        Nruns_phys_dom(b1,b2)  = sum(n_phys);
                        Nruns_psych_dom(b1,b2) = sum(n_psych);
                        Nruns_R2_dom(b1,b2)    = sum(n_R2);
                    end
            end
        end
    end

    % 将优势半球结果存入 results_fix，字段名保持为 *_bilat 以兼容下游脚本
    results_fix.(condLabel).beta_PPI_roi   = beta_PPI_fix_roi;
    results_fix.(condLabel).beta_phys_roi  = beta_phys_fix_roi;
    results_fix.(condLabel).beta_psych_roi = beta_psych_fix_roi;
    results_fix.(condLabel).R2_roi         = R2_fix_roi;

    results_fix.(condLabel).Nruns_PPI_roi   = Nruns_PPI_roi;
    results_fix.(condLabel).Nruns_phys_roi  = Nruns_phys_roi;
    results_fix.(condLabel).Nruns_psych_roi = Nruns_psych_roi;
    results_fix.(condLabel).Nruns_R2_roi    = Nruns_R2_roi;

    results_fix.(condLabel).beta_PPI_bilat   = beta_PPI_dom;
    results_fix.(condLabel).beta_phys_bilat  = beta_phys_dom;
    results_fix.(condLabel).beta_psych_bilat = beta_psych_dom;
    results_fix.(condLabel).R2_bilat         = R2_dom;

    results_fix.(condLabel).Nruns_PPI_bilat   = Nruns_PPI_dom;
    results_fix.(condLabel).Nruns_phys_bilat  = Nruns_phys_dom;
    results_fix.(condLabel).Nruns_psych_bilat = Nruns_psych_dom;
    results_fix.(condLabel).Nruns_R2_bilat    = Nruns_R2_dom;

    results_fix.(condLabel).nRun_cond = nRun_cond;
end

%% 6.2 保存优势半球 FFX 结果
out_name_ffx = fullfile(cfg.out_dir, ...
    sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral_domHemi.mat', cfg.subj));

save(out_name_ffx, 'results_fix', '-v7.3');

fprintf('\n>>> 已保存被试级固定效应（优势半球网络）结果到: %s\n', out_name_ffx);
fprintf('================ ROI-ROI PPI (FFX, within-subject, dominant hemisphere) 完成 ================\n');

%% 7. 可视化：基于优势半球 β 的单被试 ROI-ROI PPI 脑网络（含 FDR 校正）
% 使用：
%   - 边强度：results_fix.(cond).beta_PPI_bilat（现在表示优势半球 base ROI 网络）
%   - 统计显著性：仅基于该优势半球的 run-wise beta_PPI（不再合并 L/R）
%     例如若 domHemi_global = 'L'，则只用 L(B1)->L(B2) 的所有 run-wise β 做 t 检验 vs 0；
%   - 再对 base ROI-pair 的 p 做 FDR (BH)，并按正/负向 + (未)校正显著性着色

fprintf('\n[Step 7] 可视化优势半球固定效应后的 ROI-ROI PPI 脑网络（含 FDR 校正）...\n');

alpha_unc = 0.05;   % 未校正显著性阈值
q_fdr     = 0.05;   % FDR 控制水准（Benjamini-Hochberg）

% 颜色定义
col_red    = [1   0     0  ];  % FDR 正向显著
col_blue   = [0   0.6   1  ];  % FDR 负向显著
col_orange = [1   0.5   0  ];  % 未校正正向显著但 FDR 不显著
col_green  = [0   1     0  ];  % 未校正负向显著但 FDR 不显著
col_white  = [1   1     1  ];  % 完全不显著

condLevels_fix = results_fix.condLevels;
baseROI_list   = results_fix.baseROI_list;
idxL           = results_fix.idxL;
idxR           = results_fix.idxR;
nBase          = numel(baseROI_list);

for ci = 1:numel(condLevels_fix)
    condLabel = condLevels_fix{ci};

    if ~isfield(results_fix, condLabel)
        warning('[Vis] results_fix 中没有条件 %s，跳过可视化。', condLabel);
        continue;
    end

    % 优势半球固定效应后的有向 β 矩阵（base ROI 空间）
    if ~isfield(results_fix.(condLabel), 'beta_PPI_bilat') || ...
       isempty(results_fix.(condLabel).beta_PPI_bilat)
        warning('[Vis] 条件 %s 的 beta_PPI_bilat 为空，跳过可视化。', condLabel);
        continue;
    end
    beta_fix_dom = results_fix.(condLabel).beta_PPI_bilat;  % [nBase × nBase]

    % run-wise β 来源：原始 results.(condLabel).beta_PPI（ROI 空间）
    if ~isfield(results, condLabel) || ~isfield(results.(condLabel), 'beta_PPI')
        warning('[Vis] results 中缺少条件 %s 的 run-wise beta_PPI，跳过可视化。', condLabel);
        continue;
    end
    beta_run_roi = results.(condLabel).beta_PPI;   % [nROI × nROI × nRun_cond]
    nRun_cond    = size(beta_run_roi,3);

    % ---------- 对每个 base ROI-pair，基于优势半球的 run-wise β 做 t 检验 ----------
    p_bilat_dir  = NaN(nBase, nBase);   % 有向 p 值
    t_bilat_dir  = NaN(nBase, nBase);   % t 统计量
    n_obs_dir    = NaN(nBase, nBase);   % 实际 run 数

    for b1 = 1:nBase
        for b2 = 1:nBase
            if b1 == b2
                continue;
            end

            iL = idxL(b1);  jL = idxL(b2);
            iR = idxR(b1);  jR = idxR(b2);

            % 选取优势半球对应的 run-wise β
            switch domHemi_global
                case 'L'
                    if iL <= 0 || jL <= 0
                        continue;
                    end
                    beta_run_vec = squeeze(beta_run_roi(iL, jL, :));  % [nRun_cond × 1]
                case 'R'
                    if iR <= 0 || jR <= 0
                        continue;
                    end
                    beta_run_vec = squeeze(beta_run_roi(iR, jR, :));
                case 'LR'
                    % 回退模式：合并 L 和 R 的 run-wise β
                    beta_run_vec = [];
                    if iL > 0 && jL > 0
                        beta_run_vec = [beta_run_vec; squeeze(beta_run_roi(iL, jL, :))];
                    end
                    if iR > 0 && jR > 0
                        beta_run_vec = [beta_run_vec; squeeze(beta_run_roi(iR, jR, :))];
                    end
            end

            beta_run_vec = beta_run_vec(:);
            beta_run_vec = beta_run_vec(~isnan(beta_run_vec));

            if numel(beta_run_vec) < 2
                continue;
            end

            [~, p_val, ~, stats] = ttest(beta_run_vec, 0, 'Alpha', alpha_unc);
            p_bilat_dir(b1,b2) = p_val;
            t_bilat_dir(b1,b2) = stats.tstat;
            n_obs_dir(b1,b2)   = numel(beta_run_vec);
        end
    end

    % ---------- 构建对称 β / p 矩阵，用于无向网络可视化 ----------
    beta_sym = NaN(nBase, nBase);
    p_sym    = NaN(nBase, nBase);

    for i = 1:nBase
        for j = i+1:nBase
            % β：对 B1->B2 和 B2->B1 的 β 取平均
            pair_beta  = [beta_fix_dom(i,j), beta_fix_dom(j,i)];
            valid_beta = pair_beta(~isnan(pair_beta));
            if ~isempty(valid_beta)
                b_ij = mean(valid_beta);
                beta_sym(i,j) = b_ij;
                beta_sym(j,i) = b_ij;
            end

            % p：对双向 p 取最小值（更保守地选择更显著的一边）
            pair_p  = [p_bilat_dir(i,j), p_bilat_dir(j,i)];
            valid_p = pair_p(~isnan(pair_p));
            if ~isempty(valid_p)
                p_ij = min(valid_p);
                p_sym(i,j) = p_ij;
                p_sym(j,i) = p_ij;
            end
        end
        beta_sym(i,i) = 0;
        p_sym(i,i)    = NaN;
    end

    % 若全部 NaN，则不画图
    if all(isnan(beta_sym(:)))
        warning('[Vis] 条件 %s 所有 β 为 NaN，跳过绘图。', condLabel);
        continue;
    end

    % ---------- 计算 FDR 校正（对每个条件单独做） ----------
    mask_ut = triu(~isnan(p_sym), 1);  % 只看上三角
    p_vec   = p_sym(mask_ut);
    m_edge  = numel(p_vec);

    sig_fdr_mat = false(size(p_sym));
    p_crit = NaN;

    if m_edge > 0 && any(~isnan(p_vec))
        [h_fdr_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_mat(mask_ut) = h_fdr_vec;
        sig_fdr_mat = sig_fdr_mat | sig_fdr_mat.';  % 对称化

        fprintf('[Vis] 条件 %s：FDR q=%.3f，%d/%d 条边通过 FDR (p<=%.3g)\n', ...
            condLabel, q_fdr, sum(h_fdr_vec), m_edge, p_crit);
    else
        fprintf('[Vis] 条件 %s：没有可用于 FDR 的有效 p 值。\n', condLabel);
    end

    % 未校正显著性矩阵（两侧 p<alpha_unc）
    sig_unc_mat = (p_sym < alpha_unc) & ~isnan(p_sym);
    sig_unc_mat = sig_unc_mat | sig_unc_mat.';

    % ---------- 节点坐标：均匀放在一个圆上 ----------
    theta = linspace(0, 2*pi, nBase+1);
    theta = theta(1:nBase);
    R     = 1;

    x = R * cos(theta);
    y = R * sin(theta);

    % ---------- 计算最大 |beta| 用于设置线宽 ----------
    valid_abs_beta = abs(beta_sym(~isnan(beta_sym) & (triu(true(size(beta_sym)),1))));
    max_abs_beta   = max(valid_abs_beta);
    if isempty(max_abs_beta) || max_abs_beta == 0
        max_abs_beta = 1;
    end

    % ---------- 画图 ----------
    fig_name = sprintf('%s - %s PPI network (FFX dominant hemi, Psych_diff, FDR)', ...
        cfg.subj, condLabel);
    figure('Name', fig_name);
    set(gcf, 'Color', [0 0 0]);
    ax = gca;
    set(ax, 'Color', [0 0 0]);
    hold on; axis equal off;

    % 画边（只画上三角）
    for i = 1:nBase-1
        for j = i+1:nBase
            if isnan(beta_sym(i,j))
                continue;
            end

            b_ij = beta_sym(i,j);
            p_ij = p_sym(i,j);

            is_fdr = sig_fdr_mat(i,j);
            is_unc = sig_unc_mat(i,j);

            % 默认颜色：白色（整体不显著）
            edge_color = col_white;

            if is_fdr
                % FDR 显著
                if b_ij > 0
                    edge_color = col_red;
                elseif b_ij < 0
                    edge_color = col_blue;
                end
            elseif is_unc
                % 未经校正显著但 FDR 不显著
                if b_ij > 0
                    edge_color = col_orange;
                elseif b_ij < 0
                    edge_color = col_green;
                end
            else
                % 未经校正也不显著：保持白色
                edge_color = col_white;
            end

            % 线宽随 |beta| 变化
            w = 1.5 + 3 * (abs(b_ij) / max_abs_beta);

            % 画无向边
            line([x(i) x(j)], [y(i) y(j)], ...
                'Color', edge_color, 'LineWidth', w);

            % 在线的中点标 β 值
            xm = (x(i) + x(j)) / 2;
            ym = (y(i) + y(j)) / 2;
            text(xm, ym, sprintf('%.3f', b_ij), ...
                'Color', edge_color, ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'none');
        end
    end

    % 画节点
    node_size = 300;
    for i = 1:nBase
        scatter(x(i), y(i), node_size, ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5);
        text(x(i)*1.2, y(i)*1.2, baseROI_list{i}, ...
            'Color', [1 1 1], ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end

    title(sprintf('%s - %s: PPI network (FFX dominant hemi, numerosity vs baseline, FDR)', ...
        cfg.subj, condLabel), ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('\n[Step 7] 优势半球固定效应后的脑网络可视化完成。\n');



%% ========== 子函数 ==========

function [psych_mat, psych_names] = build_psych(cfg)
% 从 cfg.events_tsv 构造 Psych_diff (numerosity vs baseline)

assert(exist(cfg.events_tsv,'file')==2, 'events.tsv 不存在: %s', cfg.events_tsv);

Te = readtable(cfg.events_tsv, 'FileType','text', 'Delimiter','\t');
fprintf('    events.tsv: %d 行, %d 列\n', height(Te), width(Te));

req_cols = {'onset','duration','trial_type'};
missing  = setdiff(req_cols, Te.Properties.VariableNames);
assert(isempty(missing), 'events.tsv 缺少列: %s', strjoin(missing, ', '));

t = (0:cfg.NTR_target-1)' * cfg.TR;

% HRF
if exist('spm_hrf','file')==2
    fprintf('    使用 spm_hrf.\n');
    hrf = spm_hrf(cfg.TR);
else
    fprintf('    找不到 spm_hrf，使用简易 gamma HRF.\n');
    tt  = (0:1:32)';
    hrf = tt.^8 .* exp(-tt/0.9);
    hrf = hrf / sum(hrf);
    hrf = interp1(tt, hrf, 0:cfg.TR:32, 'pchip', 0)';
end

cname_num  = cfg.cond_names{1};
cname_base = cfg.cond_names{2};

% numerosity
mask_num = strcmp(Te.trial_type, cname_num);
assert(any(mask_num), 'trial_type 中找不到 "%s"', cname_num);
onset_num= Te.onset(mask_num);
dur_num  = Te.duration(mask_num);

box_num = zeros(cfg.NTR_target,1);
for k = 1:numel(onset_num)
    t_on  = onset_num(k);
    t_off = onset_num(k) + dur_num(k);
    box_num = box_num + (t >= t_on & t < t_off);
end

% baseline
mask_base = strcmp(Te.trial_type, cname_base);
assert(any(mask_base), 'trial_type 中找不到 "%s"', cname_base);
onset_base= Te.onset(mask_base);
dur_base  = Te.duration(mask_base);

box_base = zeros(cfg.NTR_target,1);
for k = 1:numel(onset_base)
    t_on  = onset_base(k);
    t_off = onset_base(k) + dur_base(k);
    box_base = box_base + (t >= t_on & t < t_off);
end

diff_box  = box_num - box_base;
conv_diff = conv(diff_box, hrf);
conv_diff = conv_diff(1:cfg.NTR_target);

% z-score
conv_diff = conv_diff - mean(conv_diff, 'omitnan');
sd_diff   = std(conv_diff, 0, 'omitnan');
if sd_diff == 0
    warning('    Psych_diff 方差为 0，置为 0。');
    conv_diff(:) = 0;
else
    conv_diff = conv_diff / sd_diff;
end

psych_mat   = conv_diff;
psych_names = {'Psych_diff'};

fprintf('    已生成 Psych_diff（行数=%d）。\n', numel(conv_diff));
end

function [Xc, cnames, meta_out] = build_confounds_one_run(conf_tsv, bold_nii, ...
    wm_nii, csf_nii, cfg)
% 为单个 run 构建混杂项矩阵 Xc （motion+dmotion+FD_pulse+aCompCor）
% 现在：aCompCor 直接使用 fMRIPrep tsv 中已有的 w_comp_cor_* / c_comp_cor_*
% bold_nii, wm_nii, csf_nii 不再参与 aCompCor 计算，仅保留接口一致性

% --- 读 confounds 表 ---
Tc = readtable(conf_tsv, 'FileType','text', 'Delimiter','\t');
N_full = height(Tc);

% 确定时间裁剪索引
if N_full == cfg.NTR_target
    idx = (1:cfg.NTR_target)';
elseif N_full >= cfg.n_dummy + cfg.NTR_target
    idx = (cfg.n_dummy+1) : (cfg.n_dummy + cfg.NTR_target);
else
    error('confounds 行数(%d) 太少，无法通过 dummy(%d) 获得 %d 帧', ...
        N_full, cfg.n_dummy, cfg.NTR_target);
end

Tc_use = Tc(idx, :);   % 裁剪后的 confounds 时序

% --------- 运动参数 6 + 导数 6 ---------
motNames = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
hasAll = all(ismember(motNames, Tc_use.Properties.VariableNames));
assert(hasAll, 'confounds 中缺少运动列 (trans/rot)');

mot = table2array(Tc_use(:, motNames));   % [NTR × 6]
dmot = [zeros(1,6); diff(mot,1,1)];       % 一阶导数

% --------- FD 脉冲列 ---------
if ismember('framewise_displacement', Tc_use.Properties.VariableNames)
    FD = Tc_use.framewise_displacement;
    FD(isnan(FD)) = 0;
    FD_pulse = double(FD > cfg.fd_thr);
else
    warning('confounds 中没有 framewise_displacement 列，FD 脉冲列用 0 填充');
    FD_pulse = zeros(cfg.NTR_target,1);
end

% --------- 直接使用 fMRIPrep 的 WM / CSF aCompCor ---------
allVars = Tc_use.Properties.VariableNames;

% WM 成分：w_comp_cor_00, w_comp_cor_01, ...
wmColsAll = allVars(startsWith(allVars, 'w_comp_cor_'));
wmColsAll = sort(wmColsAll);   % 有前导 0，字典序就是正确顺序

% CSF 成分：c_comp_cor_00, c_comp_cor_01, ...
csfColsAll = allVars(startsWith(allVars, 'c_comp_cor_'));
csfColsAll = sort(csfColsAll);

% 取前 cfg.nComp_WM / cfg.nComp_CSF 个
if isempty(wmColsAll)
    warning('confounds 中没有 w_comp_cor_* 列，WM aCompCor 将用 0 填充');
    score_wm = zeros(cfg.NTR_target, cfg.nComp_WM);
    wmCols = arrayfun(@(k) sprintf('WM_aCompCor_%02d',k), 1:cfg.nComp_WM, 'uni', false);
else
    nW = min(cfg.nComp_WM, numel(wmColsAll));
    wmCols = wmColsAll(1:nW);
    score_wm = table2array(Tc_use(:, wmCols));
end

if isempty(csfColsAll)
    warning('confounds 中没有 c_comp_cor_* 列，CSF aCompCor 将用 0 填充');
    score_csf = zeros(cfg.NTR_target, cfg.nComp_CSF);
    csfCols = arrayfun(@(k) sprintf('CSF_aCompCor_%02d',k), 1:cfg.nComp_CSF, 'uni', false);
else
    nC = min(cfg.nComp_CSF, numel(csfColsAll));
    csfCols = csfColsAll(1:nC);
    score_csf = table2array(Tc_use(:, csfCols));
end

% ==== 合并所有列：motion + dmotion + FD_pulse + WM/CSF aCompCor ====
Xc_raw = [mot, dmot, FD_pulse, score_wm, score_csf];

cnames_raw = [ ...
    motNames, ...
    strcat('d_', motNames), ...
    {'FD_pulse'}, ...
    wmCols, ...
    csfCols ...
];

% 去均值
Xc_raw = Xc_raw - mean(Xc_raw,1,'omitnan');

% 检查常数（std = 0）的列，并丢掉
colStd = std(Xc_raw, 0, 1, 'omitnan');
keep   = colStd > 0;

if any(~keep)
    fprintf('    [Confounds] 本 run 中下列混杂列为常数/全0，将被丢弃:\n');
    badIdx = find(~keep);
    for ii = 1:numel(badIdx)
        fprintf('        %s (std = %.4g)\n', cnames_raw{badIdx(ii)}, colStd(badIdx(ii)));
    end
end

Xc     = Xc_raw(:, keep);
cnames = cnames_raw(keep);

% meta 信息（记录原始情况）
meta_out = struct();
meta_out.N_full     = N_full;
meta_out.idx        = idx;
meta_out.motNames   = motNames;
meta_out.fd_thr     = cfg.fd_thr;
meta_out.wmCols     = wmCols;
meta_out.csfCols    = csfCols;
meta_out.source     = 'fMRIPrep_confounds_tsv';
meta_out.constColsDropped = cnames_raw(~keep);  % 被丢弃的常数列名（可能为空）
end

%% 安全的秩感知 GLM 求解函数
function beta_full = safe_glm_rankaware(y, X)
    % 安全求解 beta，自动处理秩亏（QR + 列主元）
    % y: T × 1 时间序列
    % X: T × P 设计矩阵（含 PPI, psych, phys, nuisance, 常数等）

    [Q,R,perm] = qr(X,0);  % 带列主元 QR
    d   = abs(diag(R));
    if isempty(d)
        beta_full = zeros(size(X,2),1);
        return;
    end

    tol = max(size(X)) * eps(max(d));
    r   = sum(d > tol);   % 有效秩

    if r == 0
        % 所有列都几乎为 0，返回 0 解
        warning('safe_glm_rankaware: design matrix rank = 0，返回全 0 beta。');
        beta_full = zeros(size(X,2),1);
        return;
    end

    if r < size(X,2)
        warning('safe_glm_rankaware: design matrix rank deficient: P = %d, rank = %d, drop %d col(s).',...
                size(X,2), r, size(X,2)-r);
    end

    % 在秩 r 的子空间上求解
    Q1 = Q(:,1:r);
    R1 = R(1:r,1:r);
    beta_red = R1 \ (Q1' * y);

    % 把 beta 映射回原始列顺序
    beta_full = zeros(size(X,2),1);
    beta_full(perm(1:r)) = beta_red;
end

%% 本地 FDR 函数（Benjamini-Hochberg）
function [h, p_crit] = fdr_bh_local(p_vals, q)
    % p_vals: 向量形式的 p 值（可含 NaN）
    % q     : FDR 控制水平，例如 0.05
    %
    % 返回：
    %   h      : 与 p_vals 同维度的逻辑向量，true 表示通过 FDR
    %   p_crit : 阈值 p*，即最大的 p，使得 p(k) <= (k/m)*q

    p = p_vals(:);
    m = numel(p);
    h = false(size(p));

    valid = ~isnan(p);
    pv    = p(valid);
    if isempty(pv)
        p_crit = NaN;
        h(~valid) = false;
        h = reshape(h, size(p_vals));
        return;
    end

    [pv_sorted, ~] = sort(pv);
    m_valid = numel(pv_sorted);
    crit    = (1:m_valid)' * (q / m_valid);

    cmp = pv_sorted <= crit;
    if any(cmp)
        k_max   = find(cmp, 1, 'last');
        p_crit  = pv_sorted(k_max);
        sig_idx = pv <= p_crit & valid;
        h(sig_idx) = true;
    else
        p_crit = 0;
    end

    h = reshape(h, size(p_vals));
end