%% 2_ppi_roi_roi_runwise_fromOriginal_cmd.m
% 从 Gray/Original 的 run-wise ROI 时序 + fMRIPrep confounds 做 ROI-ROI PPI
%   - 每个 run 单独建模（run-wise）
%   - confounds 使用 fMRIPrep 输出（motion, FD, aCompCor WM/CSF）
%   - 本版本：增加左右半球固定效应合并 + 秩感知 QR 求解 + confounds 对 Z 的正交化：
%       * 在 ROI(含 L/R) 空间中做 run-wise PPI
%       * 设计矩阵 X = [Z, N_ortho]，其中 Z = [1, psych, phys, ppi]，
%         N_ortho 为 confounds 对 Z 正交化后的部分（不改变 PPI 本身）
%       * 用 safe_glm_rankaware 在秩 r 子空间求解 GLM，自动处理秩亏
%       * 在被试内部先对每个半球单独跨 run 平均，得到 β_L、β_R
%       * 再对 β_L、β_R 简单平均，得到 base ROI–base ROI 级别的 bilateral β

clear; clc;
fprintf('\n================ ROI-ROI PPI (run-wise, from Original, bilateral FFX) ================\n');

%% 0. 基本配置（只改数据/路径/映射/命名）
cfg = struct();

% ------------------ 双 ID：读取 sub-03，保存/命名 sub-08 ------------------
cfg.data_subj  = 'sub-03';   % 实际数据所在被试
cfg.save_subj  = 'sub-08';   % 输出记录/命名用的被试号（不影响读取）
cfg.subj       = cfg.data_subj;  % 保持兼容（脚本内部读取仍用 cfg.subj）

% mrVista 工程根（Gray/Original 在这里）
cfg.prf_ses    = 'ses-01';  % 旧实验 Gray/Original 合并的会话名（mrVista 里）
cfg.root_prf   = '/home/xue/daifenpei/sub-03/results1';

% 输出 derivatives 的根（你想放在 results3 下）
cfg.root_bids  = '/home/xue/daifenpei/sub-03/results3';

% 统一的“被试根目录”（用于 resultsX ↔ ses-0X 自动映射）
% 例如：/home/xue/daifenpei/sub-03
cfg.root_base  = fileparts(cfg.root_prf);  % /home/xue/daifenpei/sub-03

% T1 分割概率图优先使用的 session（如果缺失，会自动在 ses-01..ses-04 里找）
cfg.t1_ses     = 'ses-03';

% ROI multi-run 原始时序 .mat（旧实验专用 ROI 提取脚本生成）
cfg.roi_ts_dir = 'results_ppi_roi_tseries_fromOriginal';
cfg.roi_ts_file= sprintf('%s_%s_Original_roiTS_multiRun_dualPRF.mat', ...
                         cfg.data_subj, cfg.prf_ses);

% events（所有 run 共用）
cfg.events_tsv = '/media/xue/new_B1/glm/events_output/events_bids.tsv';

% 时间信息
cfg.TR         = 1.95;
cfg.n_dummy    = 6;         % 原始 BOLD / confounds 前面的 dummy 数量
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

% 输出目录（放在 results3/derivatives 下，用 save_subj 命名）
cfg.out_dir = fullfile(cfg.root_bids, 'derivatives', 'ppi_roi_runwise', cfg.save_subj);
if ~exist(cfg.out_dir, 'dir')
    fprintf('  [Info] 创建输出目录: %s\n', cfg.out_dir);
    mkdir(cfg.out_dir);
end

fprintf('  [Config] data_subj=%s, save_subj=%s\n', cfg.data_subj, cfg.save_subj);
fprintf('  [Config] root_prf=%s\n', cfg.root_prf);
fprintf('  [Config] root_bids(for output)=%s\n', cfg.root_bids);
fprintf('  [Config] TR=%.3f, NTR_target=%d, n_dummy=%d\n', ...
    cfg.TR, cfg.NTR_target, cfg.n_dummy);
fprintf('  [Config] cond_names: %s\n', strjoin(cfg.cond_names, ', '));
fprintf('  [Config] use_confounds = %d\n', cfg.use_confounds);

%% 1. 载入 ROI 时序 (run-wise, Original) 并解析 L/R + base ROI 信息
fprintf('\n[Step 1] 载入 run-wise ROI 时序 (Original)...\n');

sess_dir  = fullfile(cfg.root_prf, cfg.data_subj, cfg.prf_ses);
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

runIdx_connect   = meta.runIdx_connect;
runIdx_unconnect = meta.runIdx_unconnect;

fprintf('    [Info] ROI 数 = %d\n', nROI);
fprintf('    [Info] 总扫描数 (scanInfo) = %d\n', nScan);
fprintf('    [Info] runIdx_connect   = %s\n', mat2str(runIdx_connect));
fprintf('    [Info] runIdx_unconnect = %s\n', mat2str(runIdx_unconnect));

% 解析 ROI 的 hemis/base 名
roiHemi = cell(nROI,1);
roiBase = cell(nROI,1);
for i = 1:nROI
    rn = roiNames{i};
    if isempty(rn)
        roiHemi{i} = '';
        roiBase{i} = '';
        continue;
    end

    % 只看第一个字符，大小写不敏感：'L'/'l' 视为左半球，'R'/'r' 视为右半球
    h = upper(rn(1));
    if h == 'L'
        roiHemi{i} = 'L';
        roiBase{i} = rn(2:end);  % 去掉首字母
    elseif h == 'R'
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
% sub-03 映射规则（按你给的合并逻辑）：
%   Scan1..8   -> ses-01, run1..8
%   Scan9..17  -> ses-02, run1..9
%   Scan18..27 -> ses-03, run1..10
%   Scan28..36 -> ses-04, run1..9   （本次 PPI 不用，但映射函数仍可支持）
%
% 注意：你的 meta.scanInfo 里 nScan=18，但 scanNums 是“原始 Gray Scan 号”，不一定连续
runMeta = struct('scanIdx', cell(1,nScan), ...
                 'scanNum', [], 'sessionIdx', [], 'cond', [], ...
                 'runInSes', [], 'sesLabel', [], 'runLabel', []);

fprintf('\n[Step 1.1] 构造 runMeta（Gray Scan → BIDS ses/run）...\n');

for i = 1:nScan
    sNum   = scanInfo.scanNums(i);   % Gray Scan 号（例如 2,3,4,5,...,27）
    cond_i = scanInfo.cond{i};       % 'connect' / 'unconnect'

    rm = struct();
    rm.scanIdx    = i;
    rm.scanNum    = sNum;
    rm.cond       = cond_i;
    rm.sessionIdx = NaN;
    rm.runInSes   = NaN;
    rm.sesLabel   = '';
    rm.runLabel   = '';

    [sesIdx, runInSes] = scanNum_to_sesrun_sub03style(sNum);

    if isnan(sesIdx) || isnan(runInSes)
        fprintf('    [Warn] Scan%02d 无法映射到 ses/run（超出 1..36 或不在规则范围）。\n', sNum);
        runMeta(i) = rm;
        continue;
    end

    rm.sessionIdx = sesIdx;
    rm.runInSes   = runInSes;
    rm.sesLabel   = sprintf('ses-%02d', sesIdx);
    rm.runLabel   = sprintf('run-%03d', runInSes); % 仅用于日志

    runMeta(i) = rm;

    fprintf('    Scan%2d -> %s, %s, cond=%s\n', ...
        sNum, rm.sesLabel, rm.runLabel, cond_i);
end

% 确认 meta.runIdx_connect/unconnect 中引用的 scan 都成功映射
used_idx = unique([runIdx_connect(:); runIdx_unconnect(:)]);
for uu = used_idx(:)'
    assert(uu >= 1 && uu <= nScan, 'runIdx (%d) 超出 scan 范围 1..%d', uu, nScan);
    rm_chk = runMeta(uu);
    assert(~isnan(rm_chk.runInSes), ...
        'scanIdx=%d (scanNum=%d) 被用于 PPI，但未成功映射到 sessionIdx/runInSes。', ...
        uu, rm_chk.scanNum);
end
fprintf('    [Check] 所有用于 PPI 的 runIdx_connect/unconnect 均成功映射到 ses-01/02/03(/04) & run。\n');

%% 2. 构造心理回归量 Psych_diff（numerosity vs baseline），每个 run 共用
fprintf('\n[Step 2] 从 events_bids.tsv 构造 Psych_diff...\n');
[psych_mat, psych_names] = build_psych(cfg);
NTR_psych = size(psych_mat,1);
assert(NTR_psych == cfg.NTR_target, ...
    'Psych_diff 长度(%d) != cfg.NTR_target(%d)', NTR_psych, cfg.NTR_target);

%% 3. 为每个 run 构建 fMRIPrep confounds（motion + FD + aCompCor）
fprintf('\n[Step 3] 为每个 run 构建 fMRIPrep 混杂项...\n');

% WM / CSF 概率图路径（T1 space）
[wm_nii, csf_nii] = resolve_wm_csf_probseg(cfg);

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
        i = run_idx_cond(k);      % 全部 scanInfo 条目中的索引
        rm = runMeta(i);

        % func 根目录：resultsX 对应 ses-0X
        sesFuncRoot = get_func_root_for_session(rm.sesLabel, cfg);
        sesDir      = fullfile(sesFuncRoot, cfg.data_subj, rm.sesLabel, 'func');

        % ---- 关键修复：兼容 task-run10 vs task-run010（以及 task-run-XX） ----
        [conf_tsv, bold_nii, baseStub_used] = resolve_fmriprep_run_files( ...
            cfg.data_subj, rm.sesLabel, rm.runInSes, sesDir);

        fprintf('        Run %d/%d: %s, %s\n', k, nRun_cond, rm.sesLabel, rm.runLabel);
        fprintf('            stub     : %s\n', baseStub_used);
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
n_dct = max(n_dct, 1);

if exist('spm_dctmtx','file') == 2
    X0 = spm_dctmtx(N, n_dct);
else
    t = (0:N-1)';
    k = 0:(n_dct-1);
    X0 = cos( (t + 0.5) * pi / N * k );
end

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
results.confounds   = confounds;

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

    beta_PPI  = nan(nROI, nROI, nRun_cond);
    beta_phys = nan(nROI, nROI, nRun_cond);
    beta_psych= nan(nROI, nROI, nRun_cond);
    R2_mat    = nan(nROI, nROI, nRun_cond);

    for k = 1:nRun_cond
        X_conf = confounds.(condLabel).X{k};

        Y_all = struct();
        for r = 1:nROI
            rn = roiNames{r};
            Y_all.(rn) = roiSeed.(condLabel).(rn)(:, k);
        end

        for ii = 1:nROI
            rn_seed = roiNames{ii};
            phys    = Y_all.(rn_seed);

            phys_z = zscore(phys);

            for jj = 1:nROI
                if ii == jj
                    continue;
                end

                rn_targ = roiNames{jj};
                y       = Y_all.(rn_targ);

                psych = psych_mat;

                ppi = phys_z .* psych;

                Z = [ones(cfg.NTR_target,1), psych, phys_z, ppi];

                if cfg.use_confounds && ~isempty(X_conf)
                    N0 = X_conf;
                    N_ortho = N0 - Z * (Z \ N0);
                    X = [Z, N_ortho];
                else
                    X = Z;
                end

                y_filt = H * y;
                X_filt = H * X;

                beta = safe_glm_rankaware(y_filt, X_filt);

                yhat = X_filt * beta;
                res  = y_filt - yhat;
                R2   = 1 - var(res,0,1)/var(y_filt,0,1);

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

%% 5. 保存 run-wise 结果（用 save_subj 命名）
out_name = fullfile(cfg.out_dir, ...
    sprintf('%s_ppi_roi_roi_runwise_fromOriginal.mat', cfg.save_subj));
save(out_name, 'results', '-v7.3');
fprintf('\n>>> 已保存 run-wise ROI-ROI PPI 结果到: %s\n', out_name);

%% 6. 被试级固定效应（跨 run & 半球的双层 FFX：先 hemi，再 L/R 合并）
fprintf('\n[Step 6] 计算被试级固定效应 (跨 run + 左右半球合并)...\n');

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
results_fix.source_runwise_file = out_name;

for ci = 1:numel(condLevels)
    condLabel = condLevels{ci};

    if ~isfield(results, condLabel)
        warning('[FFX] results 中没有条件 %s，跳过', condLabel);
        continue;
    end

    beta_PPI  = results.(condLabel).beta_PPI;
    beta_phys = results.(condLabel).beta_phys;
    beta_psych= results.(condLabel).beta_psych;
    R2_mat    = results.(condLabel).R2;

    if isempty(beta_PPI)
        warning('[FFX] 条件 %s 的 beta_PPI 为空，跳过', condLabel);
        continue;
    end

    nRun_cond = size(beta_PPI,3);

    beta_PPI_fix_roi   = mean(beta_PPI,  3, 'omitnan');
    beta_phys_fix_roi  = mean(beta_phys, 3, 'omitnan');
    beta_psych_fix_roi = mean(beta_psych,3, 'omitnan');
    R2_fix_roi         = mean(R2_mat,    3, 'omitnan');

    Nruns_PPI_roi   = sum(~isnan(beta_PPI),  3);
    Nruns_phys_roi  = sum(~isnan(beta_phys), 3);
    Nruns_psych_roi = sum(~isnan(beta_psych),3);
    Nruns_R2_roi    = sum(~isnan(R2_mat),    3);

    beta_PPI_bilat = nan(nBase, nBase);
    Nr_bilat       = zeros(nBase, nBase);

    for b1 = 1:nBase
        iL = idxL(b1);
        iR = idxR(b1);
        if isnan(iL) || isnan(iR)
            continue;
        end
        for b2 = 1:nBase
            jL = idxL(b2);
            jR = idxR(b2);
            if isnan(jL) || isnan(jR)
                continue;
            end

            b_L_runs = squeeze(beta_PPI(iL, jL, :));
            b_L_runs = b_L_runs(~isnan(b_L_runs));
            if ~isempty(b_L_runs)
                beta_L_fix = mean(b_L_runs);
            else
                beta_L_fix = NaN;
            end

            b_R_runs = squeeze(beta_PPI(iR, jR, :));
            b_R_runs = b_R_runs(~isnan(b_R_runs));
            if ~isempty(b_R_runs)
                beta_R_fix = mean(b_R_runs);
            else
                beta_R_fix = NaN;
            end

            hemi_vals = [beta_L_fix, beta_R_fix];
            hemi_vals = hemi_vals(~isnan(hemi_vals));

            if isempty(hemi_vals)
                beta_PPI_bilat(b1,b2) = NaN;
                Nr_bilat(b1,b2)       = 0;
            else
                beta_PPI_bilat(b1,b2) = mean(hemi_vals);
                Nr_bilat(b1,b2)       = numel(b_L_runs) + numel(b_R_runs);
            end
        end
    end

    results_fix.(condLabel).beta_PPI_fix_roi   = beta_PPI_fix_roi;
    results_fix.(condLabel).beta_phys_fix_roi  = beta_phys_fix_roi;
    results_fix.(condLabel).beta_psych_fix_roi = beta_psych_fix_roi;
    results_fix.(condLabel).R2_fix_roi         = R2_fix_roi;

    results_fix.(condLabel).Nruns_PPI_roi      = Nruns_PPI_roi;
    results_fix.(condLabel).Nruns_phys_roi     = Nruns_phys_roi;
    results_fix.(condLabel).Nruns_psych_roi    = Nruns_psych_roi;
    results_fix.(condLabel).Nruns_R2_roi       = Nruns_R2_roi;

    results_fix.(condLabel).beta_PPI_bilat     = beta_PPI_bilat;
    results_fix.(condLabel).Nr_bilat           = Nr_bilat;
    results_fix.(condLabel).nRun_cond          = nRun_cond;

    fprintf('    [FFX] 条件 %-9s: nRun_cond = %d，已完成 ROI 级 & bilateral 级固定效应聚合。\n', ...
        condLabel, nRun_cond);
end

out_name_ffx = fullfile(cfg.out_dir, ...
    sprintf('%s_ppi_roi_roi_FFX_fromOriginal_bilateral.mat', cfg.save_subj));
save(out_name_ffx, 'results_fix', '-v7.3');

fprintf('\n>>> 已保存被试级固定效应（ROI + bilateral）结果到: %s\n', out_name_ffx);
fprintf('================ ROI-ROI PPI (FFX, within-subject, bilateral) 完成 ================\n');

%% 7. 可视化：基于 bilateral β 的单被试 ROI-ROI PPI 脑网络（含 FDR 校正）
fprintf('\n[Step 7] 可视化 bilateral 固定效应后的 ROI-ROI PPI 脑网络（含 FDR 校正）...\n');

alpha_unc = 0.05;
q_fdr     = 0.05;

col_red    = [1   0     0  ];
col_blue   = [0   0.6   1  ];
col_orange = [1   0.5   0  ];
col_green  = [0   1     0  ];
col_white  = [1   1     1  ];

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

    if ~isfield(results_fix.(condLabel), 'beta_PPI_bilat') || ...
       isempty(results_fix.(condLabel).beta_PPI_bilat)
        warning('[Vis] 条件 %s 的 beta_PPI_bilat 为空，跳过可视化。', condLabel);
        continue;
    end
    beta_fix_bilat = results_fix.(condLabel).beta_PPI_bilat;

    beta_run = [];
    if isfield(results, condLabel) && isfield(results.(condLabel), 'beta_PPI')
        beta_run = results.(condLabel).beta_PPI;
        if size(beta_run,1) ~= numel(roiNames) || size(beta_run,2) ~= numel(roiNames)
            warning('[Vis] 条件 %s 的 beta_PPI 尺寸与 roiNames 不符，跳过可视化。', condLabel);
            continue;
        end
    else
        warning('[Vis] 条件 %s 中没有 run-wise beta_PPI，p 值将全部为 NaN。', condLabel);
    end

    p_bilat_dir = NaN(nBase, nBase);

    if ~isempty(beta_run)
        for b1 = 1:nBase
            iL = idxL(b1); iR = idxR(b1);
            if isnan(iL) || isnan(iR), continue; end
            for b2 = 1:nBase
                if b1 == b2, continue; end
                jL = idxL(b2); jR = idxR(b2);
                if isnan(jL) || isnan(jR), continue; end

                b_L_runs = squeeze(beta_run(iL, jL, :));
                b_R_runs = squeeze(beta_run(iR, jR, :));

                b_all = [b_L_runs(:); b_R_runs(:)];
                b_all = b_all(~isnan(b_all));
                n_rep = numel(b_all);

                if n_rep < 2
                    continue;
                end

                m  = mean(b_all);
                sd = std(b_all, 0);
                if sd == 0
                    if m == 0, p_val = 1; else, p_val = eps; end
                else
                    t_val = m / (sd / sqrt(n_rep));
                    p_val = 2 * tcdf(-abs(t_val), n_rep - 1);
                end
                p_bilat_dir(b1,b2) = p_val;
            end
        end
    end

    beta_sym = NaN(nBase, nBase);
    p_sym    = NaN(nBase, nBase);

    for i = 1:nBase
        for j = i+1:nBase
            pair_beta  = [beta_fix_bilat(i,j), beta_fix_bilat(j,i)];
            valid_beta = pair_beta(~isnan(pair_beta));
            if ~isempty(valid_beta)
                b_ij = mean(valid_beta);
                beta_sym(i,j) = b_ij;
                beta_sym(j,i) = b_ij;
            end

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

    if all(isnan(beta_sym(:)))
        warning('[Vis] 条件 %s 所有 bilateral β 为 NaN，跳过绘图。', condLabel);
        continue;
    end

    mask_ut = triu(~isnan(p_sym), 1);
    p_vec   = p_sym(mask_ut);
    m_edge  = numel(p_vec);

    sig_fdr_mat = false(size(p_sym));
    p_crit = NaN;

    if m_edge > 0 && any(~isnan(p_vec))
        [h_fdr_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_mat(mask_ut) = h_fdr_vec;
        sig_fdr_mat = sig_fdr_mat | sig_fdr_mat.';
        fprintf('[Vis] 条件 %s：FDR q=%.3f，%d/%d 条边通过 FDR (p<=%.3g)\n', ...
            condLabel, q_fdr, sum(h_fdr_vec), m_edge, p_crit);
    else
        fprintf('[Vis] 条件 %s：没有可用于 FDR 的有效 p 值。\n', condLabel);
    end

    sig_unc_mat = (p_sym < alpha_unc) & ~isnan(p_sym);
    sig_unc_mat = sig_unc_mat | sig_unc_mat.';

    theta = linspace(0, 2*pi, nBase+1);
    theta = theta(1:nBase);
    R     = 1;
    x = R * cos(theta);
    y = R * sin(theta);

    valid_abs_beta = abs(beta_sym(~isnan(beta_sym) & (triu(true(size(beta_sym)),1))));
    max_abs_beta   = max(valid_abs_beta);
    if isempty(max_abs_beta) || max_abs_beta == 0
        max_abs_beta = 1;
    end

    fig_name = sprintf('%s - %s PPI network (FFX bilateral, Psych_diff, FDR)', ...
        cfg.save_subj, condLabel);
    figure('Name', fig_name);
    set(gcf, 'Color', [0 0 0]);
    ax = gca;
    set(ax, 'Color', [0 0 0]);
    hold on; axis equal off;

    for i = 1:nBase-1
        for j = i+1:nBase
            if isnan(beta_sym(i,j)), continue; end

            b_ij = beta_sym(i,j);
            is_fdr = sig_fdr_mat(i,j);
            is_unc = sig_unc_mat(i,j);

            edge_color = col_white;
            if is_fdr
                if b_ij > 0, edge_color = col_red;
                elseif b_ij < 0, edge_color = col_blue;
                end
            elseif is_unc
                if b_ij > 0, edge_color = col_orange;
                elseif b_ij < 0, edge_color = col_green;
                end
            end

            w = 1.5 + 3 * (abs(b_ij) / max_abs_beta);

            line([x(i) x(j)], [y(i) y(j)], ...
                'Color', edge_color, 'LineWidth', w);

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

    title(sprintf('%s - %s: PPI network (FFX bilateral, numerosity vs baseline, FDR)', ...
        cfg.save_subj, condLabel), ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('\n[Step 7] bilateral 固定效应后的脑网络可视化完成。\n');

%% ========================= 子函数 =========================

function [sesIdx, runInSes] = scanNum_to_sesrun_sub03style(scanNum)
% Scan1..8   -> ses-01, run1..8
% Scan9..17  -> ses-02, run1..9
% Scan18..27 -> ses-03, run1..10
% Scan28..36 -> ses-04, run1..9
sesIdx   = NaN;
runInSes = NaN;

if scanNum >= 1 && scanNum <= 8
    sesIdx   = 1;
    runInSes = scanNum;
elseif scanNum >= 9 && scanNum <= 17
    sesIdx   = 2;
    runInSes = scanNum - 8;
elseif scanNum >= 18 && scanNum <= 27
    sesIdx   = 3;
    runInSes = scanNum - 17;
elseif scanNum >= 28 && scanNum <= 36
    sesIdx   = 4;
    runInSes = scanNum - 27;
end
end

function root = get_func_root_for_session(sesLabel, cfg)
% 通用规则：resultsX 对应 ses-0X（你的 sub-03/sub-02 都符合这个规律）
sesNum = sscanf(sesLabel, 'ses-%d');
assert(~isempty(sesNum), '无法从 sesLabel 解析数字: %s', sesLabel);

root = fullfile(cfg.root_base, sprintf('results%d', sesNum));
assert(exist(root,'dir')==7, 'func 根目录不存在: %s', root);
end

function [wm_nii, csf_nii] = resolve_wm_csf_probseg(cfg)
% 优先用 cfg.t1_ses；若不存在则在 ses-01..ses-04 里自动寻找
ses_try = {};
ses_try{end+1} = cfg.t1_ses;

% 加入其它候选（去重）
for s = 1:4
    ses_try{end+1} = sprintf('ses-%02d', s);
end
ses_try = unique(ses_try, 'stable');

wm_nii  = '';
csf_nii = '';

for i = 1:numel(ses_try)
    ses = ses_try{i};
    sesNum = sscanf(ses, 'ses-%d');
    if isempty(sesNum), continue; end

    root_i = fullfile(cfg.root_base, sprintf('results%d', sesNum));
    anatDir = fullfile(root_i, cfg.data_subj, ses, 'anat');

    wm1  = fullfile(anatDir, sprintf('%s_%s_label-WM_probseg.nii',  cfg.data_subj, ses));
    wm2  = [wm1 '.gz'];
    csf1 = fullfile(anatDir, sprintf('%s_%s_label-CSF_probseg.nii', cfg.data_subj, ses));
    csf2 = [csf1 '.gz'];

    if isempty(wm_nii)
        if exist(wm1,'file')==2, wm_nii = wm1;
        elseif exist(wm2,'file')==2, wm_nii = wm2;
        end
    end
    if isempty(csf_nii)
        if exist(csf1,'file')==2, csf_nii = csf1;
        elseif exist(csf2,'file')==2, csf_nii = csf2;
        end
    end

    if ~isempty(wm_nii) && ~isempty(csf_nii)
        fprintf('    [Info] 使用 T1 概率图: %s (WM), %s (CSF)\n', wm_nii, csf_nii);
        return;
    end
end

assert(~isempty(wm_nii), '找不到任何 session 的 WM_probseg（检查 results1..4 / anat）');
assert(~isempty(csf_nii), '找不到任何 session 的 CSF_probseg（检查 results1..4 / anat）');
end

function [conf_tsv, bold_nii, baseStub_used] = resolve_fmriprep_run_files(subj, sesLabel, runInSes, sesDir)
% 兼容 run 编号可能是 2 位(run02) 或 3 位(run010) 的情况
% 也兼容 task-run02 / task-run-02 两种写法

runFmtList = {sprintf('%02d', runInSes), sprintf('%03d', runInSes)};
stubList = {};
for i = 1:numel(runFmtList)
    rstr = runFmtList{i};
    stubList{end+1} = sprintf('%s_%s_task-run%s',  subj, sesLabel, rstr);  %#ok<AGROW>
    stubList{end+1} = sprintf('%s_%s_task-run-%s', subj, sesLabel, rstr);  %#ok<AGROW>
end

conf_tsv = '';
bold_nii = '';
baseStub_used = '';

tried = {};

for s = 1:numel(stubList)
    baseStub = stubList{s};

    c = fullfile(sesDir, [baseStub '_desc-confounds_timeseries.tsv']);
    tried{end+1} = c; %#ok<AGROW>
    if exist(c,'file')==2
        b1 = fullfile(sesDir, [baseStub '_space-T1w_desc-preproc_bold.nii']);
        b2 = [b1 '.gz'];
        tried{end+1} = b1; %#ok<AGROW>
        tried{end+1} = b2; %#ok<AGROW>

        if exist(b1,'file')==2
            conf_tsv = c; bold_nii = b1; baseStub_used = baseStub; return;
        elseif exist(b2,'file')==2
            conf_tsv = c; bold_nii = b2; baseStub_used = baseStub; return;
        end
    end
end

msg = sprintf('找不到 confounds/BOLD（已尝试 2位/3位 run 编号 & task-run / task-run- 形式）。\n示例尝试路径(部分)：\n');
nShow = min(numel(tried), 12);
for i = 1:nShow
    msg = sprintf('%s  - %s\n', msg, tried{i});
end
error('%s', msg);
end

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

function [Xc, cnames, meta_out] = build_confounds_one_run(conf_tsv, bold_nii, wm_nii, csf_nii, cfg)
% 为单个 run 构建混杂项矩阵 Xc （motion+dmotion+FD_pulse+aCompCor）
% 现在：aCompCor 直接使用 fMRIPrep tsv 中已有的 w_comp_cor_* / c_comp_cor_*
% bold_nii, wm_nii, csf_nii 不再参与 aCompCor 计算，仅保留接口一致性

Tc = readtable(conf_tsv, 'FileType','text', 'Delimiter','\t');
N_full = height(Tc);

if N_full == cfg.NTR_target
    idx = (1:cfg.NTR_target)';
elseif N_full >= cfg.n_dummy + cfg.NTR_target
    idx = (cfg.n_dummy+1) : (cfg.n_dummy + cfg.NTR_target);
else
    error('confounds 行数(%d) 太少，无法通过 dummy(%d) 获得 %d 帧', ...
        N_full, cfg.n_dummy, cfg.NTR_target);
end

Tc_use = Tc(idx, :);

motNames = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
hasAll = all(ismember(motNames, Tc_use.Properties.VariableNames));
assert(hasAll, 'confounds 中缺少运动列 (trans/rot)');

mot = table2array(Tc_use(:, motNames));
dmot = [zeros(1,6); diff(mot,1,1)];

if ismember('framewise_displacement', Tc_use.Properties.VariableNames)
    FD = Tc_use.framewise_displacement;
    FD(isnan(FD)) = 0;
    FD_pulse = double(FD > cfg.fd_thr);
else
    warning('confounds 中没有 framewise_displacement 列，FD 脉冲列用 0 填充');
    FD_pulse = zeros(cfg.NTR_target,1);
end

allVars = Tc_use.Properties.VariableNames;

wmColsAll = allVars(startsWith(allVars, 'w_comp_cor_'));
wmColsAll = sort(wmColsAll);

csfColsAll = allVars(startsWith(allVars, 'c_comp_cor_'));
csfColsAll = sort(csfColsAll);

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

Xc_raw = [mot, dmot, FD_pulse, score_wm, score_csf];

cnames_raw = [ ...
    motNames, ...
    strcat('d_', motNames), ...
    {'FD_pulse'}, ...
    wmCols, ...
    csfCols ...
];

Xc_raw = Xc_raw - mean(Xc_raw,1,'omitnan');

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

meta_out = struct();
meta_out.N_full     = N_full;
meta_out.idx        = idx;
meta_out.motNames   = motNames;
meta_out.fd_thr     = cfg.fd_thr;
meta_out.wmCols     = wmCols;
meta_out.csfCols    = csfCols;
meta_out.source     = 'fMRIPrep_confounds_tsv';
meta_out.constColsDropped = cnames_raw(~keep);
end

function beta_full = safe_glm_rankaware(y, X)
% 安全求解 beta，自动处理秩亏（QR + 列主元）
[Q,R,perm] = qr(X,0);
d   = abs(diag(R));
if isempty(d)
    beta_full = zeros(size(X,2),1);
    return;
end

tol = max(size(X)) * eps(max(d));
r   = sum(d > tol);

if r == 0
    warning('safe_glm_rankaware: design matrix rank = 0，返回全 0 beta。');
    beta_full = zeros(size(X,2),1);
    return;
end

if r < size(X,2)
    warning('safe_glm_rankaware: design matrix rank deficient: P = %d, rank = %d, drop %d col(s).',...
            size(X,2), r, size(X,2)-r);
end

Q1 = Q(:,1:r);
R1 = R(1:r,1:r);
beta_red = R1 \ (Q1' * y);

beta_full = zeros(size(X,2),1);
beta_full(perm(1:r)) = beta_red;
end

function [h, p_crit] = fdr_bh_local(p_vals, q)
% Benjamini-Hochberg (local)
p = p_vals(:);
m = numel(p);
h = false(size(p));

valid = ~isnan(p);
pv    = p(valid);
if isempty(pv)
    p_crit = NaN;
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
