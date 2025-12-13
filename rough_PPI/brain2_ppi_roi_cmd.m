% 2brain_ppi_roi_cmd.m
% 仅使用 ROI-ROI 的 PPI 管线骨架（6 个 ROI 两两做 PPI）
% ----------------------------------------------------------
clear; clc;
fprintf('\n================ ROI-ROI PPI Pipeline (skeleton) ================\n');

%% 0. 配置参数（你只要改这里的路径/名字）
cfg = struct();

% 基本信息
cfg.subj        = 'sub-10';            % TODO: 被试 ID
cfg.ses       = 'ses-01'; 

cfg.root_bids   = '/home/xue/data/BIDS_prep';     % 仅用于输出路径组织

% 事件文件（所有 run 共用）
cfg.events_tsv  = '/media/xue/new_B1/glm/events_output/events_bids.tsv';

% 时间参数
cfg.TR          = 1.95;      % 你的真实 TR
cfg.n_dummy     = 6;         % 记录在案，不在本脚本中裁剪
cfg.NTR_target  = 176;       % 你希望 ROI 时序长度必须是 176

% 心理条件名（trial_type 中要包含，约定：1= numerosity, 2= baseline）
cfg.cond_names  = {'numerosity', 'baseline'};

% ROI blackbox .mat —— 两个 dataType 条件
cfg.roi_mat_files = struct();
base_dir = '/home/xue/data/prf_result';
cfg.roi_mat_files.unconnect = ...
    fullfile(base_dir, cfg.subj, cfg.ses, ...
        'results_prf_roi_tseries_cmd', ...
        [cfg.subj '_Averages_unconnect_layer_01_roiSeed.mat']);

cfg.roi_mat_files.connect = ...
    fullfile(base_dir, cfg.subj, cfg.ses, ...
        'results_prf_roi_tseries_cmd', ...
        [cfg.subj '_Averages_connect_layer_01_roiSeed.mat']);



% 输出目录
cfg.out_dir = fullfile(cfg.root_bids, 'derivatives', 'ppi_roi', cfg.subj);
if ~exist(cfg.out_dir, 'dir')
    fprintf('  [Info] 创建输出目录: %s\n', cfg.out_dir);
    mkdir(cfg.out_dir);
end

fprintf('  [Config] subj=%s\n', cfg.subj);
fprintf('  [Config] TR=%.3f, NTR_target=%d\n', cfg.TR, cfg.NTR_target);
fprintf('  [Config] cond_names: %s\n', strjoin(cfg.cond_names, ', '));

%% 1. 读取 events.tsv 并构造 Psych_diff (numerosity vs baseline)
fprintf('\n[Step 1] 从 events_bids.tsv 构造心理量 Psych_diff (numerosity vs baseline)...\n');

assert(exist(cfg.events_tsv, 'file') == 2, ...
    'events.tsv 不存在: %s', cfg.events_tsv);

Te = readtable(cfg.events_tsv, 'FileType', 'text', 'Delimiter', '\t');
fprintf('    events.tsv: %d 行, %d 列\n', height(Te), width(Te));

% 必要列：onset, duration, trial_type
req_cols = {'onset','duration','trial_type'};
missing_cols = setdiff(req_cols, Te.Properties.VariableNames);
assert(isempty(missing_cols), 'events.tsv 缺少列: %s', strjoin(missing_cols, ', '));

% 时间轴：按目标 NTR 构造
t = (0:cfg.NTR_target-1)' * cfg.TR;   % [NTR_target x 1]

% HRF：优先用 SPM 的 canonical HRF，没有的话用一个简单 gamma 近似
use_spm_hrf = (exist('spm_hrf', 'file') == 2);
if use_spm_hrf
    fprintf('    使用 spm_hrf 生成 canonical HRF.\n');
    hrf = spm_hrf(cfg.TR);
else
    fprintf('    找不到 spm_hrf，使用简易 gamma HRF 近似.\n');
    tt = (0:1:32)';
    hrf = tt.^8 .* exp(-tt/0.9);
    hrf = hrf / sum(hrf);
    hrf = interp1(tt, hrf, 0:cfg.TR:32, 'pchip', 0)';
end

% --- numerosity / baseline boxcar ---
assert(numel(cfg.cond_names) == 2, ...
    '当前 Psych_diff 假定 cond_names 只有两个：numerosity & baseline');

cname_num  = cfg.cond_names{1};
cname_base = cfg.cond_names{2};

% numerosity
mask_num = strcmp(Te.trial_type, cname_num);
assert(any(mask_num), 'events.tsv 中找不到 trial_type="%s" 的事件', cname_num);
onset_num = Te.onset(mask_num);
dur_num   = Te.duration(mask_num);
fprintf('    条件 %-10s: 事件数 = %d\n', cname_num, numel(onset_num));

box_num = zeros(cfg.NTR_target, 1);
for k = 1:numel(onset_num)
    t_on  = onset_num(k);
    t_off = onset_num(k) + dur_num(k);
    box_num = box_num + (t >= t_on & t < t_off);
end

% baseline
mask_base = strcmp(Te.trial_type, cname_base);
assert(any(mask_base), 'events.tsv 中找不到 trial_type="%s" 的事件', cname_base);
onset_base = Te.onset(mask_base);
dur_base   = Te.duration(mask_base);
fprintf('    条件 %-10s: 事件数 = %d\n', cname_base, numel(onset_base));

box_base = zeros(cfg.NTR_target, 1);
for k = 1:numel(onset_base)
    t_on  = onset_base(k);
    t_off = onset_base(k) + dur_base(k);
    box_base = box_base + (t >= t_on & t < t_off);
end

% 差异 boxcar: numerosity(+1) vs baseline(-1)
diff_box = box_num - box_base;

% HRF 卷积 + 截断
conv_diff = conv(diff_box, hrf);
conv_diff = conv_diff(1:cfg.NTR_target);

% z-score 处理得到 Psych_diff
conv_diff = conv_diff - mean(conv_diff, 'omitnan');
sd_diff   = std(conv_diff, 0, 'omitnan');
if sd_diff == 0
    warning('    Psych_diff 方差为 0，置为零列。');
    conv_diff(:) = 0;
else
    conv_diff = conv_diff / sd_diff;
end

psych_mat   = conv_diff;          % [NTR x 1]
psych_names = {'Psych_diff'};

fprintf('    已生成差异心理回归量 Psych_diff（行数=%d）。\n', size(psych_mat,1));

%% 2. 对每个 dataType 条件做 ROI-ROI PPI（L/R 仍单独计算）
fprintf('\n[Step 2] 对每个 dataType 条件进行 ROI-ROI PPI（原始 L/R ROI 层面）...\n');

cond_levels = fieldnames(cfg.roi_mat_files);
results = struct();

for ci = 1:numel(cond_levels)
    cond_label = cond_levels{ci};
    roi_mat    = cfg.roi_mat_files.(cond_label);

    fprintf('\n  === DataType 条件: %s ===\n', cond_label);
    fprintf('      ROI ts .mat: %s\n', roi_mat);
    assert(exist(roi_mat, 'file') == 2, ...
        'ROI ts .mat 不存在: %s', roi_mat);

    S = load(roi_mat);
    assert(isfield(S, 'roiSeed') && isfield(S, 'meta'), ...
        'ROI .mat 中缺少 roiSeed 或 meta');

    % meta.NTR 与 cfg.NTR_target 对齐检查
    if isfield(S, 'meta') && isfield(S.meta, 'NTR')
        NTR_meta = S.meta.NTR;
        fprintf('      meta.NTR = %d\n', NTR_meta);
        assert(NTR_meta == cfg.NTR_target, ...
            'meta.NTR (%d) != cfg.NTR_target (%d)', NTR_meta, cfg.NTR_target);
    else
        fprintf('      [Warning] meta 中没有 NTR 字段，仅用 ROI 长度做检查。\n');
    end

    % ROI 名列表
    if isfield(S.meta, 'roiNames')
        roi_all = S.meta.roiNames(:);
    else
        roi_all = fieldnames(S.roiSeed);
        fprintf('      [Warning] meta.roiNames 不存在，用 roiSeed 字段名代替。\n');
    end

    fprintf('      .mat 中包含 ROI: %s\n', strjoin(roi_all', ', '));
    nROI = numel(roi_all);
    fprintf('      本条件下参与 PPI 的 ROI 数: %d\n', nROI);
    fprintf('      共有 %d 对有序 (seed,target) 组合（不含自连接）。\n', nROI*(nROI-1));

    % 预先对 Psych_diff 再做一次 z-score（稳一点）
    n_psych = size(psych_mat, 2);   % 这里应为 1
    psych_z = zeros(size(psych_mat));
    for c = 1:n_psych
        pc = psych_mat(:, c);
        pc = pc - mean(pc, 'omitnan');
        sd_pc = std(pc, 0, 'omitnan');
        if sd_pc == 0
            warning('      Psych 列 %s 方差为 0，置为零列。', psych_names{c});
            pc(:) = 0;
        else
            pc = pc / sd_pc;
        end
        psych_z(:, c) = pc;
    end

    res_cond = struct();

    % ===== 外层循环：每个 ROI 轮流当 seed =====
    for si = 1:nROI
        seed_name = roi_all{si};

        assert(isfield(S.roiSeed, seed_name), ...
            'seed ROI "%s" 不在 roiSeed 中', seed_name);

        seed = S.roiSeed.(seed_name);
        seed = seed(:);
        assert(numel(seed) == cfg.NTR_target, ...
            'Seed ROI %s 时序长度 (%d) != cfg.NTR_target (%d)', ...
            seed_name, numel(seed), cfg.NTR_target);

        % z-score seed
        seed = seed - mean(seed, 'omitnan');
        sd_seed = std(seed, 0, 'omitnan');
        assert(sd_seed > 0, 'Seed ROI %s 的时序方差为 0', seed_name);
        seed_z = seed / sd_seed;

        % 构造 PPI_diff（对当前 seed）
        ppi_raw = seed_z .* psych_z(:,1);
        ppi_raw = ppi_raw - mean(ppi_raw, 'omitnan');
        sd_ppi  = std(ppi_raw, 0, 'omitnan');
        if sd_ppi == 0
            warning('      [Seed %s] PPI_diff 方差为 0，置为零列。', seed_name);
            ppi_z = zeros(cfg.NTR_target, 1);
        else
            ppi_z = ppi_raw / sd_ppi;
        end

        % 设计矩阵 X： [const, Seed, Psych_diff, PPI_diff]
        X = [ones(cfg.NTR_target, 1), seed_z, psych_z, ppi_z];
        col_names = {'const', ...
                     sprintf('Phys_%s', seed_name), ...
                     'Psych_diff', ...
                     sprintf('PPI_Psych_diff_x_%s', seed_name)};

        fprintf('      [Seed %s] 设计矩阵尺寸: [%d 时间点 x %d 回归量]\n', ...
            seed_name, size(X,1), size(X,2));

        % 秩与共线性检查（对当前 seed 一次）
        rX = rank(X);
        if rX < size(X,2)
            warning('      [Seed %s] 设计矩阵非满秩: rank=%d, nCols=%d', ...
                seed_name, rX, size(X,2));
        else
            fprintf('      [Seed %s] 设计矩阵满秩: rank=%d。\n', seed_name, rX);
        end

        R = corrcoef(X);
        R_no_diag = R - diag(diag(R));
        max_corr = max(abs(R_no_diag(:)));
        fprintf('      [Seed %s] 最大 |相关系数| = %.3f\n', seed_name, max_corr);
        if max_corr > 0.9
            warning('      [Seed %s] 存在高度共线回归量 (>|0.90|)，建议检查。', seed_name);
        end

        % 为了效率，XtX 只求一次
        XtX_inv = inv(X' * X);

        % 为当前 seed 创建子结构
        res_cond.(seed_name) = struct();

        fprintf('      [Seed %s] 开始对所有 target ROI 拟合 GLM...\n', seed_name);

        % ===== 内层循环：每个其他 ROI 当 target =====
        for ti = 1:nROI
            if ti == si
                continue;   % 不做 self-PPI
            end

            target_name = roi_all{ti};
            fprintf('        [Seed %s -> Target %s] ... ', seed_name, target_name);

            assert(isfield(S.roiSeed, target_name), ...
                'target ROI "%s" 不在 roiSeed 中', target_name);

            y = S.roiSeed.(target_name);
            y = y(:);
            assert(numel(y) == cfg.NTR_target, ...
                'ROI %s 时序长度 (%d) != cfg.NTR_target (%d)', ...
                target_name, numel(y), cfg.NTR_target);

            % GLM：y = X * beta + eps
            beta = X \ y;
            y_hat = X * beta;
            resid = y - y_hat;

            df = cfg.NTR_target - rX;
            assert(df > 0, '[Seed %s -> Target %s] 自由度 df<=0，检查设计矩阵。', ...
                seed_name, target_name);

            s2   = sum(resid.^2) / df;
            se_b = sqrt(diag(s2 * XtX_inv));
            t_b  = beta ./ se_b;
            p_b  = 2 * tcdf(-abs(t_b), df);  % 双侧检验

            % PPI_diff 在最后一列
            idx_ppi = numel(col_names);
            beta_ppi = beta(idx_ppi);
            t_ppi    = t_b(idx_ppi);
            p_ppi    = p_b(idx_ppi);

            % 打印结果
            fprintf('PPI_diff: beta=%.4f, t=%.2f, p=%.3g\n', ...
                beta_ppi, t_ppi, p_ppi);

            % 保存结果
            rr = struct();
            rr.seed      = seed_name;
            rr.target    = target_name;
            rr.beta      = beta;
            rr.se_beta   = se_b;
            rr.t_beta    = t_b;
            rr.p_beta    = p_b;
            rr.beta_ppi  = beta_ppi;
            rr.t_ppi     = t_ppi;
            rr.p_ppi     = p_ppi;
            rr.df        = df;
            rr.col_names = col_names;

            res_cond.(seed_name).(target_name) = rr;
        end
    end

    % 把本 dataType 条件的结果写入总结构
    results.(cond_label) = res_cond;

    % （保持原有保存逻辑不改，如果需要可以额外再保存 bilateral 结果）
    out_name_cond = fullfile(cfg.out_dir, ...
        sprintf('%s_ppi_roi_roi_%s_network.mat', cfg.subj, cond_label));
    save(out_name_cond, 'cfg', 'cond_label', 'psych_mat', 'results', '-v7.3');

    fprintf('      [保存] 本条件结果已写入: %s\n', out_name_cond);
end

%% 2.5 在被试内合并左右半脑 —— 得到 bilateral ROI-ROI PPI 结果
fprintf('\n[Step 2.5] 在被试内对 L/R 网络做固定效应合并，得到 bilateral ROI-ROI PPI 结果...\n');

results_bilat = struct();
cond_levels = fieldnames(results);

for ci = 1:numel(cond_levels)
    cond_label = cond_levels{ci};
    if ~isfield(results, cond_label)
        continue;
    end

    res_cond = results.(cond_label);
    roi_all  = fieldnames(res_cond);

    % 1) 解析 ROI 名，找出同时有 L/R 的 base ROI
    base_info = struct();
    for ri = 1:numel(roi_all)
        roi_name = roi_all{ri};
        [base_name, hemi] = parse_roi_name(roi_name);
        if isempty(hemi)
            % 非 L/R 开头的 ROI 不参与 bilateral 合并
            continue;
        end
        if ~isfield(base_info, base_name)
            tmp = struct();
            tmp.hemi = {};
            base_info.(base_name) = tmp;
        end
        if ~ismember(hemi, base_info.(base_name).hemi)
            base_info.(base_name).hemi{end+1} = hemi;
        end
    end

    if isempty(fieldnames(base_info))
        warning('  [HemiCombine] 条件 %s 中没有可识别的 L/R ROI，跳过。', cond_label);
        continue;
    end

    base_names_all = fieldnames(base_info);
    bilat_nodes = {};
    for bi = 1:numel(base_names_all)
        bname = base_names_all{bi};
        hemis = base_info.(bname).hemi;
        if ismember('L', hemis) && ismember('R', hemis)
            bilat_nodes{end+1} = bname; %#ok<AGROW>
        end
    end

    if isempty(bilat_nodes)
        warning('  [HemiCombine] 条件 %s 中没有同时具有 L/R 的 ROI，跳过。', cond_label);
        continue;
    end

    fprintf('  [HemiCombine] 条件 %s：%d 个 bilateral 节点: %s\n', ...
        cond_label, numel(bilat_nodes), strjoin(bilat_nodes, ', '));

    % 2) 针对每对 base ROI，合并 L/R 的 seed->target 边
    res_bilat_cond = struct();
    nBase = numel(bilat_nodes);

    for si = 1:nBase
        seed_base = bilat_nodes{si};
        res_bilat_cond.(seed_base) = struct();  % 初始化 seed 结构

        for ti = 1:nBase
            if ti == si
                continue;
            end
            target_base = bilat_nodes{ti};

            beta_vec = [];
            p_vec    = [];

            % 左半球 L<seed> -> L<target>
            seed_L   = ['L' seed_base];
            target_L = ['L' target_base];
            if isfield(res_cond, seed_L) && isfield(res_cond.(seed_L), target_L)
                rrL = res_cond.(seed_L).(target_L);
                if isfield(rrL, 'beta_ppi') && isfield(rrL, 'p_ppi') && ~isnan(rrL.beta_ppi)
                    beta_vec(end+1) = rrL.beta_ppi; %#ok<AGROW>
                    p_vec(end+1)    = rrL.p_ppi;   %#ok<AGROW>
                end
            end

            % 右半球 R<seed> -> R<target>
            seed_R   = ['R' seed_base];
            target_R = ['R' target_base];
            if isfield(res_cond, seed_R) && isfield(res_cond.(seed_R), target_R)
                rrR = res_cond.(seed_R).(target_R);
                if isfield(rrR, 'beta_ppi') && isfield(rrR, 'p_ppi') && ~isnan(rrR.beta_ppi)
                    beta_vec(end+1) = rrR.beta_ppi; %#ok<AGROW>
                    p_vec(end+1)    = rrR.p_ppi;   %#ok<AGROW>
                end
            end

            if isempty(beta_vec)
                % 该 base 边在任一半球都没成功拟合，跳过
                continue;
            end

            % 固定效应合并：β 取均值；p 取最小值（可视化用）
            beta_bilat = mean(beta_vec, 'omitnan');
            p_bilat    = min(p_vec);

            rr_bilat = struct();
            rr_bilat.seed     = seed_base;
            rr_bilat.target   = target_base;
            rr_bilat.beta_ppi = beta_bilat;
            rr_bilat.p_ppi    = p_bilat;

            res_bilat_cond.(seed_base).(target_base) = rr_bilat;
        end
    end

    results_bilat.(cond_label) = res_bilat_cond;
end

% 如果需要，可以额外保存 bilateral 结果
out_name_bilat = fullfile(cfg.out_dir, ...
    sprintf('%s_ppi_roi_roi_bilateral_allconds.mat', cfg.subj));
save(out_name_bilat, 'cfg', 'psych_mat', 'results', 'results_bilat', '-v7.3');
fprintf('  [HemiCombine] bilateral 结果已保存到: %s\n', out_name_bilat);

%% 3. 简单脑网络可视化（基于 bilateral 结果 + 双向 PPI 的对称 β + FDR 校正）
% 说明：
%   - 使用上面 Step 2.5 得到的 results_bilat 结构；
%   - 对每个 dataType 条件（unconnect / connect）分别画一张图；
%   - 节点：bilateral base ROI（去掉 L/R 前缀）；
%   - 边的 β：seed->target 和 target->seed 的 β_ppi 取平均；
%   - 边的 p：取双向 p_ppi 的最小值；
%   - 线颜色逻辑同原脚本。

fprintf('\n[Step 3] 可视化 bilateral ROI-ROI PPI 脑网络（unconnect / connect，含 FDR 校正）...\n');

alpha_unc = 0.05;   % 未校正显著性阈值
q_fdr     = 0.05;   % FDR 控制水准（Benjamini-Hochberg）

% 颜色定义
col_red    = [1   0     0  ];  % FDR 正向显著
col_blue   = [0   0.6   1  ];  % FDR 负向显著
col_orange = [1   0.5   0  ];  % 未校正正向显著但 FDR 不显著
col_green  = [0   1     0  ];  % 未校正负向显著但 FDR 不显著
col_white  = [1   1     1  ];

% 这里的 cond_levels 来自 bilateral 结果
cond_levels = fieldnames(results_bilat);

for ci = 1:numel(cond_levels)
    cond_label = cond_levels{ci};
    if ~isfield(results_bilat, cond_label)
        warning('  [Warning] bilateral 结果中不存在条件 %s，跳过可视化。', cond_label);
        continue;
    end

    % 取该条件下的所有 ROI 名（作为 seed 的字段名）
    res_cond = results_bilat.(cond_label);
    roi_all  = fieldnames(res_cond);
    nROI     = numel(roi_all);

    if nROI < 2
        warning('  [Warning] 条件 %s 中 bilateral ROI 数 < 2，无法画网络。', cond_label);
        continue;
    end

    fprintf('  [Network] 条件 %s：共 %d 个 bilateral ROI，构建对称 β / p 矩阵...\n', ...
        cond_label, nROI);

    % 建立有向矩阵（seed->target）
    beta_dir = NaN(nROI, nROI);
    p_dir    = NaN(nROI, nROI);

    for si = 1:nROI
        seed_name = roi_all{si};
        res_seed  = res_cond.(seed_name);

        targ_names = fieldnames(res_seed);
        for ti = 1:numel(targ_names)
            target_name = targ_names{ti};
            idx_t = find(strcmp(roi_all, target_name));
            if isempty(idx_t)
                continue;
            end

            beta_dir(si, idx_t) = res_seed.(target_name).beta_ppi;
            p_dir(si, idx_t)    = res_seed.(target_name).p_ppi;
        end
    end

    % 构建对称 β / p 矩阵（只用于可视化，不改变统计本身）
    beta_sym = NaN(nROI, nROI);
    p_sym    = NaN(nROI, nROI);

    for i = 1:nROI
        for j = i+1:nROI
            pair_beta = [beta_dir(i,j), beta_dir(j,i)];
            pair_p    = [p_dir(i,j),    p_dir(j,i)];

            % β：取非 NaN 的平均值
            valid_beta = pair_beta(~isnan(pair_beta));
            if ~isempty(valid_beta)
                beta_sym(i,j) = mean(valid_beta);
                beta_sym(j,i) = beta_sym(i,j);
            end

            % p：取非 NaN 的最小值
            valid_p = pair_p(~isnan(pair_p));
            if ~isempty(valid_p)
                p_sym(i,j) = min(valid_p);
                p_sym(j,i) = p_sym(i,j);
            end
        end
        beta_sym(i,i) = 0;
        p_sym(i,i)    = NaN;
    end

    % 打印一个 β 矩阵的 summary（对称矩阵的一半）
    fprintf('  [Network] 条件 %s：对称 β 矩阵（上三角）\n', cond_label);
    for i = 1:nROI
        for j = i+1:nROI
            if isnan(beta_sym(i,j)), continue; end
            fprintf('    %s - %s : beta=%.4f, p=%.3g\n', ...
                roi_all{i}, roi_all{j}, beta_sym(i,j), p_sym(i,j));
        end
    end

    % 若全部 NaN，则不画图
    if all(isnan(beta_sym(:)))
        warning('  [Warning] 条件 %s 所有 β 为 NaN，跳过绘图。', cond_label);
        continue;
    end

    % ---------- 计算 FDR 校正（对每个条件单独做） ----------
    % 只对上三角（i<j）非 NaN 的 p 做 FDR
    mask_ut  = triu(~isnan(p_sym), 1);
    p_vec    = p_sym(mask_ut);
    m        = numel(p_vec);

    sig_fdr_mat = false(size(p_sym));
    p_crit = NaN;

    if m > 0
        [h_fdr_vec, p_crit] = fdr_bh_local(p_vec, q_fdr);
        sig_fdr_mat(mask_ut) = h_fdr_vec;
        sig_fdr_mat = sig_fdr_mat | sig_fdr_mat.';  % 对称化
        fprintf('  [Network] 条件 %s：FDR q=%.3f，%d/%d 条边通过 FDR (p<=%.3g)\n', ...
            cond_label, q_fdr, sum(h_fdr_vec), m, p_crit);
    else
        fprintf('  [Network] 条件 %s：没有可用于 FDR 的边。\n', cond_label);
    end

    % 未校正显著性矩阵（两侧 p<alpha_unc）
    sig_unc_mat = (p_sym < alpha_unc) & ~isnan(p_sym);
    sig_unc_mat = sig_unc_mat | sig_unc_mat.';

    % 节点坐标：均匀放在一个圆上
    theta = linspace(0, 2*pi, nROI+1);
    theta = theta(1:nROI);
    R     = 1;  % 半径

    x = R * cos(theta);
    y = R * sin(theta);

    % 计算最大 |beta| 用于设置线宽
    valid_abs_beta = abs(beta_sym(~isnan(beta_sym)));
    max_abs_beta   = max(valid_abs_beta);
    if isempty(max_abs_beta) || max_abs_beta == 0
        max_abs_beta = 1;
    end

    % 画图
    figure('Name', sprintf('%s - %s PPI network (bilateral, Psych_diff, FDR)', cfg.subj, cond_label));
    set(gcf, 'Color', [0 0 0]);      % 黑色背景
    ax = gca;
    set(ax, 'Color', [0 0 0]);
    hold on; axis equal off;

    % 画所有边（只画上三角）
    for i = 1:nROI-1
        for j = i+1:nROI
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
                    edge_color = col_red;    % FDR 正向显著
                elseif b_ij < 0
                    edge_color = col_blue;   % FDR 负向显著
                else
                    edge_color = col_white;
                end
            elseif is_unc
                % 未经校正显著但 FDR 不显著
                if b_ij > 0
                    edge_color = col_orange; % 正向：橙色
                elseif b_ij < 0
                    edge_color = col_green;  % 负向：绿色
                else
                    edge_color = col_white;
                end
            else
                % 未经校正就不显著：保持白色
                edge_color = col_white;
            end

            % 线宽随 |beta| 变化
            w = 1.5 + 3 * (abs(b_ij) / max_abs_beta);

            % 画边
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

    % 再画节点
    node_size = 300;
    for i = 1:nROI
        scatter(x(i), y(i), node_size, ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', 1.5);
        text(x(i)*1.2, y(i)*1.2, roi_all{i}, ...
            'Color', [1 1 1], ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
    end

    title(sprintf('%s - %s: bilateral PPI network (numerosity vs baseline, FDR)', ...
        cfg.subj, cond_label), ...
        'Color', [1 1 1], 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('\n[Step 3] 脑网络可视化完成（bilateral + FDR 校正）。\n');

%% 本地函数：解析 ROI 名的半球信息（首字母 L/R）
function [base_name, hemi] = parse_roi_name(roi_name)
    % roi_name: 字符串，例如 'LNPO', 'RNPO', 'NPO' 等
    % 返回：
    %   base_name: 去掉 L/R 前缀后的名字（若无 L/R，则原样）
    %   hemi     : 'L' / 'R' / ''（无半球信息）
    if isempty(roi_name) || ~ischar(roi_name)
        base_name = roi_name;
        hemi      = '';
        return;
    end
    c1 = upper(roi_name(1));
    if c1 == 'L' || c1 == 'R'
        hemi = c1;
        if numel(roi_name) > 1
            base_name = roi_name(2:end);
        else
            base_name = '';
        end
    else
        hemi = '';
        base_name = roi_name;
    end
end

%% 本地 FDR 函数（Benjamini-Hochberg）
function [h, p_crit] = fdr_bh_local(p_vals, q)
    % p_vals: 列向量或行向量的 p 值（可含 NaN，不会影响结果）
    % q     : FDR 控制水平，例如 0.05
    %
    % 返回：
    %   h      : 与 p_vals 同维度的逻辑向量，true 表示通过 FDR
    %   p_crit : 阈值 p*，即最大的 p，使得 p(k) <= (k/m)*q

    % 去掉 NaN（内部用），但保留外部长度
    p = p_vals(:);
    m = numel(p);
    h = false(size(p));

    valid = ~isnan(p);
    pv    = p(valid);
    if isempty(pv)
        p_crit = NaN;
        h(~valid) = false;
        return;
    end

    [pv_sorted, sort_idx] = sort(pv); %#ok<ASGLU>
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
        % 所有 h 都是 false
    end

    % reshape 回原来形状
    h = reshape(h, size(p_vals));
end
