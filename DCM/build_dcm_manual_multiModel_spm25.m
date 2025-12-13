function build_dcm_manual_multiModel_spm25(subjects)
% 使用 Gray/Original 的 roiSeed 时序，手动构造 *多模型* DCM 结构并估计（SPM25）
% - 节点：6 个双侧数量 map (NF, NPC1-3, NPO, NTO) -> bilateral average
% - 输入：1 个 numerosity block（baseline 为隐式）
% - A：全连接（除对角线）
% - C：NPO 和 NTO 接受 driving input
% - B：7 个候选模型（global / hub / FF / FB / parietal-only / hierarchy / null）
%
% 需要：
%   - SPM12/25 在 path 中（spm_dctmtx, spm_dcm_fmri_priors, spm_dcm_estimate）
%   - 你的 roiSeed .mat（之前 PPI 用的那个）
%   - events_bids.tsv（包含 numerosity / baseline 的 block 信息）

clc;

%% ------------ 0. 全局配置：路径 / 被试列表 / 条件等 ------------

root_prf   = '/home/xue/data/prf_result';   % prf 结果 & roi ts 所在根目录
prf_ses    = 'ses-01';
TR         = 1.95;
NTR_target = 176;

% 条件列表：分别为 connect / unconnect 各建一批 DCM（每个 run × 每个 model 一个）
condList = {'connect','unconnect'};

% *** 在这里设定要跑的被试列表 ***
if nargin < 1 || isempty(subjects)
    subjects = { 'sub-01','sub-02','sub-08'};
end

% events.tsv（block 设计）
events_tsv = '/media/xue/new_B1/glm/events_output/events_bids.tsv';
name_num   = 'numerosity';   % trial_type 中的 numerosity block 名
name_base  = 'baseline';     % baseline 的 trial_type 名（这里仅做 sanity check）

% 模型数 + 名字
nModels = 7;
modelNames = { ...
    'B1_globalGain', ...
    'B2_NPOhub', ...
    'B3_feedforward', ...
    'B4_feedback', ...
    'B5_parietalOnly', ...
    'B6_hierarchy', ...
    'B0_null'};

fprintf('\n=========== 手动构造多模型 DCM (SPM25): ses=%s ===========\n', prf_ses);
fprintf(' root_prf = %s\n', root_prf);
fprintf(' events   = %s\n', events_tsv);

%% ------------ 遍历被试 ------------

for si = 1:numel(subjects)
    subj = subjects{si};
    fprintf('\n########################################################\n');
    fprintf(' >>> SUBJECT: %s (%d/%d)\n', subj, si, numel(subjects));
    fprintf('########################################################\n');

    % roiSeed 时序文件（之前 PPI 用的那个）
    roi_ts_dir = 'results_ppi_roi_tseries_fromOriginal';
    roi_ts_file = sprintf('%s_%s_Original_roiTS_multiRun_dualPRF.mat', subj, prf_ses);
    roi_ts_mat = fullfile(root_prf, subj, prf_ses, roi_ts_dir, roi_ts_file);

    % 输出目录：保存 DCM 结构和拟合结果
    out_dir = fullfile(root_prf, subj, prf_ses, 'DCM_manual_spm25_multiModel');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    fprintf(' ROI ts file: %s\n', roi_ts_mat);
    fprintf(' DCM out dir: %s\n', out_dir);

    %% ------------ 1. 载入 roiSeed 时序，整理成 6 个双侧 ROI ------------

    assert(exist(roi_ts_mat,'file')==2, '找不到 roiSeed 文件: %s', roi_ts_mat);
    S = load(roi_ts_mat);
    roiSeed = S.roiSeed;
    meta = S.meta;

    roiNames = meta.roiNames(:);
    nROI = numel(roiNames);
    NTR = meta.NTR;

    fprintf(' [Info] NTR (per run) = %d, 目标 NTR = %d\n', NTR, NTR_target);
    assert(NTR == NTR_target, 'NTR(%d) != NTR_target(%d)', NTR, NTR_target);

    roiHemi = cell(nROI,1);
    roiBase = cell(nROI,1);
    for i = 1:nROI
        rn = roiNames{i};
        if isempty(rn)
            roiHemi{i} = '';
            roiBase{i} = '';
            continue;
        end
    
        % 只看第一个字符，统一转成大写再判断：'L'/'l' 都当作左半球，'R'/'r' 都当作右半球
        h = upper(rn(1));
        if h == 'L'
            roiHemi{i} = 'L';
            roiBase{i} = rn(2:end);   % 去掉首字母，后面保留原样（可以是小写 nf / npc1 等）
        elseif h == 'R'
            roiHemi{i} = 'R';
            roiBase{i} = rn(2:end);
        else
            % 没有 L/R 前缀的 ROI，当作“无半球标签”
            roiHemi{i} = '';
            roiBase{i} = rn;
        end
    end


    % 所有 base ROI（去掉没有 hemi 的）
    baseROI_all = unique(roiBase(~cellfun(@isempty, roiBase)));

    % 只保留 L/R 都存在的 base ROI
    hasL = false(size(baseROI_all));
    hasR = false(size(baseROI_all));
    for b = 1:numel(baseROI_all)
        bn = baseROI_all{b};
        hasL(b) = any(strcmp(roiBase,bn) & strcmp(roiHemi,'L'));
        hasR(b) = any(strcmp(roiBase,bn) & strcmp(roiHemi,'R'));
    end
    baseROI_bilat = baseROI_all(hasL & hasR);

    % 希望的顺序：NF, NPC1-3, NPO, NTO
    desired_order = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};
    baseROI_list = intersect(desired_order, baseROI_bilat, 'stable');
    nBase = numel(baseROI_list);

    assert(nBase == 6, '双侧 base ROI 数不是 6，请检查 roiNames / baseROI_list');

    % 找出每个 base ROI 对应的 L / R 索引
    idxL = nan(nBase,1);
    idxR = nan(nBase,1);
    for b = 1:nBase
        bn = baseROI_list{b};
        idxL(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'L'), 1);
        idxR(b) = find(strcmp(roiBase,bn) & strcmp(roiHemi,'R'), 1);
        fprintf(' base=%-5s L=%2d (%s) R=%2d (%s)\n', ...
            bn, idxL(b), roiNames{idxL(b)}, idxR(b), roiNames{idxR(b)});
    end

    % 再取出各 ROI 的 index，便于写 B 矩阵
    idx_NF = find(strcmp(baseROI_list, 'NF'));
    idx_NPC1 = find(strcmp(baseROI_list, 'NPC1'));
    idx_NPC2 = find(strcmp(baseROI_list, 'NPC2'));
    idx_NPC3 = find(strcmp(baseROI_list, 'NPC3'));
    idx_NPO = find(strcmp(baseROI_list, 'NPO'));
    idx_NTO = find(strcmp(baseROI_list, 'NTO'));

    %% ------------ 2. 构造 numerosity 的 0/1 输入（不卷积 HRF） ------------

    U_num = make_U_numerosity(events_tsv, NTR_target, TR, name_num, name_base);
    fprintf(' [Info] 已生成 numerosity 输入（DCM.U），长度=%d\n', size(U_num.u,1));

    %% ------------ 3. 高通滤波基函数 X0（DCT） ------------

    hp_cutoff = 256; % 秒
    n_dct = fix(2 * NTR_target * TR / hp_cutoff) + 1;
    n_dct = max(n_dct, 1);

    if exist('spm_dctmtx','file') == 2
        X0 = spm_dctmtx(NTR_target, n_dct); % [NTR × n_dct]
    else
        warning('未找到 spm_dctmtx，高通滤波基函数 X0 将为空');
        X0 = [];
    end

    %% ------------ 4. A / C / D + 7 个 B 模型（不含 run 和条件） ------------

    n = nBase;

    % A: baseline 有效连接（全部 1，除了对角线 0）
    A = ones(n) - eye(n);

    % C: driving input -> NPO, NTO
    C = zeros(n, 1);
    C(idx_NPO, 1) = 1;
    C(idx_NTO, 1) = 1;

    % D: 无非线性调制
    D = zeros(n, n, 0);

    % B_all(:,:,m): 第 m 个模型的 B 矩阵
    B_all = zeros(n, n, nModels);
    allIdx = 1:n;
    parietalIdx = [idx_NF, idx_NPC1, idx_NPC2, idx_NPC3, idx_NPO];
    parietalNoNPO = [idx_NF, idx_NPC1, idx_NPC2, idx_NPC3];

    % ---- B1: global gain（全网统一增益） ----
    B1 = ones(n) - eye(n);
    B_all(:,:,1) = B1;

    % ---- B2: NPO hub（NPO <-> 其他 ROI） ----
    B2 = zeros(n);
    others = setdiff(allIdx, idx_NPO);
    for t = others
        % NPO -> t: 连接 j->i, j=idx_NPO, i=t
        B2(t, idx_NPO) = 1;
        % t -> NPO: j=t, i=idx_NPO
        B2(idx_NPO, t) = 1;
    end
    B_all(:,:,2) = B2;

    % ---- B3: feedforward（NTO -> 所有 + NPO -> 上层 parietal） ----
    B3 = zeros(n);
    % NTO -> (NF, NPC1-3, NPO)
    targets_from_NTO = [parietalNoNPO, idx_NPO];
    for t = targets_from_NTO
        if t ~= idx_NTO
            B3(t, idx_NTO) = 1;
        end
    end
    % NPO -> (NF, NPC1-3)
    for t = parietalNoNPO
        B3(t, idx_NPO) = 1;
    end
    B_all(:,:,3) = B3;

    % ---- B4: feedback（parietal -> NTO） ----
    B4 = zeros(n);
    for s = parietalIdx
        if s ~= idx_NTO
            % s -> NTO: j=s, i=idx_NTO
            B4(idx_NTO, s) = 1;
        end
    end
    B_all(:,:,4) = B4;

    % ---- B5: parietal-only（parietal 内部全连接，NTO 不参与） ----
    B5 = zeros(n);
    for ii = 1:numel(parietalIdx)
        for jj = 1:numel(parietalIdx)
            if parietalIdx(ii) ~= parietalIdx(jj)
                B5(parietalIdx(ii), parietalIdx(jj)) = 1;
            end
        end
    end
    B_all(:,:,5) = B5;

    % ---- B6: hierarchy chain（NTO -> NPO -> NPC3 -> NPC2 -> NPC1 -> NF） ----
    B6 = zeros(n);
    chain = [idx_NTO, idx_NPO, idx_NPC3, idx_NPC2, idx_NPC1, idx_NF];
    for k = 1:(numel(chain)-1)
        src = chain(k);
        trg = chain(k+1);
        % src -> trg: j=src, i=trg
        B6(trg, src) = 1;
    end
    B_all(:,:,6) = B6;

    % ---- B7: null model（全 0） ----
    B7 = zeros(n);
    B_all(:,:,7) = B7;

    fprintf(' [Info] 已构造 7 个 B 模型.\n');

    %% ------------ 5. 对每个条件 × 每个 run × 每个 model 构造并估计 DCM ------------

    for ci = 1:numel(condList)
        condLabel = condList{ci};
        fprintf('\n===== 条件: %s =====\n', condLabel);

        % 检查该条件是否存在
        if ~isfield(roiSeed, condLabel)
            warning(' roiSeed 中没有条件 %s，跳过', condLabel);
            continue;
        end

        % 确定该条件下有多少个 run（用其中一个 ROI 的列数判断）
        example_roi_name = roiNames{idxL(1)}; % 比如 LNF
        this_ts = roiSeed.(condLabel).(example_roi_name);
        nRun_cond = size(this_ts, 2);
        fprintf(' [Info] 条件 %s 下检测到 nRun_cond = %d\n', condLabel, nRun_cond);

        for rr = 1:nRun_cond
            fprintf('\n----- 构造 DCM: subj=%s, cond=%s, run=%d/%d -----\n', ...
                subj, condLabel, rr, nRun_cond);

            % --------- 5.1 构造 6 个 bilateral ROI 的 BOLD 时序矩阵 Y ---------
            % Y: [NTR × 6]，每列一个 base ROI（L/R 平均）
            Y = nan(NTR_target, nBase);
            for b = 1:nBase
                iL = idxL(b);
                iR = idxR(b);
                rnL = roiNames{iL};
                rnR = roiNames{iR};

                yL = roiSeed.(condLabel).(rnL)(:, rr); % [NTR × 1]
                yR = roiSeed.(condLabel).(rnR)(:, rr); % [NTR × 1]
                Y(:,b) = mean([yL, yR], 2, 'omitnan');
            end

            % --------- 5.2 构造不含 B 的 DCM 基础结构（模板） ---------
            DCM_base = struct();
            DCM_base.n = nBase; % 节点数
            DCM_base.v = NTR_target; % 时间点数

            % 数据 Y
            DCM_base.Y.y = Y; % [v × n]
            DCM_base.Y.dt = TR;
            DCM_base.Y.name = baseROI_list; % 节点名称
            DCM_base.Y.X0 = X0; % 高通滤波基函数

            % 输入 U: numerosity block, baseline 隐式
            DCM_base.U = U_num; % make_U_numerosity 返回的结构
            DCM_base.U.dt = TR;

            % 选项（贴近 SPM DCM GUI 的 fMRI 设置）
            DCM_base.delays = repmat(0.5, nBase, 1); % hemodynamic delay (0.5 * TR)
            DCM_base.TE = 0.03; % 30 ms

            DCM_base.options = struct();
            DCM_base.options.nonlinear = 0;
            DCM_base.options.two_state = 0;
            DCM_base.options.stochastic = 0;
            DCM_base.options.centre = 1;
            DCM_base.options.induced = 0;
            DCM_base.options.maxnodes = nBase;
            DCM_base.options.endogenous = 0;
            DCM_base.options.analysis = 'fMRI';

            % A / C / D 固定
            DCM_base.a = A;
            DCM_base.c = C;
            DCM_base.d = D;

            % --------- 5.3 针对 7 个 B 模型分别构建 & 估计 ---------
            for mm = 1:nModels
                DCM = DCM_base;

                DCM.b = B_all(:,:,mm);

                % 一些 meta 信息，方便后续 QC / 解析
                DCM.subj = subj;
                DCM.cond = condLabel;
                DCM.run = rr;
                DCM.model_id = mm;
                DCM.model_name = modelNames{mm};

                % 文件名：加上 modelXX
                dcm_name = fullfile(out_dir, ...
                    sprintf('%s_%s_DCM_%s_run%02d_model%02d_spm25.mat', ...
                    subj, prf_ses, condLabel, rr, mm));

                % 先保存未估计版本（以防估计报错）
                save(dcm_name, 'DCM', '-v7.3');
                fprintf(' [Save] 未估计 DCM: %s (model %d: %s)\n', ...
                    dcm_name, mm, modelNames{mm});

                % ---- 估计 DCM ----
                if exist('spm_dcm_fmri_priors','file') == 2 && exist('spm_dcm_estimate','file') == 2
                    try
                        % fMRI DCM 先验
                        [pE, pC] = spm_dcm_fmri_priors(DCM.a, DCM.b, DCM.c, DCM.d, DCM.options);

                        DCM.M = struct();
                        DCM.M.pE = pE;
                        DCM.M.pC = pC;

                        fprintf(' [Est] spm_dcm_estimate... (model %d: %s)\n', ...
                            mm, modelNames{mm});
                        DCM = spm_dcm_estimate(DCM);

                        save(dcm_name, 'DCM', '-v7.3');
                        fprintf(' [Est] 完成并写回: %s\n', dcm_name);

                    catch ME
                        warning(' [Est] 估计 DCM model %d 时出错: %s\n 保留未估计版本。', ...
                            mm, ME.message);
                    end
                else
                    warning(' [Est] 未找到 spm_dcm_fmri_priors 或 spm_dcm_estimate，仅保存 DCM 结构，不估计。');
                end

            end % models

        end % runs
    end % condList

    fprintf('\n=========== SUBJECT %s: 多模型 DCM 构造 & 估计完成 ===========\n', subj);
end % subjects

fprintf('\n=========== 所有被试的多模型 DCM 构造 & 估计完成 ===========\n');

end % function build_dcm_manual_multiModel_spm25


%% ================== 辅助函数 ==================
function U = make_U_numerosity(events_tsv, NTR, TR, name_num, name_base)
% 为 DCM 构造一个单一输入：numerosity block 的 0/1 boxcar（不卷积 HRF）
% baseline 作为隐式对照，这里只用 numerosity 的 on/off 信息
%
% U.u : [NTR × 1]，每个 TR 一个采样
% U.name = {'numerosity'}

assert(exist(events_tsv,'file')==2, 'events.tsv 不存在: %s', events_tsv);

T = readtable(events_tsv, 'FileType','text', 'Delimiter','\t');
req_cols = {'onset','duration','trial_type'};
missing = setdiff(req_cols, T.Properties.VariableNames);
assert(isempty(missing), 'events.tsv 缺少列: %s', strjoin(missing, ', '));

t = (0:NTR-1)' * TR; % 每个 TR 的时间点

% numerosity block
mask_num = strcmp(T.trial_type, name_num);
assert(any(mask_num), 'trial_type 中找不到 "%s"', name_num);

onset_num = T.onset(mask_num);
dur_num   = T.duration(mask_num);

u_num = zeros(NTR,1);
for k = 1:numel(onset_num)
    t_on  = onset_num(k);
    t_off = onset_num(k) + dur_num(k);
    u_num = u_num + (t >= t_on & t < t_off);
end

% baseline 只是做个 sanity check，不构造显式输入
if any(strcmp(T.trial_type, name_base))
    fprintf(' [make_U] 检测到 baseline 事件 (%s)，但不建显式输入（隐式 baseline）。\n', name_base);
else
    fprintf(' [make_U] 未检测到 baseline 行，将其视为隐式 baseline。\n');
end

U = struct();
U.u = u_num;
U.name = {'numerosity'};
end

