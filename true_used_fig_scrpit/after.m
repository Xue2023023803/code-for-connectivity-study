%% plot_BMRLQ_3panel_CURVED_SWALLOWTAIL_EPS.m
% ============================================================
% 3子图（1×3）横向：BMR L/Q/LB 的“弯曲燕尾箭头”网络（EPS 矢量）
%
% A: Linear (Pp_incl_Linear >= PP_THR) 的有向网络
% B: Quadratic (Pp_incl_Quadratic >= PP_THR) 的有向网络
% C: 联合下界 LB = max(0, P(L)+P(Q)-1) >= PP_THR 的有向网络
%
% 可视化参考 BrainNet 脚本思路：
%   - 曲线：三次 Bézier，曲率随跨度变化，双向边额外分离
%   - 箭头：燕尾箭头 patch（带燕尾缺口）
%   - 简单避碰：若与已画曲线距离过近，逐步增大外弯曲率
%   - EPS：painters，且同步 PaperSize 防裁剪
%
% 方向约定（与 BMR mat 一致）：
%   (i,j) 表示 j -> i（source=j, target=i）
% ============================================================

clear; clc; close all;

%% ===================== 用户设置（只改这里） =====================
root_prf      = '/home/xue/data/prf_result';
prf_ses       = 'ses-01';

MATRIX_MODE   = 'AB';    % 'A' / 'B' / 'AB'
FIELD_TO_PLOT = 'A';     % 'A' 或 'B'
PP_THR        = 0.95;    % inclusion 概率阈值（标题不显示）

% ROI 顺序（必须与矩阵索引一致）
roi_labels = {'NF','NPC1','NPC2','NPC3','NPO','NTO'};

% 用于“跨度->曲率分层”的顺时针顺序（仅用于计算曲率；不改变节点摆放顺序）
roi_order_curve = {'NTO','NPO','NPC1','NPC2','NPC3','NF'};

% 颜色（单色网络）
col_L  = [0.00 0.45 0.74];   % Linear
col_Q  = [0.85 0.33 0.10];   % Quadratic
col_LB = [0.49 0.18 0.56];   % L&Q(LB)

% 节点与布局
NODE_R = 0.085;              % 节点半径（用于缩短边，避免箭头插入节点）
NODE_SIZE_SCATTER = 140;

% “弯曲燕尾箭头”几何参数（均为比例参数，更稳定）
% —— 宽度：基于 chord 长度 L，乘以一个由 Pp 映射得到的缩放因子（见下方）
BASE_W_FRAC  = 0.030;        % 箭杆基础宽度 = BASE_W_FRAC * L
HEAD_W_FRAC  = 0.060;        % 箭头基底宽度 = HEAD_W_FRAC * L
TAIL_W_FRAC  = 0.050;        % 燕尾外宽 = TAIL_W_FRAC * L
NOTCH_W_FRAC = 0.020;        % 燕尾缺口宽 = NOTCH_W_FRAC * L

HEAD_LEN_FRAC  = 0.18;       % 箭头长度 = HEAD_LEN_FRAC * L
TAIL_LEN_FRAC  = 0.12;       % 燕尾长度 = TAIL_LEN_FRAC * L
NOTCH_LEN_FRAC = 0.08;       % 缺口深度 = NOTCH_LEN_FRAC * L

% 线宽/粗细映射：用你原先 EDGE_W_MIN/MAX 的思想来控制“宽度缩放”
EDGE_W_MIN = 1.5;            % 对应 Pp=0 时的相对最小（用于缩放）
EDGE_W_MAX = 4.5;            % 对应 Pp=1 时的相对最大（用于缩放）

% 曲率：随“跨度”变化（更外弯更不重叠）
CURV_BASE_FRAC = 0.03;       % 基础弯曲幅度（占弦长比例）
CURV_STEP_FRAC = 0.06;       % 随跨度增加的弯曲增量（占弦长比例）
DIR_SEP_FRAC   = 0.02;       % 同一对 ROI 的反向边额外分离（占弦长比例）

% 避碰（2D 距离单位为当前坐标系，圆半径≈1）
AVOID_ENABLE      = true;
AVOID_MIN_DIST    = 0.10;    % 新曲线与已画曲线的最小距离阈值
AVOID_MAX_ITERS   = 10;      % 每条边最多迭代次数
AVOID_STEP_FRAC   = 0.04;    % 每次增加的曲率步长（占弦长比例）
AVOID_MAX_CURV_FR = 0.35;    % 曲率最大上限（占弦长比例）

% 画布与标题/图例布局
FIG_W_CM  = 46;
FIG_H_CM  = 16;
TITLE_DY  = 0.08;            % 标题上移（normalized）
TL_OUTERPOS = [0.02 0.16 0.96 0.72];  % [left bottom width height]，留足上方空白

%% ===================== 自动定位 BMR mat =====================
MATRIX_MODE   = upper(string(MATRIX_MODE));
FIELD_TO_PLOT = upper(char(FIELD_TO_PLOT));
assert(ismember(MATRIX_MODE, ["A","B","AB"]), 'MATRIX_MODE 必须是 A/B/AB');
assert(ismember(FIELD_TO_PLOT, ['A','B']), 'FIELD_TO_PLOT 必须是 A 或 B');

bmr_dir = fullfile(root_prf, sprintf('PEB_group_AL3cond_poly%s', char(MATRIX_MODE)), ...
    prf_ses, 'BMR_designCompare_LQ');

bmr_mat = fullfile(bmr_dir, sprintf('BMR_LQ_designCompare_%s.mat', char(MATRIX_MODE)));
assert(exist(bmr_mat,'file')==2, '找不到 BMR mat：%s', bmr_mat);
S = load(bmr_mat);

assert(isfield(S,'edgeInclusion') && isfield(S.edgeInclusion, FIELD_TO_PLOT), ...
    'BMR mat 中缺少 edgeInclusion.%s', FIELD_TO_PLOT);

EI  = S.edgeInclusion.(FIELD_TO_PLOT);
PpL = EI.Pp_incl_Linear;
PpQ = EI.Pp_incl_Quadratic;

nROI = min([size(PpL,1), size(PpL,2), numel(roi_labels)]);
PpL = PpL(1:nROI,1:nROI);
PpQ = PpQ(1:nROI,1:nROI);
roi_labels = roi_labels(1:nROI);

% (i,j) 表示 j->i
needL = (PpL >= PP_THR); needL(eye(nROI)==1) = false;
needQ = (PpQ >= PP_THR); needQ(eye(nROI)==1) = false;

% 联合下界（Fréchet 下界）
PpLB = max(0, PpL + PpQ - 1);
needLB = (PpLB >= PP_THR); needLB(eye(nROI)==1) = false;

%% ===================== 绘图：1×3 横向（EPS 矢量） =====================
fig = figure('Color','w','Renderer','painters');
set(fig,'Units','centimeters','Position',[2 2 FIG_W_CM FIG_H_CM]);
set(fig,'InvertHardcopy','off');
sync_paper_to_figure(fig);

TL = tiledlayout(1,3,'TileSpacing','compact','Padding','loose');
TL.OuterPosition = TL_OUTERPOS;

% （不改参数：这里用一个固定额外偏移，避免标题压到 ROI 文字）
TITLE_DY_EXTRA = 0.03;

% A：Linear
ax1 = nexttile(TL,1);
plot_network_curved_swallowtail(ax1, nROI, roi_labels, roi_order_curve, needL, PpL, ...
    col_L, NODE_R, NODE_SIZE_SCATTER, ...
    BASE_W_FRAC, HEAD_W_FRAC, TAIL_W_FRAC, NOTCH_W_FRAC, ...
    HEAD_LEN_FRAC, TAIL_LEN_FRAC, NOTCH_LEN_FRAC, ...
    EDGE_W_MIN, EDGE_W_MAX, ...
    CURV_BASE_FRAC, CURV_STEP_FRAC, DIR_SEP_FRAC, ...
    AVOID_ENABLE, AVOID_MIN_DIST, AVOID_MAX_ITERS, AVOID_STEP_FRAC, AVOID_MAX_CURV_FR);

t1 = title(ax1, sprintf('Linear | Field %s', FIELD_TO_PLOT), 'FontWeight','bold');
move_title_up_normalized(t1, TITLE_DY + TITLE_DY_EXTRA);

% B：Quadratic
ax2 = nexttile(TL,2);
plot_network_curved_swallowtail(ax2, nROI, roi_labels, roi_order_curve, needQ, PpQ, ...
    col_Q, NODE_R, NODE_SIZE_SCATTER, ...
    BASE_W_FRAC, HEAD_W_FRAC, TAIL_W_FRAC, NOTCH_W_FRAC, ...
    HEAD_LEN_FRAC, TAIL_LEN_FRAC, NOTCH_LEN_FRAC, ...
    EDGE_W_MIN, EDGE_W_MAX, ...
    CURV_BASE_FRAC, CURV_STEP_FRAC, DIR_SEP_FRAC, ...
    AVOID_ENABLE, AVOID_MIN_DIST, AVOID_MAX_ITERS, AVOID_STEP_FRAC, AVOID_MAX_CURV_FR);

t2 = title(ax2, sprintf('Quadratic | Field %s', FIELD_TO_PLOT), 'FontWeight','bold');
move_title_up_normalized(t2, TITLE_DY + TITLE_DY_EXTRA);

% C：LB
ax3 = nexttile(TL,3);
plot_network_curved_swallowtail(ax3, nROI, roi_labels, roi_order_curve, needLB, PpLB, ...
    col_LB, NODE_R, NODE_SIZE_SCATTER, ...
    BASE_W_FRAC, HEAD_W_FRAC, TAIL_W_FRAC, NOTCH_W_FRAC, ...
    HEAD_LEN_FRAC, TAIL_LEN_FRAC, NOTCH_LEN_FRAC, ...
    EDGE_W_MIN, EDGE_W_MAX, ...
    CURV_BASE_FRAC, CURV_STEP_FRAC, DIR_SEP_FRAC, ...
    AVOID_ENABLE, AVOID_MIN_DIST, AVOID_MAX_ITERS, AVOID_STEP_FRAC, AVOID_MAX_CURV_FR);

t3 = title(ax3, sprintf('L&Q (LB) | Field %s', FIELD_TO_PLOT), 'FontWeight','bold');
move_title_up_normalized(t3, TITLE_DY + TITLE_DY_EXTRA);

%% ===================== 单一共享图例（挂在 ax2 上） =====================
hold(ax2,'on');
hL  = plot(ax2, nan, nan, '-', 'Color', col_L,  'LineWidth', 4);
hQ  = plot(ax2, nan, nan, '-', 'Color', col_Q,  'LineWidth', 4);
hLB = plot(ax2, nan, nan, '-', 'Color', col_LB, 'LineWidth', 4);
hold(ax2,'off');

lgd = legend(ax2, [hL hQ hLB], {'Linear','Quadratic','L&Q'}, ...
    'Location','southoutside', 'Orientation','horizontal');
set(lgd,'Box','off');

% （不改参数：图例再略微下移，避免和 ROI 文字重叠）
move_legend_down_normalized(lgd, 0.02);

% 加 legend 后再同步一次，避免裁剪
sync_paper_to_figure(fig);

%% ===================== 导出 EPS（矢量） =====================
out_eps = fullfile(root_prf, sprintf('BMRLQ_3panel_CURVED_SWALLOWTAIL_poly%s_Field%s_PP%03d.eps', ...
    char(MATRIX_MODE), FIELD_TO_PLOT, round(PP_THR*100)));

print(fig, out_eps, '-depsc2', '-painters');
fprintf('[SAVED] %s\n', out_eps);

%% =====================================================================
%% ========================= 本地函数区 ================================
%% =====================================================================

function plot_network_curved_swallowtail(ax, nROI, roi_labels, roi_order_curve, maskToPlot, PpMat, ...
    colEdge, NODE_R, NODE_SIZE_SCATTER, ...
    BASE_W_FRAC, HEAD_W_FRAC, TAIL_W_FRAC, NOTCH_W_FRAC, ...
    HEAD_LEN_FRAC, TAIL_LEN_FRAC, NOTCH_LEN_FRAC, ...
    EDGE_W_MIN, EDGE_W_MAX, ...
    CURV_BASE_FRAC, CURV_STEP_FRAC, DIR_SEP_FRAC, ...
    AVOID_ENABLE, AVOID_MIN_DIST, AVOID_MAX_ITERS, AVOID_STEP_FRAC, AVOID_MAX_CURV_FR)

axes(ax); cla(ax);
axis(ax,'equal'); axis(ax,'off'); hold(ax,'on');

% ---------- 节点圆形布局（从上方开始顺时针） ----------
theta = pi/2 - (0:nROI-1)*(2*pi/nROI);
x = cos(theta);
y = sin(theta);
center2d = [0 0];

% ---------- 关键修正 1：全局常数长度（与“边本身长度”无关） ----------
% 目的：彻底消灭 “BASE_W_FRAC*L” 造成的宽度差异，三个子图也天然一致
L_SCALE = compute_global_L_scale(x, y, NODE_R);    % 只由节点几何决定的常数

% ---------- 关键修正 2：固定轴范围，避免不同子图自动缩放导致“看起来不一样粗” ----------
AX_LIM = 1.45;   % 固定一点 margin，保证弯曲与文字都不裁
xlim(ax, [-AX_LIM AX_LIM]);
ylim(ax, [-AX_LIM AX_LIM]);

% ---------- ROI 顺序映射（用于跨度->曲率分层） ----------
orderPos = nan(1,nROI);
for ii = 1:nROI
    orderPos(ii) = find(strcmp(roi_order_curve, roi_labels{ii}), 1);
end
hasOrder = all(isfinite(orderPos));

% ---------- 收集要画的边，并按“强度/跨度”排序：先画更强/更长的，避碰更稳 ----------
[ii, jj] = find(maskToPlot);
E = struct('i',{},'j',{},'pp',{},'L',{},'span',{});
for k = 1:numel(ii)
    i = ii(k); j = jj(k);
    p0 = [x(j), y(j)];
    p3 = [x(i), y(i)];
    L = norm(p3 - p0);
    if ~isfinite(L) || L < 1e-6, continue; end

    pp = PpMat(i,j);
    if ~isfinite(pp), pp = 0; end

    if hasOrder
        ps = orderPos(j);
        pt = orderPos(i);
        d1 = mod(pt - ps, nROI);
        d2 = mod(ps - pt, nROI);
        span = min(d1,d2);
    else
        span = 1;
    end

    E(end+1).i = i;
    E(end).j = j;
    E(end).pp = pp;
    E(end).L = L;
    E(end).span = span;
end

if ~isempty(E)
    % 排序：先 pp 高，再 L 长（更重要的边优先占“好轨道”）
    [~, ord] = sortrows([[E.pp].', [E.L].'], [-1 -1]);
    E = E(ord);
end

% ---------- 用于避碰：缓存已放置曲线（中心线的降采样点） ----------
placedCurves2d = {};

% ---------- 画边（弯曲燕尾箭头 patch） ----------
for e = 1:numel(E)
    i  = E(e).i;
    j  = E(e).j;

    % 起点/终点（节点圆心）
    P0 = [x(j), y(j)];
    P3 = [x(i), y(i)];

    chord = P3 - P0;
    L0 = norm(chord);
    if L0 < 1e-6, continue; end
    uChord = chord / L0;

    % 缩短端点（避免进入节点）
    P0s = P0 + NODE_R * uChord;
    P3e = P3 - NODE_R * uChord;

    chord2 = P3e - P0s;
    L = norm(chord2);                    % 这是真实边长（只用于曲率/控制点）
    if L < 1e-6, continue; end
    u = chord2 / L;

    % 外法向（朝外侧）
    n2 = [-u(2), u(1)];
    mid = (P0s + P3e)/2;
    if dot(n2, mid - center2d) < 0
        n2 = -n2;
    end

    % 跨度归一化
    if hasOrder
        dNorm = E(e).span / max(1, floor(nROI/2));
    else
        dNorm = 0.5;
    end

    % 基础曲率（外弯）
    curvMag0 = L * (CURV_BASE_FRAC + CURV_STEP_FRAC * dNorm);

    % 如果存在反向边（i<->j），让两条边弯曲方向相反以错开
    if maskToPlot(j,i)
        dirFlag = sign(j - i);
        if dirFlag == 0, dirFlag = 1; end
        n2 = n2 * dirFlag;
        curvMag0 = curvMag0 + (DIR_SEP_FRAC * L);
    end

    % 避碰：若太近则迭代增加外弯
    curvMag = curvMag0;
    bestC = [];
    bestMinDist = -Inf;

    nTry = 1;
    if AVOID_ENABLE, nTry = max(1, AVOID_MAX_ITERS); end

    for it = 1:nTry
        curvMag = min(curvMag, AVOID_MAX_CURV_FR * L);

        % 三次 Bézier：两个控制点都放在“沿弦前进 + 外法向偏移”
        P1 = P0s + 0.28*chord2 + curvMag*n2;
        P2 = P0s + 0.72*chord2 + curvMag*n2;

        t = linspace(0,1,60)';
        C = bezier_cubic(P0s,P1,P2,P3e,t);

        % 计算与已放置曲线的最小距离（用降采样加速）
        Cq = downsample_curve(C, 25);
        if isempty(placedCurves2d)
            minDist = Inf;
        else
            minDist = min_distance_to_placed(Cq, placedCurves2d);
        end

        if minDist > bestMinDist
            bestMinDist = minDist;
            bestC = C;
        end

        if ~AVOID_ENABLE || minDist >= AVOID_MIN_DIST
            bestC = C;
            break;
        end

        % 太近：加大外弯
        curvMag = curvMag + (AVOID_STEP_FRAC * L);
        if curvMag >= (AVOID_MAX_CURV_FR * L)
            break;
        end
    end

    if isempty(bestC), continue; end
    C = bestC;

    % ---------- 关键：宽度完全固定（不再用每条边自己的 L） ----------
    % 不改参数：仍然沿用 EDGE_W_MIN/MAX 的“中间值”作为统一缩放
    widthScale = (EDGE_W_MIN + EDGE_W_MAX) / (2*EDGE_W_MAX);

    w_base  = (BASE_W_FRAC  * L_SCALE) * widthScale;
    w_head  = (HEAD_W_FRAC  * L_SCALE) * widthScale;
    w_tail  = (TAIL_W_FRAC  * L_SCALE) * widthScale;
    w_notch = (NOTCH_W_FRAC * L_SCALE) * widthScale;

    headLen  = HEAD_LEN_FRAC  * L_SCALE;
    tailLen  = TAIL_LEN_FRAC  * L_SCALE;
    notchLen = NOTCH_LEN_FRAC * L_SCALE;

    % 生成“燕尾箭头”多边形（2D）
    poly2d = swallowtail_arrow_polygon(C, w_base, w_head, headLen, w_tail, tailLen, w_notch, notchLen);

    % 画 patch（矢量）
    patch(ax, poly2d(:,1), poly2d(:,2), colEdge, ...
        'EdgeColor','none', 'FaceAlpha', 1.0, 'Clipping','off');

    % 记录已放置曲线（用于后续避碰）
    placedCurves2d{end+1} = downsample_curve(C, 25);
end

% ---------- 画节点（上层） ----------
scatter(ax, x, y, NODE_SIZE_SCATTER, 'w', 'filled', ...
    'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth', 1.2);

% （不改参数：ROI 标签距离略微增大）
LABEL_R = 1.20;
for n = 1:nROI
    text(ax, LABEL_R*x(n), LABEL_R*y(n), roi_labels{n}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',10);
end

hold(ax,'off');
end

function Ls = compute_global_L_scale(x, y, NODE_R)
% 只依赖“节点圆环几何”的常数长度：用于统一宽度/箭头长度（三子图一致）
nROI = numel(x);
Ls_all = zeros(nROI*(nROI-1),1);
c = 0;
for i = 1:nROI
    for j = 1:nROI
        if i==j, continue; end
        P0 = [x(j), y(j)];
        P3 = [x(i), y(i)];
        v  = P3 - P0;
        Lv = norm(v);
        if Lv < 1e-12, continue; end
        u  = v / Lv;
        P0s = P0 + NODE_R*u;
        P3e = P3 - NODE_R*u;
        c = c + 1;
        Ls_all(c) = norm(P3e - P0s);
    end
end
Ls_all = Ls_all(1:c);
% 用均值（比 min 更不容易“看起来太细”，且仍然是全局常数）
Ls = mean(Ls_all);
end

function poly = swallowtail_arrow_polygon(C, w_base, w_head, headLen, w_tail, tailLen, notchW, notchLen)
% 输入：C 为 Nx2 中心曲线点（从尾到头）
% 输出：poly 为封闭多边形（燕尾箭头 patch）
N = size(C,1);
if N < 12
    poly = C; return;
end

% 切线与法线（2D）
dC = gradient(C);
T  = unitvec_rows(dC);
Nn = [-T(:,2), T(:,1)];

% 弧长参数
ds = sqrt(sum(diff(C,1,1).^2,2));
s  = [0; cumsum(ds)];
L  = s(end); if L<=0, L=1; end
sn = s / L;

% 宽度函数：中间 w_base，尾端->w_tail，头端->w_head（指数过渡）
w = w_base * ones(N,1);
tailSigma = max(0.02, tailLen / L);
headSigma = max(0.02, headLen / L);

w = w + (w_tail - w_base) .* exp(-sn / tailSigma);
w = w + (w_head - w_base) .* exp(-(1-sn) / headSigma);

% 头部起点：距离末端 headLen
s_to_end = L - s;
iHead = find(s_to_end <= headLen, 1, 'first');
if isempty(iHead), iHead = round(0.85*N); end
iHead = max(6, min(N-4, iHead));

% 燕尾缺口位置：距离起点 notchLen
iNotch = find(s >= notchLen, 1, 'first');
if isempty(iNotch), iNotch = 4; end
iNotch = max(4, min(N-8, iNotch));

% 主体边界（到 iHead 为止）
idx = (1:iHead)';
Lft = C(idx,:) + 0.5*w(idx).*Nn(idx,:);
Rgt = C(idx,:) - 0.5*w(idx).*Nn(idx,:);

% 头部（三角）
Ptip = C(end,:);
Pb   = C(iHead,:);
Nb   = Nn(iHead,:);
hbL  = Pb + 0.5*w_head*Nb;
hbR  = Pb - 0.5*w_head*Nb;

% 燕尾尾部：外点 + 缺口点
P0 = C(1,:);
N0 = Nn(1,:);
tailL = P0 + 0.5*w_tail*N0;
tailR = P0 - 0.5*w_tail*N0;

Pn = C(iNotch,:);
Nn0 = Nn(iNotch,:);
notchL = Pn + 0.5*notchW*Nn0;
notchR = Pn - 0.5*notchW*Nn0;

% 多边形组装（顺时针/逆时针都可，patch 会自动填充）
poly = [ ...
    Lft; ...
    hbL; ...
    Ptip; ...
    hbR; ...
    flipud(Rgt); ...
    tailR; ...
    notchR; ...
    notchL; ...
    tailL ...
    ];
end

function Cq = downsample_curve(C, k)
% 均匀降采样到 k 个点（含首尾）
N = size(C,1);
if N <= k, Cq = C; return; end
idx = round(linspace(1, N, k));
idx = unique(max(1, min(N, idx)));
Cq = C(idx,:);
end

function dmin = min_distance_to_placed(Cq, placedCurves2d)
% 点集 Cq 到所有已放置曲线点集的最小距离
dmin2 = Inf;
for m = 1:numel(placedCurves2d)
    P = placedCurves2d{m};
    dx = Cq(:,1) - P(:,1)';
    dy = Cq(:,2) - P(:,2)';
    d2 = dx.^2 + dy.^2;
    dmin2 = min(dmin2, min(d2(:)));
end
dmin = sqrt(dmin2);
end

function p = bezier_cubic(p0,p1,p2,p3,t)
% 三次 Bézier 曲线
t = t(:);
u = 1 - t;
p = (u.^3)*p0 + 3*(u.^2).*t*p1 + 3*u.*(t.^2)*p2 + (t.^3)*p3;
end

function V = unitvec_rows(V)
n = sqrt(sum(V.^2,2));
n(n<1e-12) = 1;
V = V ./ n;
end

function move_title_up_normalized(th, dy)
if isempty(th) || ~isgraphics(th), return; end
set(th,'Units','normalized');
pos = get(th,'Position');
pos(2) = pos(2) + dy;
set(th,'Position',pos);
end

function move_legend_down_normalized(lgd, dy)
if isempty(lgd) || ~isgraphics(lgd), return; end
set(lgd,'Units','normalized');
pos = get(lgd,'Position');
pos(2) = pos(2) - dy;
set(lgd,'Position',pos);
end

function sync_paper_to_figure(fig)
% 同步 PaperSize/PaperPosition，避免 EPS 被裁剪
set(fig,'Units','centimeters');
pos = get(fig,'Position');              % [x y w h]
set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize', pos(3:4));
set(fig,'PaperPosition', [0 0 pos(3:4)]);
set(fig,'PaperPositionMode','manual');
end


