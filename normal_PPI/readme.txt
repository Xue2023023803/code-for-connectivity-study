ROI–ROI PPI 四步脚本说明（readme）

本 readme 说明以下四个脚本如何配合完成从单被试到组水平的 ROI–ROI PPI 分析与可视化：

ppi_roi_tseries_fromOriginal_cmd.m – 从 Gray/Original 提取多 run ROI 时序（双 pRF 条件，严格版）

ppi_roi_roi_runwise_cmd_tough.m – 单被试 run-wise ROI–ROI PPI + 左右半球固定效应合并（严格版）

ppi_graph_analysis_group.m – 组水平图论指标分析（strength/degree/clustering）

ppi_group_network_MNI_sagittal.m – 在 SPM avg152T1 MNI 模板上可视化组水平网络（connect / unconnect / 差值）

一、依赖环境与数据准备

软件与工具：

mrVista（VOLUME/dataTYPES/Gray/Original/pRF dataType）
MATLAB（脚本运行环境）
SPM25（spm_hrf, spm_dctmtx, canonical/avg152T1.nii 等）
fMRIPrep 预处理结果（*_desc-confounds_timeseries.tsv、*_space-T1w_desc-preproc_bold.nii[.gz] 等）
必要数据：
mrVista 会话目录 subj/ses/Gray：
coords.mat
ROIs/*.mat（仅使用以 L/R 开头并包含数量 map 的 ROI，如 LNPC1, RNPC1 等）
pRF dataType：
Averages_unconnect_layer_01（unconnect pRF）
Averages_connect_layer_01（connect pRF）
原始 BOLD dataType：Original
BIDS_prep 根目录（fMRIPrep 输出）
统一的事件文件：events_bids.tsv（包含 numerosity / baseline 两类 trial_type，用于构造 Psych_diff = numerosity − baseline）

二、整体分析流程概览

分析分为四步，从 ROI 时序提取 → 单被试 PPI → 组水平图论 → 组水平 MNI 可视化：
Step 1（单被试 / 单会话）
ppi_roi_tseries_fromOriginal_cmd.m
从 Gray/Original 提取多 run 的 ROI 时序（connect/unconnect 两套 pRF 并行）
生成用于 PPI 的 ROI seed（平均时序）与体素级时序

Step 2（单被试）
ppi_roi_roi_runwise_cmd_tough.m
结合 fMRIPrep confounds 和 events_bids.tsv，对每个 run 构建 ROI–ROI PPI GLM（numerosity vs baseline 的 Psych_diff）
跨 run 固定效应 + 左右半球固定效应合并，得到被试级 bilateral ROI–ROI PPI 矩阵（connect / unconnect）

Step 3（组水平）
ppi_graph_analysis_group.m

汇总所有被试的 bilateral PPI 矩阵，计算每个条件下每个 ROI 的：

节点强度 node strength

度数 degree（基于 |β| 归一化后的阈值）

聚类系数 clustering coefficient

输出 group-level mean ± SEM，并绘制组平均网络图（圆形布局）

Step 4（组水平 MNI 可视化）
ppi_group_network_MNI_sagittal.m

使用 Harvey & Dumoulin (2017) 的 MNI 坐标（左右平均）放置 6 个数量 map 节点（NTO/NPO/NPC1–3/NF）

在 SPM avg152T1.nii 的矢状切片上绘制：

connect 条件网络

unconnect 条件网络

connect − unconnect 差异网络

对每条边做 group-level t 检验 + FDR (BH)，按显著性与符号着色
