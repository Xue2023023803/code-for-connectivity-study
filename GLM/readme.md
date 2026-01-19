# fMRI GLM Analysis & ROI Generation Pipeline

这是一个基于 Python (`nilearn`, `nibabel`) 和 FreeSurfer 的 fMRI 数据分析流水线。文档重点标注了必须修改的配置（如实验设计逻辑、路径）以及每个脚本背后的计算原理，方便他人接手或复用代码。

---

## 📖 主要功能
* **一级分析 (First Level)**：针对每个 Run 进行 GLM 建模。
* **固定效应分析 (Fixed Effects)**：合并 Session 内的多个 Run，生成高信噪比的激活图。
* **ROI 定义**：基于灰质掩膜、FDR 校正和簇大小限制，自动生成功能性 ROI。
* **可视化**：提供 3D 体素和皮层表面 (Surface) 的交互式网页可视化。

---

## ⚙️ 环境依赖
* **Python 3.8+**
* **Python 库**: `numpy`, `pandas`, `nibabel`, `nilearn`, `scipy`
* **外部软件**: [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) (仅 `vis_surface.py` 需要，用于表面投影)

---

## 📂 数据目录结构要求
在运行代码前，请确保目录结构如下（需在脚本中修改 `ROOT` 路径）：

```plaintext
/your/project/root/ (ROOT)
├── BIDS_prep/                  # 预处理数据 (fMRIPrep 输出或 BIDS 格式)
│   ├── sub-01/
│   │   ├── ses-01/func/        # 必须包含 space-T1w_desc-preproc_bold.nii.gz
│   │   └── anat/               # 必须包含 label-GM_probseg.nii.gz (灰质概率图)
├── events_output/
│   └── events_bids.tsv         # 事件文件 (必须包含 onset, duration, trial_type)
├── first_step/                 # [自动生成] 一级分析输出
└── roi_step/                   # [自动生成] 固定效应与ROI输出
```

---

## 🚀 脚本使用指南

### 1. run_firstlevel.py (一级分析)
**功能**：对每个 Run 单独进行 GLM 建模，计算 Beta 和 Z 统计图。

**🧮 计算细节**：
* **预处理**：内存中裁切前 `DROP_FIRST` 个 TR，保留随后的 `KEEP_N` 个 TR。
* **设计矩阵**：使用 glover HRF 模型，包含 cosine 漂移项 (High-pass 1/128Hz)。自动加载 Confounds (头动等)。
* **噪声模型**：AR(1)。
* **对比度**：计算 `numerosity` vs 隐含基线。

**⚠️ 必须修改的配置 (User Config)**：
* `ROOT`: 修改为你的项目根目录。
* `TR`, `DROP_FIRST`, `KEEP_N`: 请根据你的扫描参数修改。
* `select_run(session, run_num)` 函数：修改逻辑以决定处理哪些数据。若要跑所有数据，请让该函数直接返回 `True`。

---

### 2. make_roi_fixedfx.py (固定效应与 ROI)
**功能**：合并同一 Session 下的多个 Run，并生成二值化 ROI。

**🧮 计算细节**：
* **合并算法**：逆方差加权 (Inverse Variance Weighted Fixed-effects)。
* **公式**：$\beta_{fixed} = \frac{\sum (w_i \cdot \beta_i)}{\sum w_i}$，其中权重 $w_i = 1 / Var_i$。
* **统计推断**：
    * **FDR**：Benjamini-Hochberg 过程，默认 $q=0.05$。
    * **聚类**：剔除小于 `MIN_CLUSTER` (默认 30) 个体素的孤立簇。

**⚠️ 必须修改的配置**：
* `CON_NAME`: 必须与第一步生成的对比度名称一致。
* `SIDE`: 检验方向 ("pos", "neg", 或 "two")。

---

### 3. vis.py (3D 体素可视化)
**功能**：生成交互式 HTML 文件，展示 Z 图、Beta 图、ROI Mask 和 Cluster Labels。
* **运行示例**：
    ```bash
    python3 vis.py --sub sub-01 --ses ses-01 --map all --open
    ```

---

### 4. vis_surface.py (皮层表面可视化)
**功能**：将 3D 体素数据投影到 FreeSurfer 的皮层表面并生成 HTML。

**🧮 计算细节**：
* 调用 FreeSurfer 的 `mri_vol2surf`。
* **连续变量 (Z/Beta)**：推荐使用 `frac-avg` 插值。
* **二值变量 (ROI)**：推荐使用 `frac-max` 配合 `nearest` 插值，确保边界清晰。

---

## 📝 输出文件说明
执行完所有步骤后，`roi_step/sub-XX/ses-XX/` 目录下将包含：
* `*_fixedfx_z.nii.gz`: 合并后的 Z 统计图。
* `*_fixedfx_beta.nii.gz`: 合并后的 Beta 估计图。
* `*_thr_fdrq0.05_...nii.gz`: 最终 ROI Mask (二值)。
* `*_labels.nii.gz`: ROI 簇标签图 (不同整数代表不同脑区)。
* `*_clusters.tsv`: 详细的簇坐标、大小和峰值信息表。
* `*.html`: 可视化报告文件。
