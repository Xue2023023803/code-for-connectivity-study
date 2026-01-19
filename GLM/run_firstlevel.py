#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, re, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from nilearn import image
from nilearn.glm.first_level import FirstLevelModel, make_first_level_design_matrix

# ========== 配置 ==========
ROOT = Path("/media/xue/new_B1/glm")
BIDS = ROOT / "BIDS_prep"
EVENTS_FILE = ROOT / "events_output" / "events_bids.tsv"
OUTDIR = ROOT / "first_step"          # 输出目录
TR = 1.95                             # 固定 TR
SESSIONS = ["ses-01", "ses-02"]

# BOLD 裁切策略：去掉前6个TR，只保留176个TR（不改原始文件，只在内存中）
DROP_FIRST = 6
KEEP_N = 176

# 可用就加；不存在就跳过
CONFOUND_COLS = [
    "trans_x","trans_y","trans_z","rot_x","rot_y","rot_z",
    "csf","white_matter","global_signal","framewise_displacement"
]

# ========== 小工具 ==========
def select_run(session: str, run_num: int) -> bool:
    """ses-01 选偶数；ses-02 选奇数"""
    if session == "ses-01":
        return run_num % 2 == 0
    if session == "ses-02":
        return run_num % 2 == 1
    return False

def prefer_gz(p: Path) -> Path:
    """当 .nii 与 .nii.gz 并存时，优先 .nii.gz"""
    if p.suffix == ".gz":
        return p
    gz = Path(str(p) + ".gz")
    return gz if gz.exists() else p

def find_files_for_run(bold_nii: Path):
    """
    给定 *_space-T1w_desc-preproc_bold.nii.gz，返回同 run 的 mask、confounds
    confounds 兼容两种命名：有/无 _space-T1w
    """
    base = bold_nii.name.replace("_space-T1w_desc-preproc_bold.nii.gz", "")
    func_dir = bold_nii.parent

    mask = func_dir / f"{base}_space-T1w_desc-brain_mask.nii.gz"
    mask = prefer_gz(mask)

    # 你的目录里更常见是不带 space-T1w 的 confounds：
    conf = func_dir / f"{base.replace('_space-T1w', '')}_desc-confounds_timeseries.tsv"
    if not conf.exists():
        conf2 = func_dir / f"{base}_desc-confounds_timeseries.tsv"
        if conf2.exists():
            conf = conf2

    return mask, conf

def check_events_compatible(evt_df: pd.DataFrame, n_scans: int, tr: float, tol_ratio: float = 0.25):
    """
    不修改 events，只检查它与裁切后的 BOLD 是否一致：
      - onset >= 0 - tol, offset <= totalT + tol
      - duration > 0
    仅在不一致时给出 WARNING，不中断。
    """
    required = {"onset", "duration"}
    if not required.issubset(evt_df.columns):
        warnings.warn("events 文件缺少必须列：onset/duration（仅告警）")
        return

    tol = tol_ratio * tr
    totalT = n_scans * tr
    onset = evt_df["onset"].to_numpy(dtype=float)
    dur   = evt_df["duration"].to_numpy(dtype=float)
    offset = onset + dur

    if np.any(dur <= 0):
        k = int(np.sum(dur <= 0))
        warnings.warn(f"[events 校验] 发现 {k} 个 duration <= 0 的条目（仅告警）")

    k1 = int(np.sum(onset < -tol))
    k2 = int(np.sum((onset >= -tol) & (onset < 0)))
    k3 = int(np.sum(offset > totalT + tol))
    k4 = int(np.sum((offset > totalT) & (offset <= totalT + tol)))

    if k1 or k3 or k2 or k4:
        msg = "[events 校验] 与裁切后时间轴可能不完全一致："
        if k1: msg += f" onset<-{tol:.3f}s 的有 {k1}；"
        if k2: msg += f" onset ∈ [-{tol:.3f},0) 的有 {k2}（视作浮点误差）；"
        if k3: msg += f" offset>totalT+{tol:.3f}s 的有 {k3}；"
        if k4: msg += f" offset ∈ (totalT,totalT+{tol:.3f}] 的有 {k4}（视作浮点误差）。"
        warnings.warn(msg)

def build_design(frame_times: np.ndarray, events_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """统一 events，只取 numerosity（baseline 为隐含基线），构建设计矩阵；返回 (design, 用到的events)"""
    ev = pd.read_csv(events_file, sep="\t")
    if "trial_type" not in ev.columns:
        raise RuntimeError(f"事件文件缺少 trial_type 列: {events_file}")
    ev_num = ev[ev["trial_type"].str.lower() == "numerosity"].copy()
    ev_num["trial_type"] = "numerosity"

    design = make_first_level_design_matrix(
        frame_times,
        events=ev_num,
        hrf_model="glover",
        drift_model="cosine",
        high_pass=1/128.0
    )
    return design, ev_num

def load_confounds(conf_path: Path, n_rows: int) -> pd.DataFrame:
    """读取混杂项并对齐设计矩阵行数；不存在则返回空表"""
    if not conf_path.exists():
        return pd.DataFrame(index=range(n_rows))
    conf = pd.read_csv(conf_path, sep="\t").fillna(0.0)
    keep = [c for c in CONFOUND_COLS if c in conf.columns]
    conf = conf[keep].reset_index(drop=True)
    # 对齐长度（此处以裁切后的 KEEP_N 为准）
    if len(conf) >= n_rows:
        conf = conf.iloc[:n_rows, :].copy()
    else:
        extra = pd.DataFrame(0.0, index=range(n_rows - len(conf)), columns=conf.columns)
        conf = pd.concat([conf, extra], axis=0, ignore_index=True)
    return conf

def save_design_png(design: pd.DataFrame, out_png: Path):
    """保存设计矩阵热图（PNG）"""
    import matplotlib.pyplot as plt
    plt.figure(figsize=(max(8, 0.35*len(design.columns)), 6))
    plt.imshow(design.values, aspect='auto', interpolation='nearest')
    plt.colorbar(shrink=0.6)
    plt.xticks(range(len(design.columns)), design.columns, rotation=90, fontsize=6)
    plt.xlabel("Regressors")
    plt.ylabel("Time (scans)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def try_save_residuals(model: FirstLevelModel, out_file: Path) -> bool:
    """
    尝试保存 4D 残差图：不同 nilearn 版本实现不同，这里做了尽力而为的兼容。
    成功返回 True；失败返回 False（仅告警，不中断主流程）。
    """
    try:
        # 1) 常见：residuals 是 list[Niimg] 或 list[np.ndarray]
        resid_attr = getattr(model, "residuals", None)
        if resid_attr is None:
            resid_attr = getattr(model, "residuals_", None)

        if resid_attr is not None:
            # list 情况
            if isinstance(resid_attr, (list, tuple)) and len(resid_attr) >= 1:
                r0 = resid_attr[0]
                # 若已是 Nifti 图像
                if hasattr(r0, "to_filename"):
                    r0.to_filename(str(out_file))
                    return True
                # 若是 2D 时空矩阵（T x V），回投到 4D
                if hasattr(model, "masker_"):
                    from nilearn.image import new_img_like
                    res_img = model.masker_.inverse_transform(r0)
                    res_img.to_filename(str(out_file))
                    return True
            # 单对象情况
            if hasattr(resid_attr, "to_filename"):
                resid_attr.to_filename(str(out_file))
                return True

        # 2) 兜底：通过 results_ 取各 label 的 residuals（部分版本可用）
        labels = getattr(model, "labels_", None)
        results = getattr(model, "results_", None)
        if labels is not None and results is not None and hasattr(model, "masker_"):
            # results 可能是 list 或 dict，尽量用整数 label 访问
            uniq = np.unique(labels)
            # 预取形状
            # 用任意一个 label 的 residuals 来推断时间维度和该 parcel 的体素数
            any_lab = [lab for lab in uniq if lab != -1]
            if any_lab:
                lab0 = int(any_lab[0]) if not isinstance(any_lab[0], (int, np.integer)) else any_lab[0]
                res0 = results[lab0].residuals  # (T x n_vox_in_lab)
                T = res0.shape[0]
                V = len(labels)
                R = np.zeros((T, V), dtype=res0.dtype)
                for lab in uniq:
                    if lab == -1:  # 背景
                        continue
                    ilab = int(lab) if not isinstance(lab, (int, np.integer)) else lab
                    idx = np.where(labels == lab)[0]
                    R[:, idx] = results[ilab].residuals
                res_img = model.masker_.inverse_transform(R)
                res_img.to_filename(str(out_file))
                return True
    except Exception as e:
        warnings.warn(f"[WARN] 保存残差失败：{e}")
    return False

# ========== 主流程 ==========
def main():
    assert EVENTS_FILE.exists(), f"事件文件不存在: {EVENTS_FILE}"
    subjects = sorted([p for p in BIDS.iterdir() if p.is_dir() and p.name.startswith("sub-")])

    out_recs = []

    for sub_dir in subjects:
        sub = sub_dir.name
        for ses in SESSIONS:
            func_dir = sub_dir / ses / "func"
            if not func_dir.exists():
                continue

            # 自动解析 task（用 task-* 通配）
            bolds = sorted(func_dir.glob(f"{sub}_{ses}_task-*_run-*_space-T1w_desc-preproc_bold.nii.gz"))
            if not bolds:
                continue

            for bold in bolds:
                m_run  = re.search(r"run-(\d{3})", bold.name)
                m_task = re.search(r"task-([A-Za-z0-9]+)", bold.name)
                if not m_run:
                    continue
                run_num = int(m_run.group(1))
                task = m_task.group(1) if m_task else "task"

                if not select_run(ses, run_num):
                    continue

                mask, conf = find_files_for_run(bold)

                # 原始 BOLD
                img = image.load_img(bold)
                n_scans = img.shape[-1]

                # 检查是否足够裁切
                need = DROP_FIRST + KEEP_N
                if n_scans < need:
                    warnings.warn(f"{bold.name} 的体积数 {n_scans} < 需要的 {need}（前{DROP_FIRST}丢弃+保留{KEEP_N}），已跳过。")
                    continue

                # === 仅在内存中裁切，不写回原数据 ===
                bold_trim = image.index_img(img, slice(DROP_FIRST, DROP_FIRST + KEEP_N))
                frame_times = np.arange(KEEP_N) * TR

                # 设计矩阵 + 混杂
                design, ev_used = build_design(frame_times, EVENTS_FILE)
                # 一致性校验（仅告警，不中断）
                check_events_compatible(ev_used, n_scans=KEEP_N, tr=TR, tol_ratio=0.25)

                conf_df = load_confounds(conf, n_rows=len(design))
                design = pd.concat([design.reset_index(drop=True),
                                    conf_df.reset_index(drop=True)], axis=1)

                # 拟合（提供 design_matrices 时，t_r 被忽略，属正常提示）
                model = FirstLevelModel(
                    t_r=TR,
                    mask_img=str(mask) if mask.exists() else None,
                    noise_model="ar1",
                    standardize=False,
                    minimize_memory=False
                )
                model = model.fit(bold_trim, design_matrices=design)

                # 对比：numerosity（vs baseline 隐含）
                # 支持列名或向量形式；优先列名（更直观）
                if "numerosity" not in design.columns:
                    raise RuntimeError("设计矩阵中未找到 'numerosity' 列，请检查 events 的 trial_type 命名。")
                zmap = model.compute_contrast("numerosity", output_type="z_score")
                beta = model.compute_contrast("numerosity", output_type="effect_size")

                # 输出路径
                out_dir = OUTDIR / sub / ses / f"run-{run_num:03d}"
                out_dir.mkdir(parents=True, exist_ok=True)
                prefix = f"{sub}_{ses}_task-{task}_run-{run_num:03d}"

                zmap_file    = out_dir / f"{prefix}_con-numerosity_z.nii.gz"
                beta_file    = out_dir / f"{prefix}_con-numerosity_beta.nii.gz"
                resid_file   = out_dir / f"{prefix}_residuals.nii.gz"     # 4D 残差
                design_tsv   = out_dir / f"{prefix}_design_matrix.tsv"
                design_png   = out_dir / f"{prefix}_design.png"
                meta_file    = out_dir / f"{prefix}_meta.json"

                zmap.to_filename(str(zmap_file))
                beta.to_filename(str(beta_file))
                design.to_csv(design_tsv, sep="\t", index=False)
                save_design_png(design, design_png)

                # 残差图（尽力保存，失败仅告警）
                saved_resid = try_save_residuals(model, resid_file)
                if not saved_resid:
                    warnings.warn(f"[WARN] {prefix}: 残差图未保存（已跳过）。")

                # 记录元数据
                meta = {
                    "subject": sub, "session": ses, "run": f"{run_num:03d}", "task": task,
                    "bold": str(bold),
                    "mask": str(mask) if mask.exists() else None,
                    "confounds": str(conf) if conf.exists() else None,
                    "events": str(EVENTS_FILE),
                    "TR": TR, "n_scans_original": int(n_scans),
                    "drop_first": DROP_FIRST, "keep_n": KEEP_N,
                    "contrast": "numerosity (vs implicit baseline)",
                    "hrf_model": "glover", "drift_model": "cosine", "high_pass": 1/128.0
                }
                meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")

                out_recs.append([sub, ses, f"{run_num:03d}", task, str(zmap_file), str(beta_file),
                                 str(resid_file if saved_resid else "NA")])
                print(f"[OK] {sub} {ses} run-{run_num:03d} -> z:{zmap_file.name}  beta:{beta_file.name}  resid:{resid_file.name if saved_resid else 'NA'}")

    # 汇总
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if out_recs:
        summary = pd.DataFrame(out_recs, columns=["subject","session","run","task","zmap","beta","residuals"])
        summary_path = OUTDIR / "firstlevel_outputs.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"\n[SUMMARY] Saved: {summary_path}")
    else:
        print("[WARNING] 没有发现符合筛选规则的 run；请检查目录与命名。")

if __name__ == "__main__":
    main()

