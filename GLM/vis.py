#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
from pathlib import Path
import webbrowser
from nilearn import plotting

# ===== 路径（与你现有结构保持一致）=====
ROOT    = Path("/media/xue/new_B1/glm")
BIDS    = ROOT / "BIDS_prep"
DERIV   = BIDS / "derivatives"           # 如果有 fmriprep，会放这里
ROI_DIR = ROOT / "roi_step"              # 你的 fixedfx / ROI 输出目录

def _glob_first(*patterns):
    """依次 glob 匹配 pattern，返回第一个存在的 Path 或 None。"""
    for pat in patterns:
        hits = sorted(glob.glob(str(pat)))
        if hits:
            return Path(hits[0])
    return None

def find_t1(sub: str, ses: str):
    """
    优先级：
    1) derivatives/fmriprep/sub/ses/anat/sub_ses_desc-preproc_T1w.nii*
    2) derivatives/fmriprep/sub/anat/sub_desc-preproc_T1w.nii*
    3) BIDS/sub/ses/anat/sub_ses_*T1w.nii*
    4) BIDS/sub/anat/sub_*T1w.nii*
    5) BIDS/sub/ses/func/*_space-T1w_boldref.nii*
    6) DERIV/**/.../func/*_space-T1w_boldref.nii*
    找不到就返回 None（后面会 fallback 到功能图本身）
    """
    cand = _glob_first(
        DERIV / "fmriprep" / sub / ses / "anat" / f"{sub}_{ses}_desc-preproc_T1w.nii*",
        DERIV / "fmriprep" / sub / "anat" / f"{sub}_desc-preproc_T1w.nii*",
    )
    if cand:
        print(f"[INFO] T1 背景（fMRIPrep）：{cand}")
        return cand

    cand = _glob_first(
        BIDS / sub / ses / "anat" / f"{sub}_{ses}_*T1w.nii*",
        BIDS / sub / "anat" / f"{sub}_*T1w.nii*",
    )
    if cand:
        print(f"[INFO] T1 背景（BIDS 原始/回写）：{cand}")
        return cand

    cand = _glob_first(
        BIDS / sub / ses / "func" / f"{sub}_{ses}_task-*_run-*_space-T1w_boldref.nii*",
        BIDS / sub / ses / "func" / f"{sub}_{ses}_*space-T1w_boldref.nii*",
    )
    if cand:
        print(f"[INFO] 背景回退为 space-T1w_boldref（3D）：{cand}")
        return cand

    cand = _glob_first(
        DERIV / "**" / sub / ses / "func" / f"{sub}_{ses}_*space-T1w_boldref.nii*",
        DERIV / "**" / sub / "func" / f"{sub}_*space-T1w_boldref.nii*",
    )
    if cand:
        print(f"[INFO] 背景回退为 derivatives/**/space-T1w_boldref：{cand}")
        return cand

    print("[WARN] 未找到 T1 或 space-T1w_boldref，将在后续用统计图自身作为背景。")
    return None

def find_fixedfx_z(sub: str, ses: str):
    """固定效应 Z 图（主效果 / 对比的方差加权合并）"""
    return ROI_DIR / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_z.nii.gz"

def find_fixedfx_beta(sub: str, ses: str):
    """固定效应 beta 图（效应大小估计）"""
    return ROI_DIR / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_beta.nii.gz"

def find_roi_mask(sub: str, ses: str):
    """
    二值 ROI 掩膜，一般是做了阈值+簇大小校正后的结果
    例子命名：sub-01_ses-01_con-numerosity_fixedfx_thr_fdrq0.05_zvoxXX_kYY.nii.gz
    我们用通配来匹配
    """
    return _glob_first(
        ROI_DIR / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_thr_*_k*.nii*"
    )

def find_roi_labels(sub: str, ses: str):
    """
    连通域标签图：把 ROI 各个簇分成 1,2,3,... 不同编号，方便区分多个cluster。
    假设命名里带 `_labels`。
    """
    return _glob_first(
        ROI_DIR / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_*labels.nii*",
        ROI_DIR / sub / ses / f"{sub}_{ses}_*labels.nii*",
    )

def view_and_save(
    stat_img,
    bg_img,
    cmap,
    out_html,
    title,
    vmin=None,
    vmax=None,
    threshold=0.0,
    black_bg=False,
    colorbar=True
):
    """
    用 nilearn.view_img 生成交互式可视化 (BrainSprite-like viewer)，并导出 HTML.
    """
    v = plotting.view_img(
        stat_map_img=str(stat_img),
        bg_img=str(bg_img) if bg_img is not None else None,
        cmap=cmap,
        colorbar=colorbar,
        vmin=vmin,
        vmax=vmax,
        threshold=threshold,
        black_bg=black_bg,
        title=title,
    )
    v.save_as_html(str(out_html))
    print(f"[OK] 保存：{out_html}")
    return out_html

def make_z_view(sub, ses, out_dir, bg_fallback):
    zfx = find_fixedfx_z(sub, ses)
    if not zfx.exists():
        print(f"[ERR] 找不到 fixedfx Z 图：{zfx}")
        return None
    bg = bg_fallback if bg_fallback is not None else zfx
    out_html = out_dir / f"{sub}_{ses}_view_z.html"
    print(f"[INFO] 使用 Z 图：{zfx}")
    return view_and_save(
        stat_img=zfx,
        bg_img=bg,
        cmap="cold_hot",
        vmin=None,
        vmax=None,
        threshold=0.0,
        title=f"{sub} {ses} fixedfx Z",
        out_html=out_html,
        colorbar=True,
        black_bg=False
    )

def make_beta_view(sub, ses, out_dir, bg_fallback):
    beta = find_fixedfx_beta(sub, ses)
    if not beta.exists():
        print(f"[WARN] 找不到 fixedfx beta 图：{beta}（跳过 beta）")
        return None
    bg = bg_fallback if bg_fallback is not None else beta
    out_html = out_dir / f"{sub}_{ses}_view_beta.html"
    print(f"[INFO] 使用 beta 图：{beta}")
    return view_and_save(
        stat_img=beta,
        bg_img=bg,
        cmap="cold_hot",
        vmin=None,
        vmax=None,
        threshold=0.0,
        title=f"{sub} {ses} fixedfx beta",
        out_html=out_html,
        colorbar=True,
        black_bg=False
    )

def make_roi_view(sub, ses, out_dir, bg_fallback):
    roimask = find_roi_mask(sub, ses)
    if roimask is None:
        print(f"[WARN] 找不到 ROI 掩膜（二值）：{out_dir}/*thr_*_k*.nii* （跳过 roi）")
        return None
    # 背景优先 T1/boldref，再退 Z
    zfx = find_fixedfx_z(sub, ses)
    bg = bg_fallback if bg_fallback is not None else (zfx if zfx.exists() else roimask)
    out_html = out_dir / f"{sub}_{ses}_view_roi.html"
    print(f"[INFO] 使用 ROI 掩膜：{roimask}")
    return view_and_save(
        stat_img=roimask,
        bg_img=bg,
        cmap="autumn",          # 单色调，突出 1
        vmin=0.0,
        vmax=1.0,
        threshold=0.5,          # 仅显示值为1的voxel/cluster
        title=f"{sub} {ses} ROI (FDR/cluster)",
        out_html=out_html,
        colorbar=True,          # 仍然给colorbar，方便看数值范围
        black_bg=False
    )

def make_labels_view(sub, ses, out_dir, bg_fallback):
    labels = find_roi_labels(sub, ses)
    if labels is None:
        print(f"[WARN] 找不到 ROI labels（连通域标签）：{out_dir}/*labels.nii* （跳过 labels）")
        return None
    # 同样用 T1/boldref 或者退 Z
    zfx = find_fixedfx_z(sub, ses)
    bg = bg_fallback if bg_fallback is not None else (zfx if zfx.exists() else labels)
    out_html = out_dir / f"{sub}_{ses}_view_labels.html"
    print(f"[INFO] 使用 ROI labels：{labels}")
    # 这里用多类别调色板；colorbar 也保留
    return view_and_save(
        stat_img=labels,
        bg_img=bg,
        cmap="tab20",
        vmin=None,
        vmax=None,
        threshold=0.0,
        title=f"{sub} {ses} ROI labels (clusters)",
        out_html=out_html,
        colorbar=True,
        black_bg=False
    )

def main():
    ap = argparse.ArgumentParser(
        description="体素空间交互可视化（Z/beta/ROI/labels -> HTML）"
    )
    ap.add_argument("--sub", required=True, help="如 sub-01")
    ap.add_argument("--ses", required=True, help="如 ses-01")
    ap.add_argument("--map",
                    choices=["z", "beta", "roi", "labels", "all"],
                    default="all",
                    help="画哪一种图；all=全部能找到的都画")
    ap.add_argument("--open", action="store_true",
                    help="保存后自动用浏览器打开第一个生成的HTML")
    args = ap.parse_args()

    sub, ses, which = args.sub, args.ses, args.map

    # 先准备输出目录
    out_dir = ROI_DIR / sub / ses
    out_dir.mkdir(parents=True, exist_ok=True)

    # 预先决定背景（T1 或 space-T1w_boldref）
    bg_fallback = find_t1(sub, ses)

    opened_any = None

    if which in ("z", "all"):
        opened_any = make_z_view(sub, ses, out_dir, bg_fallback) or opened_any

    if which in ("beta", "all"):
        opened_any = make_beta_view(sub, ses, out_dir, bg_fallback) or opened_any

    if which in ("roi", "all"):
        opened_any = make_roi_view(sub, ses, out_dir, bg_fallback) or opened_any

    if which in ("labels", "all"):
        opened_any = make_labels_view(sub, ses, out_dir, bg_fallback) or opened_any

    if args.open and opened_any is not None:
        webbrowser.open(f"file://{opened_any}")

if __name__ == "__main__":
    main()

