#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
from pathlib import Path
import subprocess as sp

import numpy as np
import nibabel as nib
from nilearn import plotting, surface

# ========= 路径配置 =========
ROOT = Path("/media/xue/new_B1/glm")
ROI_STEP = ROOT / "roi_step"

# ========= 脚本内可调参数 =========
# proj 可选字符串（会被翻译成 mri_vol2surf 的选项）：
#   "frac:0.5"            -> --projfrac 0.5               （按比例；0=white, 1=pial）
#   "frac-avg:0:1:0.1"    -> --projfrac-avg 0 1 0.1       （区间平均）
#   "frac-max:0:1:0.1"    -> --projfrac-max 0 1 0.1       （区间最大，适合二值 ROI）
#   "dist:1.5"            -> --projdist 1.5               （按毫米）
#   "dist-avg:-1:2:0.5"   -> --projdist-avg -1 2 0.5
#   "dist-max:-1:2:0.5"   -> --projdist-max -1 2 0.5
proj = "frac-max:0:1:0.1"  # 连续图推荐；二值 ROI 推荐 "frac-max:0:1:0.1"

# interp 可选："nearest" / "trilinear"（部分版本支持 "cubic"）
# 连续图（Z/beta）建议 "trilinear"；二值 ROI 建议 "nearest"
interp = "trilinear"

def build_proj_args(proj_str: str):
    if proj_str.startswith("frac:"):
        val = proj_str.split(":", 1)[1]
        return ["--projfrac", str(float(val))]
    if proj_str.startswith("frac-avg:"):
        a, b, step = proj_str.split(":", 1)[1].split(":")
        return ["--projfrac-avg", str(float(a)), str(float(b)), str(float(step))]
    if proj_str.startswith("frac-max:"):
        a, b, step = proj_str.split(":", 1)[1].split(":")
        return ["--projfrac-max", str(float(a)), str(float(b)), str(float(step))]
    if proj_str.startswith("dist:"):
        val = proj_str.split(":", 1)[1]
        return ["--projdist", str(float(val))]
    if proj_str.startswith("dist-avg:"):
        a, b, step = proj_str.split(":", 1)[1].split(":")
        return ["--projdist-avg", str(float(a)), str(float(b)), str(float(step))]
    if proj_str.startswith("dist-max:"):
        a, b, step = proj_str.split(":", 1)[1].split(":")
        return ["--projdist-max", str(float(a)), str(float(b)), str(float(step))]
    raise ValueError(f"无法解析 proj='{proj_str}'，请按注释格式填写。")

def find_first(*patterns):
    for pat in patterns:
        hits = sorted(glob.glob(str(pat)))
        if hits:
            return Path(hits[0]).resolve()
    return None

def pick_subjects_dir():
    cands = [
        ROOT / "BIDS_prep" / "sourcedata" / "freesurfer",   # 优先：完整 recon-all
        ROOT / "BIDS_prep" / "derivatives" / "freesurfer",  # 其次：只有表面
    ]
    for d in cands:
        if d.is_dir():
            return d
    return None

def map_vol_to_surf(vol_path: Path, subj: str, hemi_fs: str, out_mgh: Path):
    """调用 mri_vol2surf 把体素图投到表面（注意参数顺序，避免 'Option 0.0 unknown'）"""
    proj_args = build_proj_args(proj)
    args = [
        "mri_vol2surf",
        "--mov", str(vol_path),
        *proj_args,                 # 这里紧跟在 --mov 之后
        "--regheader", subj,
        "--hemi", hemi_fs,
        "--interp", interp,
        "--o", str(out_mgh),
    ]
    print("[CMD]", " ".join(args))
    sp.run(args, check=True)

def load_surf_mesh_and_bg(subj_dir: Path, hemi_fs: str, which: str):
    """加载 FreeSurfer 网格和 sulc 作为背景底纹"""
    surf_file = subj_dir / "surf" / f"{hemi_fs}.{which}"
    if not surf_file.exists():
        raise FileNotFoundError(f"找不到表面文件：{surf_file}")
    mesh = surface.load_surf_mesh(str(surf_file))
    sulc_file = subj_dir / "surf" / f"{hemi_fs}.sulc"
    bg = surface.load_surf_data(str(sulc_file)) if sulc_file.exists() else None
    return mesh, bg

def main():
    ap = argparse.ArgumentParser(description="把个体 fixedfx 结果/ROI 投到 FreeSurfer 表面并交互查看（保留 colorbar）")
    ap.add_argument("--sub", required=True, help="如 sub-01")
    ap.add_argument("--ses", required=True, help="如 ses-01")
    ap.add_argument("--map", choices=["z", "beta", "roi"], default="z",
                    help="选择要显示的图：z / beta / roi（二值掩膜）")
    ap.add_argument("--hemi", choices=["L", "R", "both"], default="both",
                    help="半球：L/R/both")
    ap.add_argument("--open", action="store_true", help="保存后自动用浏览器打开")
    ap.add_argument("--surfkind", default="inflated",
                    help="显示的表面几何：white/pial/midthickness/inflated（默认 inflated）")
    args = ap.parse_args()

    sdir = pick_subjects_dir()
    if sdir is None:
        print("[ERR] 找不到 FreeSurfer SUBJECTS_DIR（BIDS_prep/sourcedata 或 derivatives 下都没有 freesurfer）")
        sys.exit(1)
    os.environ["SUBJECTS_DIR"] = str(sdir)
    print(f"[INFO] SUBJECTS_DIR = {sdir}")

    sub, ses = args.sub, args.ses
    subj_dir = sdir / sub

    # 选择体素图 & 预设配色
    if args.map == "z":
        vol = ROI_STEP / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_z.nii.gz"
        if not vol.exists():
            print(f"[ERR] 找不到 Z 图：{vol}"); sys.exit(1)
        print(f"[INFO] 使用 Z 图：{vol}")
        cmap = "cold_hot"      # 连续图：冷热色
        is_roi = False
    elif args.map == "beta":
        vol = ROI_STEP / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_beta.nii.gz"
        if not vol.exists():
            print(f"[ERR] 找不到 beta 图：{vol}"); sys.exit(1)
        print(f"[INFO] 使用 beta 图：{vol}")
        cmap = "cold_hot"
        is_roi = False
    else:
        vol = find_first(ROI_STEP / sub / ses / f"{sub}_{ses}_con-numerosity_fixedfx_thr_fdrq0.05*_k*.nii.gz")
        if vol is None:
            print(f"[ERR] 找不到 ROI 掩膜：{ROI_STEP}/{sub}/{ses} 下 *_thr_fdrq0.05*_k*.nii.gz")
            sys.exit(1)
        print(f"[INFO] 使用 ROI 体素图：{vol}")
        cmap = "autumn"        # ROI：橙黄
        is_roi = True
        if interp != "nearest":
            print(f"[WARN] ROI 建议 interp='nearest' 以保持0/1。当前为 '{interp}'。")

    hemis = ["L", "R"] if args.hemi == "both" else [args.hemi]
    HEMIS = {"L": "lh", "R": "rh"}

    html_paths = []
    for H in hemis:
        hemi_fs = HEMIS[H]
        out_mgh = ROI_STEP / sub / ses / (f"{hemi_fs}.{'roi' if is_roi else args.map}.mgh")
        out_mgh.parent.mkdir(parents=True, exist_ok=True)

        # 1) 体素 -> 表面
        map_vol_to_surf(vol, sub, hemi_fs, out_mgh)

        # 2) 载网格 + sulc 当底纹
        mesh, bg = load_surf_mesh_and_bg(subj_dir, hemi_fs, which=args.surfkind)

        # 3) 载到顶点纹理
        data = np.squeeze(nib.load(str(out_mgh)).get_fdata())
        if data.ndim > 1:
            data = np.squeeze(data)
        if data.ndim != 1:
            data = data.reshape((-1,))

        # 4) 可视化参数
        if is_roi:
            # 二值 ROI：0 设为 NaN（透明），并固定 vmin/vmax，以便 colorbar 有 0-1 刻度
            n_pos = int((data > 0).sum())
            print(f"[INFO] {hemi_fs} 半球 ROI 顶点数：{n_pos}")
            if n_pos == 0:
                print(f"[WARN] {hemi_fs} 半球投影后 ROI 为空，页面只会显示背景皮层。")
            data = data.astype(float)
            data[data <= 0] = np.nan
            vmin, vmax = 0.0, 1.0
            threshold = None
        else:
            vmax = float(np.nanpercentile(np.abs(data), 99)) if np.any(np.isfinite(data)) else 1.0
            vmin = -vmax
            threshold = 0.0

        title = f"{sub} {ses} {args.map} ({hemi_fs}, proj='{proj}', interp='{interp}', surf='{args.surfkind}')"

        v = plotting.view_surf(
            surf_mesh=mesh,
            surf_map=data,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            bg_map=bg,          # 叠 sulc 作底纹
            bg_on_data=True,
            darkness=0.6,       # 控制底纹对统计颜色的影响（0=不影响，1=很暗）
            threshold=threshold,
            symmetric_cmap=False,
            colorbar=True,      # <<< 保留 colorbar 显示
            title=title,
        )

        html_out = ROI_STEP / sub / ses / f"vis_{args.map}_{hemi_fs}.html"
        v.save_as_html(str(html_out))
        print(f"[OK] 保存交互视图：{html_out}")
        html_paths.append(html_out)

    if args.open and html_paths:
        import webbrowser
        for hp in html_paths:
            webbrowser.open(f"file://{hp}")

if __name__ == "__main__":
    main()

