#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
第二步：被试内（会话内）固定效应合并 + 只在灰质做 FDR(单侧) + 体素阈值 + 最小簇大小
- 读取第一步的 run 级对比图（beta / z；若存在 var 也会用）
- 逆方差加权合并为被试×会话的 fixed-effects beta/var/z
- 构建“功能灰质 mask”：功能脑掩膜 ∧ (GM probseg ≥ 阈值)；全部重采样到合并后栅格
- 在该 mask 内做单侧 FDR(q=0.05)，再叠加体素 Z 阈值（可选）与最小簇大小
- 导出：fixed beta/var/z、阈后 ROI mask、簇标注图、簇信息 TSV

依赖：numpy, pandas, nibabel, nilearn, scipy（ndimage）均为常规神经影像栈
"""

import re, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import image
from scipy import ndimage

# ========== 配置 ==========
ROOT          = Path("/media/xue/new_B1/glm")
BIDS          = ROOT / "BIDS_prep"           # 与第一步一致
FIRSTLEVEL    = ROOT / "first_step"          # 第一步输出目录
OUTDIR        = ROOT / "roi_step"            # 本步输出目录
SESSIONS      = ["ses-01", "ses-02"]         # 会话列表
ROI_SESSIONS  = ["ses-01"]                   # 在这些会话上做 ROI（默认只用 ses-01）
CON_NAME      = "numerosity"                 # 对比名（与第一步一致）
GM_PROB_THR   = 0.5                          # 只在灰质：GM 概率阈值
USE_GM_ONLY   = True                         # 只在灰质做统计
SIDE          = "pos"                        # 单侧："pos" or "neg"；如要双侧用 "two"
Q_FDR         = 0.05                         # 体素层面 FDR q
Z_VOX         = 0                          # 额外体素 Z 阈值；0 代表不用
MIN_CLUSTER   = 30                           # 最小簇大小（体素数）

# ========== 小工具 ==========
def find_runs_for_subject_session(sub: str, ses: str):
    """罗列该被试该会话下第一步生成的 run 文件夹"""
    base = FIRSTLEVEL / sub / ses
    if not base.exists():
        return []
    return sorted([p for p in base.iterdir() if p.is_dir() and re.search(r"run-\d{3}$", p.name)])

def first_existing(paths):
    for p in paths:
        if p is not None and Path(p).exists():
            return Path(p)
    return None

def load_nii(path):
    img = nib.load(str(path))
    return img, img.get_fdata(dtype=np.float32), img.affine

def voxel_volume_mm3(affine):
    # 体素体积 = |det(affine[:3,:3])|
    return float(abs(np.linalg.det(affine[:3, :3])))

def safe_divide(a, b):
    out = np.zeros_like(a, dtype=np.float32)
    mask = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
    out[mask] = a[mask] / b[mask]
    return out

# ---- FDR（BH step-up）
def normal_cdf(z):
    # Φ(z) = 0.5 * (1 + erf(z / sqrt(2)))
    return 0.5 * (1.0 + math.erf(float(z) / math.sqrt(2.0)))

def z_to_p(z_arr, side="pos"):
    z = z_arr.astype(np.float64, copy=False)
    p = np.ones_like(z, dtype=np.float64)
    if side == "pos":
        # 单侧：H1: Z>0 → p = 1 - Φ(Z)
        p = 1.0 - 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
    elif side == "neg":
        # 单侧：H1: Z<0 → p = Φ(Z)
        p = 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
    elif side == "two":
        # 双侧：p = 2*min(Φ(Z), 1-Φ(Z))
        phi = 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
        p = 2.0 * np.minimum(phi, 1.0 - phi)
    else:
        raise ValueError("side must be 'pos', 'neg', or 'two'")
    # 数值安全
    p = np.clip(p, 1e-300, 1.0)
    return p

def erf_vec(x):
    # 向量化的 math.erf
    x = np.asarray(x, dtype=np.float64)
    return np.vectorize(math.erf, otypes=[np.float64])(x)

def fdr_bh_mask(pvals_in_mask, q=0.05):
    """输入：仅 mask 内的 p 值向量；输出：阈后布尔向量（同长度）与阈值 p*"""
    m = pvals_in_mask.size
    if m == 0:
        return np.zeros_like(pvals_in_mask, dtype=bool), 0.0
    order = np.argsort(pvals_in_mask)
    p_sorted = pvals_in_mask[order]
    # 找到最大的 k 使 p_(k) <= (k/m) q
    thresh = (np.arange(1, m+1) / float(m)) * q
    passed = p_sorted <= thresh
    if not np.any(passed):
        return np.zeros_like(pvals_in_mask, dtype=bool), 0.0
    k_star = np.where(passed)[0].max()
    p_star = p_sorted[k_star]
    keep = pvals_in_mask <= p_star
    return keep, float(p_star)

# ---- 构建“功能灰质 mask”：功能脑掩膜 ∧ (GM prob ≥ thr)，全部重采样到参考图（合并后 Z 图）上
def load_subject_masks_in_ref_grid(sub: str, ref_img: nib.Nifti1Image):
    # 功能脑掩膜：取该被试任意一张 *_desc-brain_mask（space-T1w）
    func_masks = sorted((BIDS / sub).glob("ses-*/func/*_space-T1w_desc-brain_mask.nii*"))
    func_mask_path = first_existing(func_masks)

    # 灰质概率图：优先 anat 下 *_label-GM_probseg（若有 space 标签，要求 space-T1w；避免误用 MNI）
    gm_candidates = []
    anat_dir = BIDS / sub / "anat"
    gm_candidates += list(anat_dir.glob(f"{sub}_label-GM_probseg.nii*"))
    gm_candidates += [p for p in anat_dir.glob(f"{sub}_space-T1w_label-GM_probseg.nii*")]
    # 排除非 T1w 的 space
    gm_candidates = [p for p in gm_candidates if ("space-" not in p.name) or ("space-T1w" in p.name)]
    gm_path = first_existing(sorted(gm_candidates))

    brain_mask = None
    gm_prob = None

    if func_mask_path is not None:
        bm = image.resample_to_img(str(func_mask_path), ref_img, interpolation="nearest")
        brain_mask = (bm.get_fdata(dtype=np.float32) > 0.5)

    if gm_path is not None:
        gm = image.resample_to_img(str(gm_path), ref_img, interpolation="continuous")
        gm_prob = gm.get_fdata(dtype=np.float32)
        gm_prob = np.clip(gm_prob, 0.0, 1.0)

    return brain_mask, gm_prob

# ---- 逆方差加权 fixed-effects 合并（逐 run 流式累加，省内存）
def fixedfx_combine_beta_var(run_items, ref_img=None):
    """
    run_items: 列表，每项为 { 'beta': Path, 'z': Path, 'var': Path or None }
    返回：fixed_beta, fixed_var, fixed_z 的 nib 图像（与参考栅格一致）
    """
    beta_ref_img, ref_affine, ref_shape = None, None, None

    sum_w = None         # Σ(1/var)
    sum_wbeta = None     # Σ(beta/var)

    for item in run_items:
        beta_img, beta_data, beta_aff = load_nii(item['beta'])
        if beta_ref_img is None:
            beta_ref_img = beta_img
            ref_affine = beta_aff
            ref_shape = beta_data.shape

        # 重采样（理论上第一步输出的 run 都在同一空间/分辨率，这里谨慎处理）
        if beta_data.shape != ref_shape or not np.allclose(beta_aff, ref_affine):
            beta_img = image.resample_img(beta_img, target_affine=ref_affine, target_shape=ref_shape, interpolation="continuous")
            beta_data = beta_img.get_fdata(dtype=np.float32)

        # var：优先读现成 var 文件；否则用 (beta/z)^2 反推
        var_data = None
        if item.get('var') is not None and Path(item['var']).exists():
            _, var_data, _ = load_nii(item['var'])
            if var_data.shape != ref_shape or not np.allclose(beta_aff, ref_affine):
                var_img = image.resample_img(nib.load(str(item['var'])), target_affine=ref_affine, target_shape=ref_shape, interpolation="continuous")
                var_data = var_img.get_fdata(dtype=np.float32)

        if var_data is None:
            _, z_data, _ = load_nii(item['z'])
            if z_data.shape != ref_shape or not np.allclose(beta_aff, ref_affine):
                z_img = image.resample_img(nib.load(str(item['z'])), target_affine=ref_affine, target_shape=ref_shape, interpolation="continuous")
                z_data = z_img.get_fdata(dtype=np.float32)
            with np.errstate(divide='ignore', invalid='ignore'):
                var_data = np.square(safe_divide(beta_data, z_data))
                # z≈0 的地方会得到 inf/NaN → 后面处理成权重 0

        # 权重：w = 1/var；对非正/非有限 var 置 0
        with np.errstate(divide='ignore', invalid='ignore'):
            w = 1.0 / var_data
        w[~np.isfinite(w)] = 0.0
        w[w < 0] = 0.0

        if sum_w is None:
            sum_w = np.zeros(ref_shape, dtype=np.float32)
            sum_wbeta = np.zeros(ref_shape, dtype=np.float32)

        sum_w += w
        sum_wbeta += (w * beta_data)

    # 计算 fixed beta/var/z
    beta_fixed = np.zeros(ref_shape, dtype=np.float32)
    var_fixed = np.full(ref_shape, np.inf, dtype=np.float32)
    mask = sum_w > 0
    beta_fixed[mask] = sum_wbeta[mask] / sum_w[mask]
    var_fixed[mask]  = 1.0 / sum_w[mask]
    z_fixed = np.zeros(ref_shape, dtype=np.float32)
    good = mask & np.isfinite(var_fixed) & (var_fixed > 0)
    z_fixed[good] = beta_fixed[good] / np.sqrt(var_fixed[good])

    # 打包成 NIfTI
    beta_img = nib.Nifti1Image(beta_fixed, ref_affine, header=beta_ref_img.header)
    var_img  = nib.Nifti1Image(var_fixed,  ref_affine, header=beta_ref_img.header)
    z_img    = nib.Nifti1Image(z_fixed,    ref_affine, header=beta_ref_img.header)
    return beta_img, var_img, z_img, mask

# ---- 阈值与聚类
def threshold_and_cluster(z_img, base_mask, side="pos", q=0.05, z_vox=0.0, min_cluster=30):
    z = z_img.get_fdata(dtype=np.float32)
    mask = (base_mask.astype(bool)) & np.isfinite(z)

    # 单/双侧 p 值
    if side == "pos":
        p = 1.0 - 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
        sign_mask = z > 0
    elif side == "neg":
        p = 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
        sign_mask = z < 0
    elif side == "two":
        phi = 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
        p = 2.0 * np.minimum(phi, 1.0 - phi)
        sign_mask = np.ones_like(z, dtype=bool)  # 双侧不看符号
    else:
        raise ValueError("side must be 'pos', 'neg', or 'two'")

    # 仅在 mask 内做 FDR
    p_in = p[mask]
    keep_vec, p_star = fdr_bh_mask(p_in, q=q)

    fdr_map = np.zeros_like(mask, dtype=bool)
    fdr_map[mask] = keep_vec

    # 额外体素 Z 阈值（单侧：z>=阈；双侧：|z|>=阈）
    if z_vox > 0:
        if side == "two":
            vox_mask = np.abs(z) >= z_vox
        else:
            vox_mask = z >= z_vox if side == "pos" else (-z >= z_vox)
        fdr_map &= vox_mask

    # 单侧时还要求符号一致
    fdr_map &= sign_mask

    # 最小簇大小
    labeled, nlab = ndimage.label(fdr_map.astype(np.int8), structure=np.ones((3,3,3), dtype=np.int8))
    if nlab == 0:
        return fdr_map, labeled, 0, 0.0, []

    sizes = ndimage.sum(np.ones_like(fdr_map, dtype=np.int32), labeled, index=np.arange(1, nlab+1))
    keep_labels = set([i+1 for i,s in enumerate(sizes) if s >= min_cluster])

    kept = np.isin(labeled, list(keep_labels))
    labeled_final, nlab_final = ndimage.label(kept.astype(np.int8), structure=np.ones((3,3,3), dtype=np.int8))
    return kept, labeled_final, nlab_final, p_star, sizes

def extract_cluster_table(z_img, labeled_img):
    """生成簇信息：id, size_vox, size_mm3, peak_z, peak_ijk, peak_xyz(mm)"""
    z = z_img.get_fdata(dtype=np.float32)
    lab = labeled_img.get_fdata().astype(int)
    affine = z_img.affine
    voxvol = voxel_volume_mm3(affine)

    labels = sorted(list(set(np.unique(lab)) - {0}))
    rows = []
    for L in labels:
        mask = (lab == L)
        size_vox = int(mask.sum())
        size_mm3 = float(size_vox * voxvol)
        # 峰值
        z_mask = np.where(mask, z, -np.inf)
        peak_idx = np.unravel_index(np.argmax(z_mask), z.shape)
        peak_z = float(z[peak_idx])
        ijk = np.array(peak_idx + (1,))[:3]  # just to shape
        ijk = np.array(peak_idx)
        xyz = nib.affines.apply_affine(affine, ijk[::-1][::-1])  # 保持 ijk → xyz
        # 直接用标准做法：
        ijk_arr = np.array(peak_idx)
        xyz_mm = nib.affines.apply_affine(affine, ijk_arr)
        rows.append({
            "cluster_id": int(L),
            "size_vox": size_vox,
            "size_mm3": round(size_mm3, 2),
            "peak_z": round(peak_z, 4),
            "peak_i": int(ijk_arr[0]),
            "peak_j": int(ijk_arr[1]),
            "peak_k": int(ijk_arr[2]),
            "peak_x": round(float(xyz_mm[0]), 2),
            "peak_y": round(float(xyz_mm[1]), 2),
            "peak_zmm": round(float(xyz_mm[2]), 2),
        })
    return pd.DataFrame(rows)

# ========== 主流程 ==========
def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    subjects = sorted([p.name for p in FIRSTLEVEL.iterdir() if p.is_dir() and p.name.startswith("sub-")])

    summary_rows = []

    for sub in subjects:
        for ses in SESSIONS:
            run_dirs = find_runs_for_subject_session(sub, ses)
            if not run_dirs:
                continue

            # 收集 run 级文件
            run_items = []
            for rdir in run_dirs:
                # 约定命名：*_con-{CON_NAME}_beta.nii.gz / _z.nii.gz / _var.nii.gz(若有)
                beta = first_existing(sorted(rdir.glob(f"*con-{CON_NAME}_beta.nii*")))
                zmap = first_existing(sorted(rdir.glob(f"*con-{CON_NAME}_z.nii*")))
                var  = first_existing(sorted(rdir.glob(f"*con-{CON_NAME}_var.nii*")))  # 若第一步没保存 var，这里大多 None
                if beta is None or zmap is None:
                    continue
                run_items.append({"beta": beta, "z": zmap, "var": var})

            if not run_items:
                print(f"[SKIP] {sub} {ses}: 未找到 {CON_NAME} 的 run 级文件")
                continue

            # 逆方差加权合并
            beta_img, var_img, z_img, valid_mask = fixedfx_combine_beta_var(run_items)

            # 构造“功能灰质 mask”
            final_mask = valid_mask.copy()
            if USE_GM_ONLY:
                bm, gm = load_subject_masks_in_ref_grid(sub, z_img)
                if bm is not None:
                    final_mask = final_mask & bm
                if gm is not None:
                    final_mask = final_mask & (gm >= GM_PROB_THR)
            if final_mask.sum() == 0:
                print(f"[WARN] {sub} {ses}: 功能灰质限制后无有效体素")
                continue

            # 阈值 + 聚类（单侧 FDR）
            kept, labeled, nlab, p_star, _sizes = threshold_and_cluster(
                z_img, base_mask=final_mask, side=SIDE, q=Q_FDR, z_vox=Z_VOX, min_cluster=MIN_CLUSTER
            )

            # 输出路径
            out_base = OUTDIR / sub / ses
            out_base.mkdir(parents=True, exist_ok=True)
            prefix = f"{sub}_{ses}_con-{CON_NAME}_fixedfx"

            beta_file   = out_base / f"{prefix}_beta.nii.gz"
            var_file    = out_base / f"{prefix}_var.nii.gz"
            z_file      = out_base / f"{prefix}_z.nii.gz"
            mask_file   = out_base / f"{prefix}_mask_functionalGM_thr.nii.gz"
            fdr_file    = out_base / f"{prefix}_thr_fdrq{Q_FDR}_zvox{Z_VOX}_k{MIN_CLUSTER}.nii.gz"
            label_file  = out_base / f"{prefix}_thr_fdrq{Q_FDR}_zvox{Z_VOX}_k{MIN_CLUSTER}_labels.nii.gz"
            table_file  = out_base / f"{prefix}_clusters.tsv"
            meta_file   = out_base / f"{prefix}_meta.json"

            nib.save(beta_img, str(beta_file))
            nib.save(var_img,  str(var_file))
            nib.save(z_img,    str(z_file))
            nib.save(nib.Nifti1Image(final_mask.astype(np.uint8), z_img.affine, header=z_img.header), str(mask_file))
            nib.save(nib.Nifti1Image(kept.astype(np.uint8), z_img.affine, header=z_img.header), str(fdr_file))
            nib.save(nib.Nifti1Image(labeled.astype(np.int32), z_img.affine, header=z_img.header), str(label_file))

            # 簇表
            if nlab > 0 and kept.any():
                table = extract_cluster_table(z_img, nib.Nifti1Image(labeled, z_img.affine, header=z_img.header))
                table.to_csv(table_file, sep="\t", index=False)
            else:
                pd.DataFrame(columns=["cluster_id","size_vox","size_mm3","peak_z","peak_i","peak_j","peak_k","peak_x","peak_y","peak_zmm"]).to_csv(table_file, sep="\t", index=False)

            # 元数据
            meta = {
                "subject": sub, "session": ses, "contrast": CON_NAME,
                "combine": "fixed-effects variance-weighted",
                "use_gm_only": bool(USE_GM_ONLY), "gm_prob_thr": float(GM_PROB_THR),
                "fdr_q": float(Q_FDR), "side": SIDE, "z_vox": float(Z_VOX), "min_cluster": int(MIN_CLUSTER),
                "n_runs": len(run_items), "final_mask_voxels": int(final_mask.sum()),
                "p_star_fdr": float(p_star),
                "inputs_example": {k: str(run_items[0][k]) for k in run_items[0].keys()},
                "notes": "如未提供每个 run 的方差图，本脚本用 var≈(beta/z)^2 反推；z≈0 的体素自动极小权重"
            }
            meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")

            # 仅在 ROI_SESSIONS 上生成 “用于 PPI 的 ROI” 的软链接（可选）
            if ses in ROI_SESSIONS:
                roi_link_dir = OUTDIR / sub / "roi_for_ppi"
                roi_link_dir.mkdir(parents=True, exist_ok=True)
                for src, name in [(fdr_file, "roi_mask.nii.gz"), (label_file, "roi_labels.nii.gz"), (table_file, "roi_clusters.tsv")]:
                    dst = roi_link_dir / f"{ses}_{name}"
                    try:
                        if dst.exists() or dst.is_symlink():
                            dst.unlink()
                        dst.symlink_to(src)
                    except Exception:
                        # 在不支持 symlink 的系统上就复制一份
                        import shutil
                        shutil.copyfile(src, dst)

            summary_rows.append([sub, ses, str(beta_file), str(var_file), str(z_file), str(mask_file), str(fdr_file), str(label_file), str(table_file)])
            print(f"[OK] {sub} {ses}: fixedfx 合并与 ROI 阈值完成。")

    if summary_rows:
        df = pd.DataFrame(summary_rows, columns=["subject","session","beta","var","z","mask_functionalGM","roi_mask","roi_labels","clusters_tsv"])
        df.to_csv(OUTDIR / "roi_summary.tsv", sep="\t", index=False)
        print(f"\n[SUMMARY] Saved: {OUTDIR / 'roi_summary.tsv'}")
    else:
        print("[WARNING] 未生成任何 ROI，请检查输入与阈值参数。")

if __name__ == "__main__":
    main()

