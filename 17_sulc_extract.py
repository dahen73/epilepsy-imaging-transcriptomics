# -*- coding: utf-8 -*-
# 17_sulc_extract.py -- 从提取的 surf/?h.sulc + label/?h.aparc.annot 计算 68 区平均脑沟深度
# 输出: /f/DJ/sulc_extract/sulc_regions_all.csv (每行一个被试)
import os, sys, csv
import numpy as np
from nibabel.freesurfer.io import read_morph_data, read_annot

base = r"F:/DJ/sulc_extract"
ids = [l.strip() for l in open(os.path.join(base, "_ids.txt")) if l.strip().isdigit()]

# 读取一次参考 annot 获取区名顺序 (取被试 1)
_, _, names_l = read_annot(os.path.join(base, "1", "label", "lh.aparc.annot"))
names_l = [n.decode() if isinstance(n, bytes) else n for n in names_l]

rows = []
bad = []
for sid in ids:
    try:
        row = {"ID": sid}
        for hemi in ("lh", "rh"):
            sulc = read_morph_data(os.path.join(base, sid, "surf", f"{hemi}.sulc"))
            lab, _, names = read_annot(os.path.join(base, sid, "label", f"{hemi}.aparc.annot"))
            names = [n.decode() if isinstance(n, bytes) else n for n in names]
            for k, nm in enumerate(names):
                if nm in ("unknown", "corpuscallosum", "Background+FreeSurfer_Defined_Mixed_Tissue"):
                    continue
                m = lab == k
                if m.sum() == 0:
                    row[f"{hemi}_{nm}"] = np.nan
                else:
                    row[f"{hemi}_{nm}"] = float(np.mean(sulc[m]))
        rows.append(row)
    except Exception as e:
        bad.append((sid, str(e)[:80]))
    if len(rows) % 50 == 0:
        print(f"processed {len(rows)}", flush=True)

# 统一列: 以被试1的 lh 区名为基准 (aparc 区名应一致)
cols = ["ID"] + [f"lh_{n}" for n in names_l if n not in ("unknown", "corpuscallosum")] \
       + [f"rh_{n}" for n in names_l if n not in ("unknown", "corpuscallosum")]
out = os.path.join(base, "sulc_regions_all.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})
print(f"DONE: {len(rows)} subjects written, {len(bad)} failed -> {out}")
for sid, msg in bad[:10]:
    print(f"  FAIL {sid}: {msg}")
