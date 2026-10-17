# -*- coding: utf-8 -*-
# 23_extract_aparc_stats.py -- 从自有 FreeSurfer subjects 目录提取 DK68 区值
# 用法: python 23_extract_aparc_stats.py <subjects_dir> <output_csv>
# 依赖: pip install pandas  (纯文本解析, 无需 nibabel)
# 读取每个被试的 stats/{lh,rh}.aparc.{thickness,area,volume}.stats
# 输出: 每行一被试, 列 = lh_bankssts_thickness ... rh_*.volume (同 IDEAS 格式)
import os, sys, glob
import pandas as pd

MEAS = {"thickness": "ThickAvg", "area": "SurfArea", "volume": "GrayVol"}
COLS = ["StructName", "NumVert", "SurfArea", "GrayVol", "ThickAvg", "MeanCurv", "GausCurv", "FoldInd", "CurvInd"]

def parse_stats(path, measure):
    want = MEAS[measure]
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= len(COLS) - 1:
                rows[parts[0]] = float(parts[COLS.index(want) - 1] if COLS.index(want) >= 1 else parts[1])
    # 注: StructName 占 parts[0], 数值从 parts[1] 起; COLS.index(want)-1 对应列位置
    return rows

def main(subj_dir, out_csv):
    ids = sorted([d for d in os.listdir(subj_dir)
                  if os.path.isdir(os.path.join(subj_dir, d, "stats"))])
    records, missing = [], []
    for sid in ids:
        rec = {"ID": sid}
        ok = True
        for hemi in ("lh", "rh"):
            for meas in MEAS:
                p = os.path.join(subj_dir, sid, "stats", f"{hemi}.aparc.{meas}.stats")
                if not os.path.exists(p):
                    ok = False; missing.append(f"{sid}/{hemi}.aparc.{meas}.stats"); continue
                for region, val in parse_stats(p, meas).items():
                    if region in ("unknown", "corpuscallosum"):
                        continue
                    rec[f"{hemi}_{region}_{meas}"] = val
        if ok:
            records.append(rec)
    df = pd.DataFrame(records)
    df.to_csv(out_csv, index=False)
    print(f"DONE: {len(records)}/{len(ids)} subjects -> {out_csv}")
    if missing:
        print(f"WARNING: {len(missing)} files missing, e.g. {missing[:3]}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python 23_extract_aparc_stats.py <subjects_dir> <output_csv>")
    main(sys.argv[1], sys.argv[2])
