#!/bin/bash
# 等待 tar 提取完成 -> python 计算 68 区均值 -> R z 标准化准备
ws="C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
PY="/c/Users/dahen/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
for i in $(seq 1 400); do
  grep -q "EXTRACT_DONE" "$ws/sulc_extract_log.txt" 2>/dev/null && break
  sleep 30
done
echo "[chain] extract finished at $(date +%H:%M:%S)"
"$PY" "$ws/17_sulc_extract.py"
echo "[chain] python done at $(date +%H:%M:%S)"
export R_LIBS_USER="C:/Users/dahen/Rlibs/4.3"
"/c/Program Files/R/R-4.3.2/bin/Rscript.exe" "$ws/18_prepare_sulc.R"
echo "[chain] SULC PREP ALL DONE at $(date +%H:%M:%S)"
