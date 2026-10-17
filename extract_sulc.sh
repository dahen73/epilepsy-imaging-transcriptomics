#!/bin/bash
# 从 542 个 FreeSurfer 包中仅提取 sulc + aparc.annot (每包 4 个小文件)
mkdir -p /f/DJ/sulc_extract
cd /f/DJ/freesurfer_data
ls *.tar.gz | sed 's/\.tar\.gz//' > /f/DJ/sulc_extract/_ids.txt
run_one() {
  id=$1
  d="/f/DJ/sulc_extract/$id"
  if [ -f "$d/label/lh.aparc.annot" ] && [ -f "$d/surf/lh.sulc" ]; then return 0; fi
  tar -xzf "/f/DJ/freesurfer_data/$id.tar.gz" -C /f/DJ/sulc_extract \
      --wildcards "$id/surf/?h.sulc" "$id/label/?h.aparc.annot" 2>/dev/null
}
export -f run_one
cat /f/DJ/sulc_extract/_ids.txt | xargs -P 4 -I {} bash -c 'run_one {}'
echo "EXTRACT_DONE $(ls /f/DJ/sulc_extract | grep -c '^[0-9]') subjects"
