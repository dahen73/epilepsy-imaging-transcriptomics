# 20_cross_modal_summary.R -- 四模态靶向分析最终跨模态汇总
suppressPackageStartupMessages(library(data.table))
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
mods <- c("thickness", "area", "vol", "sulc")
dirmod <- c(thickness = "full_run/targeted", area = "full_run/targeted_area",
            vol = "full_run/targeted_vol", sulc = "full_run/targeted_sulc")
preds <- c("Number_ASMs", "age_onset", "HS_vs_other", "SE", "ILAE_Y1", "FCDII_vs_other")

con <- file(file.path(ws, "cross_modal_summary.txt"), "wt")

## 1) 全景矩阵: 每模态 x 预测变量 -> n_nominal / n_fdr
cat("== 1. 全景矩阵 (86 靶向集, spin 5000) ==\n\n", file = con)
hdr <- sprintf("%-12s %-16s %10s %10s", "modality", "predictor", "p<0.05", "FDR<0.05")
cat(hdr, "\n", file = con); cat(strrep("-", 52), "\n", file = con)
for (m in mods) for (p in preds) {
  f <- file.path(ws, dirmod[[m]], sprintf("individual_targeted_%s.csv", p))
  if (!file.exists(f)) { cat(sprintf("%-12s %-16s %10s\n", m, p, "missing"), file = con); next }
  d <- fread(f)
  cat(sprintf("%-12s %-16s %10d %10d\n", m, p, sum(d$p_emp < 0.05), sum(d$p_emp_fdr < 0.05)), file = con)
}

## 2) 所有 FDR 存活者明细
cat("\n== 2. FDR<0.05 存活基因集明细 ==\n", file = con)
for (m in mods) for (p in preds) {
  f <- file.path(ws, dirmod[[m]], sprintf("individual_targeted_%s.csv", p))
  if (!file.exists(f)) next
  d <- fread(f); d <- d[p_emp_fdr < 0.05][order(p_emp)]
  if (nrow(d) == 0) next
  cat(sprintf("\n--- %s x %s (%d FDR survivors) ---\n", m, p, nrow(d)), file = con)
  for (i in seq_len(nrow(d)))
    cat(sprintf("  %-46s p_emp=%.4f fdr=%.4f beta=%+.3f t=%+.1f\n",
                substr(d$description[i], 1, 46), d$p_emp[i], d$p_emp_fdr[i], d$beta[i], d$t_obs[i]), file = con)
}

## 3) 家族相干性: 名义显著集的家族标签
cat("\n== 3. GABA/癫痫家族跨模态方向一致性 (名义显著者中) ==\n", file = con)
fam <- function(desc) {
  desc <- tolower(desc)
  ifelse(grepl("gaba", desc), "GABA",
  ifelse(grepl("mtor|tor signals|torc", desc), "mTOR",
  ifelse(grepl("seizure|epilep", desc), "EPILEPSY_DGN", NA)))
}
for (p in preds) {
  cat(sprintf("\n-- %s --\n", p), file = con)
  for (m in mods) {
    f <- file.path(ws, dirmod[[m]], sprintf("individual_targeted_%s.csv", p))
    if (!file.exists(f)) next
    d <- fread(f); d[, fam := fam(description)]
    d <- d[!is.na(fam)]
    for (ff in unique(d$fam)) {
      sub <- d[fam == ff]; sig <- sub[p_emp < 0.05]
      if (nrow(sig) > 0) {
        cat(sprintf("  %-10s %-12s n_sig=%2d/%2d  beta: %+d 负 / %+d 正\n",
                    m, ff, nrow(sig), nrow(sub), sum(sig$beta < 0), sum(sig$beta > 0)), file = con)
      }
    }
  }
}

## 4) region 级 age_onset 四模态
cat("\n== 4. region 级 age_onset 梯度 (34 区) ==\n", file = con)
for (m in mods) {
  f <- file.path(ws, dirmod[[m]], "region_age_onset.csv")
  if (!file.exists(f)) next
  d <- fread(f); d <- d[order(p)]
  cat(sprintf("%-10s: %d/34 fdr<0.05 | top3: %s\n", m, sum(d$p_fdr < 0.05),
              paste(sprintf("%s(p=%.4g)", substr(d$region, 1, 22), head(d$p, 3)), collapse = ", ")), file = con)
}
close(con)
cat("SUMMARY WRITTEN\n")
