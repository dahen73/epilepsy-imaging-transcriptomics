# 21_region_summary.R -- 四模态 region 级 age_onset 梯度汇总
suppressPackageStartupMessages(library(data.table))
mods <- c("thickness", "area", "vol", "sulc")
dirmod <- c(thickness = "full_run/targeted", area = "full_run/targeted_area",
            vol = "full_run/targeted_vol", sulc = "full_run/targeted_sulc")
con <- file("region_summary.txt", "wt")
for (m in mods) {
  d <- fread(file.path(dirmod[[m]], "region_age_onset.csv"))
  d <- d[order(p)]
  writeLines(sprintf("%-10s: %d/34 fdr<0.05 | top3: %s", m, sum(d$p_fdr < 0.05),
    paste(sprintf("%s(p=%.4g)", substr(head(d$region, 3), 1, 22), head(d$p, 3)), collapse = ", ")), con)
  # vol 的 7 个 FDR 区明细
  if (m == "vol") {
    sig <- d[p_fdr < 0.05]
    for (i in seq_len(nrow(sig)))
      writeLines(sprintf("    FDR: %-24s beta=%+.4f p=%.5g fdr=%.4g", sig$region[i], sig$beta[i], sig$p[i], sig$p_fdr[i]), con)
  }
}
close(con)
cat("OK\n")
