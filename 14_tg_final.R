## 汇总六个预测变量的靶向结果
files <- c(Number_ASMs = "Number_ASMs", age_onset = "age_onset", HS_vs_other = "HS_vs_other",
           SE = "SE", ILAE_Y1 = "ILAE_Y1", FCDII_vs_other = "FCDII_vs_other")
con <- file("tg_final_summary.txt", "wt")
for (pn in names(files)) {
  d <- read.csv(sprintf("full_run/targeted/individual_targeted_%s.csv", files[[pn]]))
  d <- d[order(d$p_emp), ]
  writeLines(sprintf("== %s == | p_emp<0.05: %d | fdr<0.05: %d", pn,
                     sum(d$p_emp < 0.05), sum(d$p_emp_fdr < 0.05)), con)
  top <- head(d, 5)
  for (i in seq_len(nrow(top)))
    writeLines(sprintf("  p_emp=%.4f fdr=%.3f beta=%+.3f | %s",
                       top$p_emp[i], top$p_emp_fdr[i], top$beta[i], top$description[i]), con)
}
close(con)
cat("DONE\n")
