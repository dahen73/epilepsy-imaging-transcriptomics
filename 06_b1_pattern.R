## 检查 B1 参数化结果中的模式：每个预测变量的信号分布
files <- list.files("full_run/individual_param", pattern = "_param_all\\.csv$", full.names = TRUE)
con <- file("b1_pattern_summary.txt", "wt")
for (f in files) {
  d <- read.csv(f)
  lib <- sub("_param_all\\.csv$", "", basename(f))
  preds <- unique(d$pred)
  for (p in preds) {
    s <- d[d$pred == p & !is.na(d$p_param), ]
    n <- nrow(s)
    n05 <- sum(s$p_param < 0.05)
    # 双侧 t 分布下检查正负系数是否失衡
    pos <- sum(s$p_param < 0.05 & s$beta > 0); neg <- sum(s$p_param < 0.05 & s$beta < 0)
    top <- s[order(s$p_param), ][1:min(3, n), ]
    writeLines(sprintf("%s | %s | n_set=%d n_p05=%d (pos=%d neg=%d)", lib, p, n, n05, pos, neg), con)
    for (i in seq_len(nrow(top)))
      writeLines(sprintf("    %.4g | beta=%.3f t=%.2f | %s", top$p_param[i], top$beta[i], top$t_obs[i], top$gene_set[i]), con)
  }
}
close(con)
cat("DONE\n")
