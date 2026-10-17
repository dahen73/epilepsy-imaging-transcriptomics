d <- read.csv("full_run/individual_param/DGN_param_all.csv")
cat("cols:", paste(names(d), collapse = ","), "\n")
con <- file("b1_top_hits.txt", "wt")
for (p in c("age_onset", "SE", "HS_vs_other")) {
  s <- d[d$pred == p & !is.na(d$p_param), ]
  s <- s[order(s$p_param), ]
  writeLines(c(paste0("== ", p, " =="), sprintf("%s | p=%.4g beta=%.3f | %s",
    head(s$gene_set, 5), head(s$p_param, 5), head(s$beta, 5), head(s$description, 5))), con)
}
close(con)
cat("DONE\n")
