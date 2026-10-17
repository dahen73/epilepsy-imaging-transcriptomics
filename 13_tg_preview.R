d <- read.csv("full_run/targeted/individual_targeted_Number_ASMs.csv")
d <- d[order(d$p_emp), ]
con <- file("tg_preview.txt", "wt")
writeLines(sprintf("Number_ASMs | n_set=%d | p_emp<0.05: %d | p_emp_fdr<0.05: %d | p_param_fdr<0.05: %d",
                   nrow(d), sum(d$p_emp < 0.05), sum(d$p_emp_fdr < 0.05), sum(d$p_param_fdr < 0.05)), con)
for (i in 1:min(10, nrow(d)))
  writeLines(sprintf("  p_emp=%.4f fdr=%.3f beta=%.3f | %s", d$p_emp[i], d$p_emp_fdr[i], d$beta[i], d$description[i]), con)
close(con)
cat("DONE\n")
