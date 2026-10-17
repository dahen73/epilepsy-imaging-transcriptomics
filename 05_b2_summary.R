d1 <- read.csv("full_run/individual/DGN_spin5000_all.csv")
d2 <- read.csv("full_run/individual/SynGO_spin5000_all.csv")
a <- rbind(cbind(lib = "DGN", d1), cbind(lib = "SynGO", d2))
sig <- a[!is.na(a$p_emp) & a$p_emp < 0.05, ]
sig <- sig[order(sig$p_emp), ]
con <- file("b2_sig_summary.txt", "wt")
writeLines(c(sprintf("B2 total rows: %d, p_emp<0.05: %d, p_emp<0.01: %d, p_emp_fdr<0.05: %d",
                     nrow(a), sum(a$p_emp < 0.05, na.rm = TRUE),
                     sum(a$p_emp < 0.01, na.rm = TRUE), sum(a$p_emp_fdr < 0.05, na.rm = TRUE)),
             "top 20 by p_emp (lib | gene_set | pred | beta | t | p_emp | p_emp_fdr | desc)"), con)
if (nrow(sig) > 0) {
  for (i in 1:min(nrow(sig), 20)) writeLines(sprintf("%s | %s | %s | beta=%.3f t=%.2f p_emp=%.4f fdr=%.4f | %s",
    sig$lib[i], sig$gene_set[i], sig$pred[i], sig$beta[i], sig$t_obs[i], sig$p_emp[i], sig$p_emp_fdr[i], sig$description[i]), con)
} else {
  writeLines("NO gene sets at p_emp<0.05 in B2 (all negative after spin test)", con)
}
close(con)
