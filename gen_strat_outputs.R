## =====================================================================
## gen_strat_outputs.R — stratified-analysis figure + tables for the
## Epilepsia-formatted manuscript (built on 41_stratified_pathology.R results)
## Inputs : D:/DJ/Epilepsia_work/04_full_run_results/...
## Outputs: C:/Epilepsia_Submission/figures_new, C:/Epilepsia_Submission/tables_new
## =====================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

res_dir <- "D:/DJ/Epilepsia_work/04_full_run_results"
tab_dir <- "D:/DJ/Epilepsia_work/06_tables"
out_fig <- "C:/Epilepsia_Submission/figures_new"
out_tab <- "C:/Epilepsia_Submission/tables_new"

strat <- fread(file.path(res_dir, "stratified_pathology/stratified_all_results.csv"))
t2    <- fread(file.path(tab_dir, "table2_clusters.csv"))

get_cell <- function(s, m, p) strat[stratum == s & modality == m & pred == p]

## ---------- panel a: FDR-sig counts, full cohort vs strata ----------
full_cells <- list(
  list(lab = "Volume x ASM burden", s = NULL,   m = "vol",  p = "Number_ASMs",
       full_file = "targeted_vol/individual_targeted_Number_ASMs.csv"),
  list(lab = "Volume x Age at onset", s = NULL, m = "vol",  p = "age_onset",
       full_file = "targeted_vol/individual_targeted_age_onset.csv"),
  list(lab = "Volume x Status epilepticus", s = NULL, m = "vol", p = "SE",
       full_file = "targeted_vol/individual_targeted_SE.csv"),
  list(lab = "Area x FCD II vs other", s = NULL, m = "area", p = "FCDII_vs_other",
       full_file = "targeted_area/individual_targeted_FCDII_vs_other.csv")
)
rows <- list()
for (cc in full_cells) {
  fd <- fread(file.path(res_dir, cc$full_file))
  n_full <- sum(fd$p_emp_fdr < 0.05)
  for (s in c("HS", "nonHS")) {
    dd <- get_cell(s, cc$m, cc$p)
    rows[[length(rows) + 1]] <- data.table(cell = cc$lab,
                                           cohort = if (s == "nonHS") "non-HS" else "HS",
                                           n_fdr = if (nrow(dd)) sum(dd$p_emp_fdr < 0.05) else 0L)
  }
  rows[[length(rows) + 1]] <- data.table(cell = cc$lab, cohort = "Full cohort", n_fdr = n_full)
}
pa_df <- rbindlist(rows)
pa_df$cohort <- factor(pa_df$cohort, levels = c("Full cohort", "HS", "non-HS"))
pa_df$cell   <- factor(pa_df$cell, levels = sapply(full_cells, `[[`, "lab"))

pal <- c("Full cohort" = "grey70", "HS" = "#C44E52", "non-HS" = "#4C72B0")
pa <- ggplot(pa_df, aes(cell, n_fdr, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = n_fdr), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(0, 36), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "FDR-significant gene sets (of 86)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1),
        legend.position = "top", panel.grid.major.x = element_blank())

## ---------- panel b: full-cohort vs stratum betas (internal replication) ----
c1 <- t2[Cluster == "1: area x FCD II", gs]
c2 <- t2[Cluster == "2: vol x Number_ASMs", gs]
nh_fcd <- get_cell("nonHS", "area", "FCDII_vs_other")
nh_asm <- get_cell("nonHS", "vol",  "Number_ASMs")
fu_fcd <- fread(file.path(res_dir, "targeted_area/individual_targeted_FCDII_vs_other.csv"))
fu_asm <- fread(file.path(res_dir, "targeted_vol/individual_targeted_Number_ASMs.csv"))

mk <- function(full_dt, str_dt, sets, cl) {
  merge(full_dt[gene_set %in% sets, .(gene_set, beta_full = beta, q_full = p_emp_fdr)],
        str_dt[gene_set %in% sets, .(gene_set, beta_str = beta, q_str = p_emp_fdr)],
        by = "gene_set")[, cluster := cl]
}
pb_df <- rbind(mk(fu_fcd, nh_fcd, c1, "Cluster 1: area x FCD II"),
               mk(fu_asm, nh_asm, c2, "Cluster 2: volume x ASM"))
r1 <- cor(pb_df[grepl("Cluster 1", cluster), .(beta_full, beta_str)])[1, 2]
r2 <- cor(pb_df[grepl("Cluster 2", cluster), .(beta_full, beta_str)])[1, 2]
cat(sprintf("beta concordance r: cluster1=%.4f cluster2=%.4f\n", r1, r2))
ann <- data.table(cluster = unique(pb_df$cluster),
                  lab = sprintf("r = %.3f", c(r1, r2)))
pb <- ggplot(pb_df, aes(beta_full, beta_str)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(shape = q_str < 0.05), size = 2, colour = "#4C72B0") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.1, vjust = 1.4, size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~cluster, scales = "free") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c("q >= 0.05 in stratum", "q < 0.05 in stratum"),
                     name = NULL) +
  labs(x = "Beta, full cohort", y = "Beta, non-HS stratum (N = 218)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.text = element_text(size = 9.5))

## ---------- panel c: HS-specific volume x status-epilepticus sets ----------
hs_se <- get_cell("HS", "vol", "SE")[p_emp_fdr < 0.05][order(p_emp_fdr)]
fam <- function(gs) fifelse(grepl("TORC1|mTOR", gs, ignore.case = TRUE), "mTOR/TORC1",
                     fifelse(grepl("GABA", gs, ignore.case = TRUE), "GABAergic", "Epilepsy (DGN)"))
hs_se[, family := fam(paste(gene_set, description))]
top <- hs_se[1:12]
top[, description := factor(description, levels = rev(description))]
fampal <- c("GABAergic" = "#4C72B0", "Epilepsy (DGN)" = "#C44E52", "mTOR/TORC1" = "#55A868")
pc <- ggplot(top, aes(-log10(p_emp), description, fill = family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("q = %.3f", p_emp_fdr)), hjust = -0.15, size = 3.1) +
  scale_fill_manual(values = fampal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = expression(-log[10](p[emp])), y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

## ---------- assemble ----------
fig <- (pa | pb) / pc +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a")
ggsave(file.path(out_fig, "fig3_stratified_pathology.png"), fig,
       width = 11.4, height = 8.6, dpi = 300)

## ---------- main Table 2 (compact cluster x stratum summary) ----------
nh_fcd_sig <- nh_fcd[p_emp_fdr < 0.05]; nh_asm_sig <- nh_asm[p_emp_fdr < 0.05]
hs_se_sig  <- hs_se; hs_ao_sig <- get_cell("HS", "vol", "age_onset")[p_emp_fdr < 0.05]
fu_ao      <- fread(file.path(res_dir, "targeted_vol/individual_targeted_age_onset.csv"))
c3n <- sum(fu_ao$p_emp_fdr < 0.05)

tab2 <- data.table(
  Cluster = c("1. Surface area x FCD II", "2. Volume x ASM burden",
              "3. Volume x age at onset", "4. Volume x status epilepticus (HS-stratified)"),
  `Full cohort (N = 425)` = c(sprintf("%d sets, q <= %.3f, all beta > 0", length(c1), 0.037),
                              sprintf("%d sets, q <= %.3f, all beta > 0", length(c2), 0.047),
                              sprintf("%d sets, q <= %.3f, all beta < 0", c3n, 0.043),
                              "0 sets (min q = 0.063)"),
  `HS stratum (N = 195-215)` = c("not tested (FCD II contrast undefined)", "0 sets (min q = 0.57)",
                                 sprintf("%d sets, q <= %.3f, all beta < 0", nrow(hs_ao_sig), max(hs_ao_sig$p_emp_fdr)),
                                 sprintf("%d sets, q <= %.3f, all beta > 0", nrow(hs_se_sig), max(hs_se_sig$p_emp_fdr))),
  `non-HS stratum (N = 206-226)` = c(sprintf("%d sets, q <= %.3f, all beta > 0", nrow(nh_fcd_sig), max(nh_fcd_sig$p_emp_fdr)),
                                     sprintf("%d sets, q <= %.3f, all beta > 0", nrow(nh_asm_sig), max(nh_asm_sig$p_emp_fdr)),
                                     "0 sets (min q = 0.16)", "0 sets (min q = 0.70)"),
  Interpretation = c("Replicates within non-HS (16/18 cluster sets recovered)",
                     "Replicates and strengthens within non-HS (all 28 sets recovered + 3); null in HS",
                     "Concentrated in HS (both cluster sets + 17 additional)",
                     "HS-specific; GABAergic + TORC1 + status-epilepticus gene sets")
)
fwrite(tab2, file.path(out_tab, "table2_clusters_stratified.csv"))

## ---------- Supplementary Table S5: all primary-stratum FDR hits ----------
s5 <- rbind(
  cbind(context = "HS x volume x status epilepticus",
        hs_se[, .(gene_set, description, N, beta, t_obs, p_emp, p_emp_fdr)]),
  cbind(context = "HS x volume x age at onset",
        hs_ao_sig[, .(gene_set, description, N, beta, t_obs, p_emp, p_emp_fdr)]),
  cbind(context = "non-HS x area x FCD II vs other",
        nh_fcd_sig[, .(gene_set, description, N, beta, t_obs, p_emp, p_emp_fdr)]),
  cbind(context = "non-HS x volume x ASM burden",
        nh_asm_sig[, .(gene_set, description, N, beta, t_obs, p_emp, p_emp_fdr)])
)
setorder(s5, context, p_emp_fdr)
fwrite(s5, file.path(out_tab, "suptable5_stratified_fdr_hits.csv"))

cat("DONE. fig + 2 tables written.\n")
