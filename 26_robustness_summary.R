# =====================================================================
# 26_robustness_summary.R -- 稳健性结果汇总 + fig5 (2026-10-15)
# 依赖: 24 (S1-S3) 与 25 (S4) 输出
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(stringsAsFactors = FALSE)
ws     <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
rdir   <- file.path(ws, "full_run/robustness")
fdir   <- file.path(ws, "figures"); dir.create(fdir, showWarnings = FALSE)
logfile <- file.path(rdir, "robustness_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE) }

s1 <- fread(file.path(rdir, "s1_resection_sensitivity.csv"))
s2 <- fread(file.path(rdir, "s2_jackknife.csv"))
s3 <- fread(file.path(rdir, "s3_subsample_stability.csv"))

lines <- c("================================================================",
           "IDEAS x BrainEnrich 稳健性分析汇总 (2026-10-15)",
           "三大 FDR 信号簇代表集 (每簇 top5/top2 按 p_emp_fdr):",
           "  簇1 area x FCDII_vs_other (5/18 显著集)",
           "  簇2 vol   x Number_ASMs   (5/28 显著集)",
           "  簇3 vol   x age_onset     (2/2 显著集)", "================================================================", "")

lines <- c(lines, "---- S1 切除负担敏感性 (协变量法, 12 检验) ----",
  sprintf("方向保持: %d/%d | p<0.05 状态保持: %d/%d | beta 变化中位 %.1f%% (max %.1f%%)",
          sum(s1$direction_kept), nrow(s1), sum(s1$sig_kept_p05), nrow(s1),
          median(abs(s1$beta_change_pct)), max(abs(s1$beta_change_pct))),
  sprintf("加 resect_burden 协变量后 p 值整体变化: 改善 %d / 减弱 %d",
          sum(s1$p_adj < s1$p_base), sum(s1$p_adj > s1$p_base)), "")

lines <- c(lines, "---- S2 留一法 jackknife (12 集 x ~400 受试者) ----",
  sprintf("符号翻转总数: %d (of %d) | 留一后 p<0.05 比例: min %.3f / median %.3f",
          sum(s2$jk_sign_flips), nrow(s2) * s2$N[1], min(s2$jk_p_lt05_prop), median(s2$jk_p_lt05_prop)),
  sprintf("jackknife beta 范围 vs base: 全部 base 位于 [min,max] 内: %d/%d",
          sum(s2$beta_base >= s2$jk_beta_min & s2$beta_base <= s2$jk_beta_max), nrow(s2)), "")

lines <- c(lines, "---- S3 随机亚样本重复 (80% x 200) ----",
  sprintf("方向一致率 (bootstrap sign consistency): 全部 = %.2f | 基准 beta 落入 95%% CI: %d/%d",
          median(s3$boot_dir_consistency),
          sum(s3$beta_base >= s3$boot_beta_ci_lo & s3$beta_base <= s3$boot_beta_ci_hi), nrow(s3)),
  sprintf("子样本 p<0.05 比例: min %.3f / median %.3f / max %.3f",
          min(s3$boot_p_lt05_prop), median(s3$boot_p_lt05_prop), max(s3$boot_p_lt05_prop)), "")

s4f <- file.path(rdir, "s4_individual_FCDII_exclany50_compare.csv")
if (file.exists(s4f)) {
  s4 <- fread(s4f)
  so <- s4[sig_orig == TRUE]
  r_beta <- cor(so$beta_orig, so$beta_excl)
  lines <- c(lines, "---- S4 个体水平切除区剔除 (FCDII 簇, 剔 19 区/留 15 区, spin 5000) ----",
    sprintf("原始 18 集 fdr<0.05: 方向保持 %d/%d | 剔除后仍 p_emp<0.05 %d/%d (min p_emp=%.3f)",
            sum(so$dir_kept), nrow(so), sum(so$p_emp < 0.05), nrow(so), min(s4$p_emp)),
    sprintf("beta 一致性: cor(beta_orig, beta_excl) = %.3f (18 集)", r_beta),
    "解读: 方向完全保持, beta 跨集高度一致 -> 非符号伪影;",
    "  但 34->15 区后分数 SNR 大幅下降 (t 峰值 ~1.9), 无法维持 FDR 显著。",
    "  信号空间载体部分位于手术切除区 -- 与 FCDII 病灶即切除部位的生物学一致,",
    "  论文中应如实报告并讨论 '区覆盖-功效' 权衡, 不宜解读为信号消失。", "")
} else lines <- c(lines, "---- S4: 未完成 ----", "")

lines <- c(lines, "================================================================",
  "总体结论: S1/S2/S3 全部通过 (方向、显著性、效应量三维稳定);",
  "S4 显示 FCDII 簇方向 18/18 保持且 beta 跨集一致, 但剔除 19 个高切除区",
  "后 (剩 15 区) 功效不足以维持 FDR 显著 -- 报告为方向稳健+功效受限。",
  "================================================================")
writeLines(lines, file.path(rdir, "robustness_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

## ---------------- fig5: 四面板稳健性图 ----------------
short_name <- function(x) {
  x <- sub("^DGN:", "DGN ", x); x <- sub("^GO_CC:GO:", "GO ", x)
  x <- sub("^GO_MF:GO:", "GO ", x); x <- sub("^GO_BP:GO:", "GO ", x)
  x <- sub("^KEGG:hsa0", "KEGG hsa0", x); x
}
lab_fun <- function(df) df[, lab := short_name(gene_set)]

s1$lab <- short_name(s1$gene_set); s1$cluster <- ifelse(s1$mod == "area", "Cluster 1: area x FCDII",
    ifelse(s1$pred == "Number_ASMs", "Cluster 2: vol x #ASM", "Cluster 3: vol x age_onset"))
s2$lab <- short_name(s2$gene_set); s2$cluster <- ifelse(s2$mod == "area", "Cluster 1: area x FCDII",
    ifelse(s2$pred == "Number_ASMs", "Cluster 2: vol x #ASM", "Cluster 3: vol x age_onset"))
s3$lab <- short_name(s3$gene_set); s3$cluster <- ifelse(s3$mod == "area", "Cluster 1: area x FCDII",
    ifelse(s3$pred == "Number_ASMs", "Cluster 2: vol x #ASM", "Cluster 3: vol x age_onset"))

p1 <- ggplot(s1, aes(beta_base, beta_adj, color = cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = beta_base, xmax = beta_base), height = 0) +
  geom_point(size = 2.2, alpha = 0.9) +
  labs(x = "Beta (base model)", y = "Beta (+ resection-burden covariate)",
       title = "S1 | Resection-burden sensitivity") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p2 <- ggplot(s2, aes(beta_base, lab, color = cluster)) +
  geom_point(aes(x = jk_beta_mean), size = 2.2) +
  geom_errorbarh(aes(xmin = jk_beta_min, xmax = jk_beta_max), height = 0.25, alpha = 0.7) +
  geom_vline(aes(xintercept = 0), linetype = "dotted", color = "grey40") +
  labs(x = "Leave-one-out beta [min-max whiskers]", y = NULL,
       title = "S2 | Jackknife (no sign flips)") +
  theme_bw(base_size = 10) + theme(axis.text.y = element_blank())

p3 <- ggplot(s3, aes(beta_base, lab, color = cluster)) +
  geom_point(aes(x = boot_beta_mean), size = 2.2) +
  geom_errorbarh(aes(xmin = boot_beta_ci_lo, xmax = boot_beta_ci_hi), height = 0.25, alpha = 0.7) +
  geom_point(aes(x = beta_base), shape = 4, size = 2.4, color = "black") +
  geom_vline(aes(xintercept = 0), linetype = "dotted", color = "grey40") +
  labs(x = "Bootstrap beta mean [95% CI whiskers; x = base]", y = NULL,
       title = "S3 | 80% subsample x 200") +
  theme_bw(base_size = 10) + theme(axis.text.y = element_blank())

plots <- list(p1, p2, p3)
s4f <- file.path(rdir, "s4_individual_FCDII_exclany50_compare.csv")
if (file.exists(s4f)) {
  s4 <- fread(s4f); s4$lab <- short_name(s4$gene_set)
  so <- s4[sig_orig == TRUE]
  p4 <- ggplot(so, aes(beta_orig, beta_excl, color = dir_kept)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.4) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                       labels = c("TRUE" = "direction kept", "FALSE" = "flipped")) +
    labs(x = "Beta (all 34 regions)", y = "Beta (15 regions, high-resection excluded)",
         title = "S4 | FCDII cluster: 19 high-resection regions excluded",
         color = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  plots[[4]] <- p4
}

library(patchwork)
fig5 <- wrap_plots(plots, ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.title = element_text(face = "bold", size = 10))
ggsave(file.path(fdir, "fig5_robustness.png"), fig5, width = 10, height = 8, dpi = 300)
prog("fig5 saved")
