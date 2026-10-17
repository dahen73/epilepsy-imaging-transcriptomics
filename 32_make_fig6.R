## 32_make_fig6.R — 节点2nd 收尾: lh vs rh 定量对比 + fig6
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ws  <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
rdir <- file.path(ws, "full_run/rh_individual"); dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
prog <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

d <- fread(file.path(ws, "full_run/rh/rh_vs_lh.csv"))
setnames(d, "gene_set", "gs")

## 簇标签 (与论文三大簇一致)
d[, cluster := fcase(
  mod == "area" & pred == "FCDII_vs_other", "Cluster 1: area x FCDII",
  mod == "vol"  & pred == "Number_ASMs",    "Cluster 2: vol x ASM",
  mod == "vol"  & pred == "age_onset",      "Cluster 3: vol x age_onset",
  default = "Other")]

prog("== beta 跨半球相关 (Spearman, 按簇) ==")
cor_tab <- d[, .(n = .N,
                 r_spearman = round(cor(beta_lh, beta_rh, method = "spearman"), 3),
                 r_pearson  = round(cor(beta_lh, beta_rh), 3)),
             by = cluster][order(-n)]
print(cor_tab)
fwrite(cor_tab, file.path(rdir, "beta_lh_rh_correlation.csv"))

## 全 86 集整体相关
r_all <- cor(d$beta_lh, d$beta_rh, method = "spearman")
prog(sprintf("全 86 集 beta 跨半球 Spearman r = %.3f", r_all))

## lh FDR 显著集在 rh 的统计
sig <- d[fdr_lh < 0.05]
prog(sprintf("lh FDR 显著 %d 集: rh 方向一致 %d (%.0f%%), rh 中位 fdr = %.3f, rh beta 中位缩至 lh 的 %.2f",
             nrow(sig), sum(sig$dir_agree), 100 * mean(sig$dir_agree),
             median(sig$fdr_rh), median(abs(sig$beta_rh) / abs(sig$beta_lh))))

## ---- fig6 ----
lab_beta_lh <- "beta, left hemisphere (34 regions)"
lab_beta_rh <- "beta, right hemisphere (32 regions)"
cols <- c("#C0392B", "#2471A3", "#1E8449")

ann <- d[cluster != "Other", .(lab = sprintf("%s\nSpearman r = %.2f | dir kept %d/%d",
        gsub("Cluster [0-9]+: ", "", unique(cluster))[1],
        cor(beta_lh, beta_rh, method = "spearman"),
        sum(sign(beta_rh) == sign(beta_lh) & fdr_lh < 0.05), sum(fdr_lh < 0.05))),
        by = cluster]

p1 <- ggplot(d[cluster != "Other"], aes(beta_lh, beta_rh, color = cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey70") +
  geom_point(alpha = 0.75, size = 1.9) +
  geom_point(data = d[cluster != "Other" & fdr_lh < 0.05],
             shape = 21, fill = NA, color = "black", size = 2.6, stroke = 0.55) +
  facet_wrap(~cluster, scales = "free") +
  scale_x_continuous(labels = scales::label_scientific(digits = 1),
                     breaks = scales::breaks_pretty(3)) +
  scale_color_manual(values = cols) +
  labs(x = expression(beta*","~left~hemisphere~"(34 regions)"),
       y = expression(beta*","~right~hemisphere~"(32 regions)"),
       color = NULL) +
  geom_text(data = ann, aes(label = lab), x = Inf, y = -Inf, hjust = 1.06, vjust = -0.28,
            size = 2.6, color = "black", inherit.aes = FALSE) +
  theme_bw(base_size = 10) + theme(legend.position = "none",
    plot.subtitle = element_text(size = 8))

## 名义 p_emp 对比 (-log10 尺度, 不截断)
p2 <- ggplot(d[cluster != "Other"],
             aes(-log10(p_emp_lh), -log10(p_emp_rh), color = cluster)) +
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "grey40") +
  geom_vline(xintercept = -log10(0.05), linetype = "solid", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey70") +
  geom_point(alpha = 0.65, size = 1.8) +
  facet_wrap(~cluster) +
  scale_color_manual(values = cols) +
  labs(x = expression(-log[10]*" "*italic(p)[emp]*", left"),
       y = expression(-log[10]*" "*italic(p)[emp]*", right"),
       subtitle = "Nominal spin p; lines: p = 0.05", color = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

fig <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1))
ggsave(file.path(ws, "figures/fig6_hemisphere_comparison.png"), fig,
       width = 10, height = 4.4, dpi = 300, bg = "white")
prog("fig6 saved")
