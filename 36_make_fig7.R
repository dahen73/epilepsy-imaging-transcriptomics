## 36_make_fig7.R - 3rd : ASM  +  + 
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ws    <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
adir  <- file.path(ws, "full_run/asm_mechanism")
prog  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

c  <- readRDS(file.path(adir, "asm_results_cache.rds"))
cmp <- as.data.table(c$cmp)
dsum <- as.data.table(c$dsum)

cols <- c("GABAergic" = "#7D3C98", "Epilepsy-DGN" = "#C0392B",
          "mTOR" = "#1E8449", "Other-GO/Reactome" = "grey62")

## (a) beta base vs extended (severity-adjusted)
p1 <- ggplot(cmp, aes(beta * 1e3, t_ext, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(shape = p_emp_fdr < 0.05), size = 2, alpha = 0.85) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "FDR q<0.05 (base)", `FALSE` = "not significant (base)"),
                     name = NULL) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = expression(beta["ASM, base model"]~" (x"~10^-3*")"),
       y = expression(t["ASM, severity-adjusted"]~"(FL 5000)"),
       subtitle = "Severity-adjusted vs base betas; 28/28 sets keep p_emp < 0.05 (n = 401)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")

## (b)  (:  GABA +  DGN)
top_gaba <- unique(dsum$variable)[grepl("1902710", unique(dsum$variable))][1]
top_dgn  <- unique(dsum$variable)[grepl("^DGN", unique(dsum$variable))][1]
dd <- dsum[variable %in% c(top_gaba, top_dgn)][!is.na(n_grp)]
dd[variable == top_gaba, lab := "GABA receptor complex\n(GO:1902710)"]
dd[variable == top_dgn,  lab := sprintf("%s\n(DGN)", sub(" \\(DGN\\)$", "", cmp[gs == top_dgn]$description))]
dd[, lab := factor(lab)]; dd[, n_grp := factor(n_grp, c("1", "2", "3", ">=4"))]
p2 <- ggplot(dd, aes(n_grp, mean, group = lab, color = lab)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12) +
  geom_line(linewidth = 0.4, alpha = 0.6) + geom_point(size = 2) +
  facet_wrap(~lab) +
  scale_color_manual(values = c("#7D3C98", "#C0392B"), guide = "none") +
  labs(x = "Number of ASMs at surgery", y = "Transcriptomic score (mean +/- SE)",
       subtitle = "Dose-response, two strongest sets") +
  theme_bw(base_size = 10)

## (c)  (base beta )
fam_lv <- c("GABAergic", "Epilepsy-DGN", "mTOR", "Other-GO/Reactome")
cmp[, family := factor(family, fam_lv)]
p3 <- ggplot(cmp, aes(family, beta * 1e3, fill = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_boxplot(outlier.size = 1.2, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.5) +
  scale_fill_manual(values = cols, guide = "none") +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, color = "black", fatten = 1.2) +
  labs(x = NULL, y = expression(beta["ASM, base model"]~" (x"~10^-3*")"),
       subtitle = "Family specificity (GABA vs rest: Wilcoxon p = 0.029)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig <- cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1.25, 0.85, 1))
ggsave(file.path(ws, "figures/fig7_asm_mechanism.png"), fig,
       width = 12.5, height = 4.4, dpi = 300, bg = "white")
prog("fig7 saved")
