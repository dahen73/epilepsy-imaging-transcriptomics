# 22_make_figures.R -- 论文级主图: 四模态跨模态信号汇总
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
figd <- file.path(ws, "figures"); dir.create(figd, showWarnings = FALSE)
mods <- c("thickness", "area", "vol", "sulc")
dirmod <- c(thickness = "full_run/targeted", area = "full_run/targeted_area",
            vol = "full_run/targeted_vol", sulc = "full_run/targeted_sulc")
preds <- c("Number_ASMs", "age_onset", "HS_vs_other", "SE", "ILAE_Y1", "FCDII_vs_other")
pred_lab <- c(Number_ASMs = "ASM burden", age_onset = "Age at onset",
              HS_vs_other = "HS vs other", SE = "Status epilepticus",
              ILAE_Y1 = "ILAE Yr1 outcome", FCDII_vs_other = "FCD II vs other")
mod_lab <- c(thickness = "Thickness", area = "Surface area", vol = "Volume", sulc = "Sulcal depth")
read_res <- function(m, p) {
  f <- file.path(ws, dirmod[[m]], sprintf("individual_targeted_%s.csv", p))
  if (file.exists(f)) fread(f) else NULL
}
theme_pub <- function(sz = 11) theme_bw(sz) + theme(panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey93", color = NA),
  strip.text = element_text(face = "bold"), legend.title = element_text(size = 9),
  legend.text = element_text(size = 8), plot.title = element_text(face = "bold", size = 12))

## ---- Fig 1: 跨模态全景热图 ----
grid <- rbindlist(lapply(mods, function(m) rbindlist(lapply(preds, function(p) {
  d <- read_res(m, p)
  data.table(modality = mod_lab[[m]], predictor = pred_lab[[p]],
             n_nom = if (is.null(d)) 0L else sum(d$p_emp < 0.05),
             n_fdr = if (is.null(d)) 0L else sum(d$p_emp_fdr < 0.05))
}))))
grid$predictor <- factor(grid$predictor, levels = pred_lab[preds])
grid$modality <- factor(grid$modality, levels = mod_lab[mods])
p1 <- ggplot(grid, aes(modality, predictor, fill = n_fdr)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = ifelse(n_fdr > 0, sprintf("%d / %d", n_fdr, n_nom),
                              ifelse(n_nom > 0, sprintf("0 / %d", n_nom), ""))),
            size = 3.1, fontface = "bold",
            color = ifelse(grid$n_fdr >= 5, "white", "grey20")) +
  scale_fill_gradient(low = "#f7f7f5", high = "#b2182b", name = "FDR<0.05\ncount", breaks = c(0, 10, 28)) +
  scale_y_discrete(limits = rev(levels(grid$predictor))) +
  labs(title = "A  Targeted imaging-transcriptomic associations across modalities",
       x = NULL, y = NULL, subtitle = "86 prior gene sets | spin 5,000 | cell: FDR / nominal") +
  theme_pub(11) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.key.height = unit(1.1, "cm"))
ggsave(file.path(figd, "fig1_overview_heatmap.png"), p1, width = 7.4, height = 4.3, dpi = 300, bg = "white")

## ---- Fig 2: 三大信号簇 top 命中 ----
famof <- function(x) ifelse(grepl("gaba", tolower(x)), "GABAergic", "Epilepsy (DGN)")
panel_bar <- function(m, p, ntop, ttl) {
  d <- read_res(m, p); d <- d[order(p_emp)][1:ntop]
  d[, fam := factor(famof(description), levels = c("GABAergic", "Epilepsy (DGN)"))]
  d[, lab := factor(sub("\\s*\\((GO_[A-Z]+|DGN|Reactome)\\)$", "", description),
                    levels = rev(sub("\\s*\\((GO_[A-Z]+|DGN|Reactome)\\)$", "", description)))]
  ggplot(d, aes(-log10(p_emp), lab, fill = fam)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("FDR=%.3g", p_emp_fdr)), hjust = -0.08, size = 2.5, color = "grey25") +
    scale_fill_manual(values = c(GABAergic = "#4575b4", `Epilepsy (DGN)` = "#d73027"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(title = ttl, x = expression(-log[10](italic(p)[emp])), y = NULL) +
    theme_pub(9.5) + theme(axis.text.y = element_text(size = 7.2))
}
p2a <- panel_bar("area", "FCDII_vs_other", 12, "B  Surface area x FCD II (18 FDR survivors)")
p2b <- panel_bar("vol",  "Number_ASMs",     12, "C  Volume x ASM burden (28 FDR survivors)")
p2c <- panel_bar("vol",  "age_onset",       10, "D  Volume x age at onset")
fig2 <- p2a / p2b / p2c + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(figd, "fig2_signal_clusters.png"), fig2, width = 8.6, height = 10.6, dpi = 300, bg = "white")

## ---- Fig 3: vol region 级早发梯度 ----
rv <- fread(file.path(ws, dirmod[["vol"]], "region_age_onset.csv"))
rv <- rv[order(p)]
rv[, lab := sub("^L_", "", region)]
rv[, lab := factor(lab, levels = rev(lab))]
rv[, sig := fcase(p_fdr < 0.05, "FDR < 0.05", p < 0.05, "nominal p < 0.05", default = "n.s.")]
p3 <- ggplot(rv, aes(beta, lab, color = sig)) +
  geom_vline(xintercept = 0, color = "grey55", linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_text(data = rv[sig == "FDR < 0.05"], aes(label = sprintf("p=%.2g", p)),
            hjust = -0.28, size = 2.6, color = "grey20", show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = "#b2182b", "nominal p < 0.05" = "#ef8a62", "n.s." = "grey70"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.18))) +
  labs(title = "E  Cortical volume deviation vs age at onset (34 regions, left hemisphere)",
       subtitle = "beta < 0: earlier onset -> greater volume reduction | corrected for age at scan + sex",
       x = expression(beta ~ (volume ~ z ~ per ~ year ~ of ~ onset)), y = NULL) +
  theme_pub(10) + theme(axis.text.y = element_text(size = 7.4))
ggsave(file.path(figd, "fig3_region_gradient.png"), p3, width = 6.8, height = 6.2, dpi = 300, bg = "white")

## ---- Fig 4: 家族方向一致性 (GABA / 癫痫基因 两家族 x 模态 x 预测变量) ----
rows <- list()
for (m in mods) for (p in preds) {
  d <- read_res(m, p); if (is.null(d)) next
  d[, fam := famof(description)]
  for (ff in c("GABAergic", "Epilepsy (DGN)")) {
    sub <- d[fam == ff]; sig <- sub[p_emp < 0.05]
    rows[[length(rows) + 1]] <- data.table(modality = mod_lab[[m]], predictor = pred_lab[[p]],
      family = ff, n_sig = nrow(sig), n_tot = nrow(sub),
      pos_share = if (nrow(sig) == 0) NA_real_ else mean(sig$beta > 0))
  }
}
fd <- rbindlist(rows)
fd[, lab := sprintf("%d/%d", n_sig, n_tot)]
fd$modality <- factor(fd$modality, levels = mod_lab[mods])
fd$predictor <- factor(fd$predictor, levels = pred_lab[preds])
p4 <- ggplot(fd[n_sig > 0], aes(modality, predictor, size = n_sig, color = pos_share)) +
  geom_point() +
  geom_text(aes(label = lab), size = 2.5, color = "black", vjust = -1.25, show.legend = FALSE) +
  facet_wrap(~family) +
  scale_color_gradient2(low = "#2166ac", mid = "grey88", high = "#b2182b",
                        midpoint = 0.5, name = "share\nbeta > 0") +
  scale_size_continuous(range = c(2.5, 7), name = "n sig.") +
  labs(title = "F  Signal family coherence across modalities",
       subtitle = "point size = nominal p<0.05 count; label = n_sig/n_family; color = direction",
       x = NULL, y = NULL) +
  theme_pub(10) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(figd, "fig4_family_coherence.png"), p4, width = 8.2, height = 5.6, dpi = 300, bg = "white")

cat("FIGURES DONE\n")
