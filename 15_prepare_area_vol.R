# 15_prepare_area_vol.R -- 为 area / vol (皮层68区) 构建 BrainEnrich 输入
# 输出: ideas_prep_area/ 与 ideas_prep_vol/ (组图 + 个体矩阵)
suppressPackageStartupMessages({ library(BrainEnrich); library(data.table) })
dj <- "F:/DJ"; ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"

gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
gene_regions <- rownames(gene_data)
strip <- function(x) {
  x <- sub("^(lh|rh|L|R)_", "", x, ignore.case = TRUE)
  x <- sub("_(thickness|thickavg|area|volume)$", "", x)
  x
}
gene_core <- strip(gene_regions)

for (mod in c("area", "vol")) {
  out <- file.path(ws, paste0("ideas_prep_", mod))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  zp <- fread(file.path(dj, sprintf("zscores/%s_zscores_patients.csv", mod)))
  zc <- fread(file.path(dj, sprintf("zscores/%s_zscores_controls.csv", mod)))
  # 只保留 lh_/rh_ 皮层区 (vol 的 aseg 皮下列剔除)
  regions68 <- grep("^(lh|rh)_", names(zp), value = TRUE)
  cat(sprintf("[%s] 患者 %dx%d (皮层区 %d) 对照 %d\n", mod, nrow(zp), ncol(zp), length(regions68), nrow(zc)))
  stopifnot(length(regions68) == 68)
  zmat_c <- as.matrix(zc[, ..regions68])
  cat(sprintf("    对照z均值范围 [%.3f, %.3f]\n", min(colMeans(zmat_c, na.rm = TRUE)), max(colMeans(zmat_c, na.rm = TRUE))))
  zmat_p <- as.matrix(zp[, ..regions68])
  effect_all <- colMeans(zmat_p, na.rm = TRUE)
  lh <- grep("^lh_", regions68, value = TRUE)
  idx <- match(strip(lh), gene_core)
  if (any(is.na(idx))) { print(lh[is.na(idx)]); stop(sprintf("%s 区名对齐失败", mod)) }
  effect_lh <- effect_all[lh]; names(effect_lh) <- gene_regions[idx]
  gdf <- data.frame(effect = effect_lh)
  gdf <- gdf[order(match(rownames(gdf), gene_regions)), , drop = FALSE]
  fwrite(data.frame(region = rownames(gdf), effect = gdf$effect), file.path(out, "brain_data_group_lh.csv"))
  saveRDS(gdf, file.path(out, "brain_data_group_lh.rds"))
  zl <- as.matrix(zp[, ..lh])
  rownames(zl) <- as.character(zp$ID); colnames(zl) <- gene_regions[idx]
  indiv <- t(zl)
  saveRDS(indiv, file.path(out, "brain_data_indiv_lh.rds"))
  cat(sprintf("[%s] 组图 34 区 (最负5: %s) | 个体矩阵 %dx%d | NA: %d\n", mod,
              paste(names(head(sort(gdf$effect), 3)), collapse = ","),
              nrow(indiv), ncol(indiv), sum(is.na(zl))))
}
cat("\nALL DONE\n")
