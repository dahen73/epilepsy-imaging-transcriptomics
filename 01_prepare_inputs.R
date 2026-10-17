# 01_prepare_inputs.R
# 从 F:/DJ 的 IDEAS 公开表构建 BrainEnrich 输入（组水平效应量图 + 个体水平矩阵 + 表型表）
# 运行前提: BrainEnrich 已安装
suppressPackageStartupMessages({
  library(BrainEnrich)
  library(data.table)
})

dj  <- "F:/DJ"
out <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/ideas_prep"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---------- 1. BrainEnrich 基因表达数据（Desikan, 左半球） ----------
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
gene_regions <- rownames(gene_data)
cat("[1] gene_data 区域数:", length(gene_regions), " 基因数:", ncol(gene_data), "\n")

## 区名归一化: 去掉 lh_/rh_ 前缀与 _thickness/_thickavg 后缀, 便于跨表对齐
strip <- function(x) {
  x <- sub("^(lh|rh|L|R)_", "", x, ignore.case = TRUE)
  x <- sub("_(thickness|thickavg)$", "", x)
  x
}
gene_core <- strip(gene_regions)

## ---------- 2. 读取 IDEAS z-score 表 ----------
zp <- fread(file.path(dj, "zscores/thick_zscores_patients.csv"))
zc <- fread(file.path(dj, "zscores/thick_zscores_controls.csv"))
cat("[2] 患者 z 表:", nrow(zp), "x", ncol(zp),
    " 对照 z 表:", nrow(zc), "x", ncol(zc), "\n")

regions68 <- setdiff(names(zp), "ID")
stopifnot(length(regions68) == 68)

## 对照质检: z 应近似均值 0、方差 1
zmat_c <- as.matrix(zc[, ..regions68])
cat("    对照 z 均值范围: [", round(min(colMeans(zmat_c, na.rm=TRUE)), 3),
    ",", round(max(colMeans(zmat_c, na.rm=TRUE)), 3), "]\n", sep = "")

## ---------- 3. 组水平效应量图 (每区患者 z 均值 = 标准化萎缩量) ----------
zmat_p <- as.matrix(zp[, ..regions68])
effect_all <- colMeans(zmat_p, na.rm = TRUE)

lh_regions <- grep("^lh_", regions68, value = TRUE)   # 左半球 34 区
idx <- match(strip(lh_regions), gene_core)
if (any(is.na(idx))) {
  cat("!! 以下区域未能对齐到 gene_data:\n")
  print(lh_regions[is.na(idx)])
  stop("区名对齐失败, 请检查 gene_data 命名格式")
}
effect_lh <- effect_all[lh_regions]
names(effect_lh) <- gene_regions[idx]

brain_data_group <- data.frame(effect = effect_lh)
brain_data_group <- brain_data_group[
  order(match(rownames(brain_data_group), gene_regions)), , drop = FALSE]
fwrite(data.frame(region = rownames(brain_data_group),
                  effect = brain_data_group$effect),
       file.path(out, "brain_data_group_lh.csv"))
saveRDS(brain_data_group, file.path(out, "brain_data_group_lh.rds"))
cat("[3] 组水平效应量图 (34 区) 已保存; 萎缩最重 5 区:\n")
print(head(sort(brain_data_group$effect), 5))

## ---------- 4. 个体水平矩阵 (regions x subjects, 左半球) ----------
zl <- as.matrix(zp[, ..lh_regions])
if (any(is.na(zl))) cat("    注意: 患者 z 矩阵含 NA 单元格:", sum(is.na(zl)), "\n")
rownames(zl) <- as.character(zp$ID)
colnames(zl) <- gene_regions[idx]
brain_data_indiv <- t(zl)          # 34 区 x N 患者
saveRDS(brain_data_indiv, file.path(out, "brain_data_indiv_lh.rds"))
cat("[4] 个体水平矩阵:", nrow(brain_data_indiv), "区 x", ncol(brain_data_indiv), "患者\n")

## ---------- 5. 表型表 (协变量数值化) ----------
md <- fread(file.path(dj, "tables_metadata/Metadata_Release_Anon.csv"))
bin_mid <- function(b) {
  sapply(as.character(b), function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    if (grepl(" to ", x, fixed = TRUE)) {
      mean(as.numeric(strsplit(x, " to ", fixed = TRUE)[[1]]))
    } else {
      suppressWarnings(as.numeric(gsub("\\+", "", x)))
    }
  }, USE.NAMES = FALSE)
}
md$age_scan  <- bin_mid(md$Binned_Age_at_Scan)
md$age_onset <- bin_mid(md$Binned_Onset_Age)
md$sex       <- factor(md$Sex)
for (v in c("FUS", "fqFUS", "FBTCS", "fqFBTCS", "SE")) {
  if (v %in% names(md)) md[[paste0("n_", v)]] <-
    as.integer(as.logical(md[[v]]))
}
md$pathology <- factor(md$Pathology)
fwrite(md, file.path(out, "phenotypes.csv"))
saveRDS(md, file.path(out, "phenotypes.rds"))
cat("[5] 表型表:", nrow(md), "患者 x", ncol(md), "列; 病理分布:\n")
print(table(md$Pathology, useNA = "ifany"))

cat("\nALL DONE -> ", out, "\n")
