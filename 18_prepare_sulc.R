# 18_prepare_sulc.R -- 脑沟深度 (sulc) 数据准备: 原始均值 -> 对照组 z 标准化
# 输入:  F:/DJ/sulc_extract/sulc_regions_all.csv  (17_sulc_extract.py 产出, 542 被试 x 68 区)
# 输出:  ideas_prep_sulc/ (组图 + 个体矩阵, 格式与 thick/area/vol 一致)
suppressPackageStartupMessages({ library(BrainEnrich); library(data.table) })
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"

gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
gene_regions <- rownames(gene_data)
strip <- function(x) sub("^(lh|rh|L|R)_", "", x, ignore.case = TRUE)
gene_core <- strip(gene_regions)

d <- fread("F:/DJ/sulc_extract/sulc_regions_all.csv", colClasses = list(character = "ID"))
d$ID <- as.character(d$ID)
regions68 <- grep("^(lh|rh)_", names(d), value = TRUE)
stopifnot(length(regions68) == 68)
cat(sprintf("sulc 表: %d 被试 x %d 区\n", nrow(d), length(regions68)))

idn <- suppressWarnings(as.numeric(d$ID))
is_ctl <- !is.na(idn) & idn >= 4000   # 对照 = 4001-4100; 患者 1-442 (含 400, 不能用前缀匹配)
cat(sprintf("对照 %d | 患者 %d\n", sum(is_ctl), sum(!is_ctl)))
stopifnot(sum(is_ctl) == 100)

mat <- as.matrix(d[, ..regions68])
rownames(mat) <- d$ID
cat(sprintf("原始值 NA 单元: %d\n", sum(is.na(mat))))

# 对照组 z 标准化 (逐区)
cm <- colMeans(mat[is_ctl, , drop = TRUE], na.rm = TRUE)
cs <- apply(mat[is_ctl, , drop = FALSE], 2, sd, na.rm = TRUE)
zmat <- sweep(sweep(mat, 2, cm, "-"), 2, cs, "/")
cat(sprintf("对照 z 均值范围 [%.3f, %.3f] | sd 范围 [%.3f, %.3f]\n",
            min(cm - cm), max(range(colMeans(zmat[is_ctl, ], na.rm = TRUE))), min(cs), max(cs)))
# 注: sulc 原值在脑沟处为负; z<0 = 该区脑沟比对照更深(增深), z>0 = 变浅

zp <- zmat[!is_ctl, , drop = FALSE]
lh <- grep("^lh_", regions68, value = TRUE)
idx <- match(strip(lh), gene_core)
if (any(is.na(idx))) { print(lh[is.na(idx)]); stop("sulc 区名对齐失败") }

out <- file.path(ws, "ideas_prep_sulc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# 组水平效应图 (患者均值)
effect_all <- colMeans(zp, na.rm = TRUE)
effect_lh <- effect_all[lh]; names(effect_lh) <- gene_regions[idx]
gdf <- data.frame(effect = effect_lh, row.names = gene_regions[idx])
gdf <- gdf[order(match(rownames(gdf), gene_regions)), , drop = FALSE]
fwrite(data.frame(region = rownames(gdf), effect = gdf$effect), file.path(out, "brain_data_group_lh.csv"))
saveRDS(gdf, file.path(out, "brain_data_group_lh.rds"))

# 个体矩阵 34 x N
zl <- zp[, lh, drop = FALSE]
colnames(zl) <- gene_regions[idx]
indiv <- t(zl)
saveRDS(indiv, file.path(out, "brain_data_indiv_lh.rds"))
cat(sprintf("组图 34 区 (最负3: %s) | 个体矩阵 %dx%d | 个体 NA: %d\n",
            paste(names(head(sort(gdf$effect), 3)), collapse = ","),
            nrow(indiv), ncol(indiv), sum(is.na(zl))))
cat("\nSULC PREP DONE\n")
