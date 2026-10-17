# =====================================================================
# 29_prepare_rh.R -- 节点2nd: 右半球敏感性分析·数据准备 (2026-10-15)
#
# 注: get_geneExp(hem="R") 仅返回 32 区 (AHBA rh 有 2 区缺采样), rh 分析以
#     32 区进行, 与 lh 34 区对比时需注明。
# 输出:
#   ideas_prep_rh/brain_data_indiv_{thick,area,vol,sulc}_rh.rds  (32 区 x N 患者, 行名 R_xxx)
#   full_run/rh/perm_id_dk_rh_5000.rds  (rotate_parcellation 自建, 32x5000)
# =====================================================================
suppressPackageStartupMessages({ library(BrainEnrich); library(data.table) })
options(stringsAsFactors = FALSE)
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
out <- file.path(ws, "ideas_prep_rh"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(ws, "prep_rh_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logf, append = TRUE); flush.console() }
prog("PREP-RH start")

md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
strip <- function(x) sub("^(lh|rh|L|R)_", "", x, ignore.case = TRUE)
gene_data_rh <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "R"))
rh_regions <- rownames(gene_data_rh)                       # R_xxx (32)
nr <- length(rh_regions)
gene_core <- strip(rh_regions)
prog(sprintf("geneExp rh: %d regions", nr))

## ---------- 1. 三模态 z 表 rh 列 ----------
build_modal <- function(pat_csv, ctl_csv, suffix, tag) {
  zp <- fread(pat_csv, check.names = FALSE); zc <- fread(ctl_csv, check.names = FALSE)
  cn_p <- grep(paste0("^rh_.+", suffix, "$"), names(zp), value = TRUE)
  cn_c <- grep(paste0("^rh_.+", suffix, "$"), names(zc), value = TRUE)
  stopifnot(length(cn_p) == 34, length(cn_c) == 34)
  mp <- data.frame(ID = as.character(zp$ID), zp[, ..cn_p], check.names = FALSE)
  mc <- data.frame(ID = as.character(zc$ID), zc[, ..cn_c], check.names = FALSE)
  core_p <- sub(paste0("_", suffix, "$"), "", strip(cn_p))
  core_c <- sub(paste0("_", suffix, "$"), "", strip(cn_c))
  names(mp) <- c("ID", paste0("R_", core_p))
  names(mc) <- c("ID", paste0("R_", core_c))
  zall <- rbind(mp, mc)
  zall <- zall[order(as.numeric(zall$ID)), ]
  keep_id <- intersect(zall$ID, md$ID)
  mat <- t(as.matrix(zall[match(keep_id, zall$ID), -1, drop = FALSE]))   # 行=区(34) 列=受试者
  colnames(mat) <- keep_id
  ki <- match(rh_regions, rownames(mat))                                  # 34 -> 32 区
  stopifnot(!any(is.na(ki)))
  mat <- mat[ki, , drop = FALSE]
  rownames(mat) <- rh_regions
  saveRDS(mat, file.path(out, sprintf("brain_data_indiv_%s_rh.rds", tag)))
  prog(sprintf("%s: %d regions x %d subjects; NA=%d", tag, nrow(mat), ncol(mat), sum(is.na(mat))))
  mat
}
build_modal("F:/DJ/zscores/thick_zscores_patients.csv", "F:/DJ/zscores/thick_zscores_controls.csv", "thickness", "thick")
build_modal("F:/DJ/zscores/area_zscores_patients.csv",  "F:/DJ/zscores/area_zscores_controls.csv",  "area",      "area")
build_modal("F:/DJ/zscores/vol_zscores_patients.csv",   "F:/DJ/zscores/vol_zscores_controls.csv",   "volume",    "vol")

## ---------- 2. sulc rh (17 号 CSV 的 rh 列, 同 18 号流程) ----------
d <- fread("F:/DJ/sulc_extract/sulc_regions_all.csv", colClasses = list(character = "ID"))
d$ID <- as.character(d$ID)
regions_rh <- grep("^rh_", names(d), value = TRUE)
stopifnot(length(regions_rh) == 34)
idn <- suppressWarnings(as.numeric(d$ID))
is_ctl <- !is.na(idn) & idn >= 4000
mat0 <- as.matrix(d[, ..regions_rh]); rownames(mat0) <- d$ID
cm <- colMeans(mat0[is_ctl, , drop = FALSE], na.rm = TRUE)
cs <- apply(mat0[is_ctl, , drop = FALSE], 2, sd, na.rm = TRUE)
zmat <- sweep(sweep(mat0, 2, cm, "-"), 2, cs, "/")
zp <- zmat[!is_ctl, , drop = FALSE]
colnames(zp) <- paste0("R_", strip(regions_rh))
ki <- match(rh_regions, colnames(zp))
stopifnot(!any(is.na(ki)))
zp <- zp[, ki, drop = FALSE]
keep_id <- intersect(rownames(zp), md$ID)
msulc <- t(zp[match(keep_id, rownames(zp)), , drop = FALSE])
msulc <- msulc[, order(as.numeric(colnames(msulc))), drop = FALSE]
rownames(msulc) <- rh_regions
saveRDS(msulc, file.path(out, "brain_data_indiv_sulc_rh.rds"))
prog(sprintf("sulc: %d regions x %d patients; NA=%d", nrow(msulc), ncol(msulc), sum(is.na(msulc))))
if (any(is.na(msulc))) print(which(is.na(msulc), arr.ind = TRUE))

## ---------- 3. rh spin 置换 (镜像坐标 + rotate_parcellation) ----------
data(coord_dk_lh)
crd <- coord_dk_lh
crd$centroid1 <- -crd$centroid1          # x 镜像 -> rh 质心
ci <- match(gene_core, rownames(crd))
if (any(is.na(ci))) { cat("missing coords:", rownames(crd)[is.na(match(gene_core, rownames(crd)))], "\n"); stop("坐标区名对齐失败") }
crd <- crd[ci, ]
t0 <- Sys.time()
pi_rh <- rotate_parcellation(coord.r = as.matrix(crd), nrot = 5000, method = "hungarian", seed = 20261015)
prog(sprintf("rotate_parcellation returned dim %s in %.1f min",
             paste(dim(pi_rh), collapse = "x"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
dir.create(file.path(ws, "full_run/rh"), showWarnings = FALSE, recursive = TRUE)
if (ncol(pi_rh) >= 5000) pi_rh <- pi_rh[, 1:5000]
stopifnot(all(dim(pi_rh) == c(nr, 5000)))
saveRDS(pi_rh, file.path(ws, "full_run/rh/perm_id_dk_rh_5000.rds"))
prog("PREP-RH DONE")
