## 34_asm_mechanism.R — 节点3rd: vol x ASM 簇解释强化 (严重度分解 + 剂量反应 + 家族特异性)
suppressPackageStartupMessages({library(data.table)})
ws    <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir <- file.path(ws, "full_run/asm_mechanism"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
prog  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
RNGkind("Mersenne-Twister"); set.seed(20261015)

## ---------- 数据 ----------
S    <- readRDS(file.path(ws, "full_run/robustness/Sobs_vol.rds"))     # 受试者 x 86 集
md   <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
idx  <- match(rownames(S), md$ID); stopifnot(!any(is.na(idx)))
mdc  <- md[idx, ]
gsn  <- colnames(S)

## 组装协变量
mdc$duration <- mdc$age_scan - mdc$age_onset
mdc$n_asm    <- as.numeric(mdc$Number_ASMs)
mdc$se_i     <- as.integer(mdc$SE)
mdc$path_g <- sub("^(HS|FCD|DNT|CAV|TUMOR).*", "\\1",
                  toupper(ifelse(is.na(mdc$Pathology), "OTHER", mdc$Pathology)))
tab <- table(mdc$path_g); rare <- names(tab)[tab < 20]
mdc$path_g[mdc$path_g %in% rare] <- "OTHER"   # 罕见病理并入 OTHER, 稳定模型
cov_base <- c("age_scan", "sex")
ext_ok   <- complete.cases(mdc[, c("duration", "fqFUS", "fqFBTCS", "se_i", "path_g", "n_asm")])
prog(sprintf("完全病例数: base=%d, extended=%d (总 %d)",
             sum(!is.na(mdc$n_asm)), sum(ext_ok), nrow(mdc)))
if (sum(ext_ok) < 350) {  # 频率缺失过多则降级
  mdc$fqFUS <- NULL; mdc$fqFBTCS <- NULL
  ext_ok <- complete.cases(mdc[, c("duration", "se_i", "path_g", "n_asm")])
  prog(sprintf("频率缺失多, extended 降级为 duration+SE+pathology: n=%d", sum(ext_ok)))
}
dat <- mdc[ext_ok, ]; Sv <- S[ext_ok, , drop = FALSE]

mk_cov <- function(df, which = c("base", "ext")) {
  if (which == "base") data.frame(age_scan = df$age_scan, sex = factor(df$sex))
  else {
    out <- data.frame(age_scan = df$age_scan, sex = factor(df$sex),
                      duration = df$duration, se_i = factor(df$se_i),
                      path_g = factor(df$path_g))
    if (!is.null(df$fqFUS))    out$fqFUS    <- df$fqFUS
    if (!is.null(df$fqFBTCS))  out$fqFBTCS  <- df$fqFBTCS
    out
  }
}
Xc_base <- as.matrix(model.matrix(~ ., mk_cov(dat, "base")))[, -1, drop = FALSE]
Xc_ext  <- as.matrix(model.matrix(~ ., mk_cov(dat, "ext")))[, -1, drop = FALSE]

## 基准结果 (16 号, cov=age_scan+sex) 用于对比
res_base <- fread(file.path(ws, "full_run/targeted_vol/individual_targeted_Number_ASMs.csv"))
setnames(res_base, "gene_set", "gs")

## Freedman-Lane: 多响应一次性算 t
## X_full = [Xc_ext, n_asm]; H0: n_asm 系数=0。简化残差 = Y - Xc_ext b0 - n_asm*0
fl_perm <- function(Y, Xc, pred, B = 5000) {
  n <- nrow(Y); k <- ncol(Y)
  X0 <- cbind(1, Xc)                       # 简化设计 (含截距)
  X1 <- cbind(1, Xc, pred)                 # 完整设计
  p1 <- ncol(X1)                           # pred 的位置
  XtXi <- chol2inv(qr.R(qr(X1)))
  b0   <- lm.fit(X0, Y)$coefficients       # 简化模型 (含截距)
  if (is.null(dim(b0))) b0 <- matrix(b0, ncol = 1)
  fitted0 <- X0 %*% b0
  Res <- Y - fitted0
  df  <- n - p1
  t_obs <- rep(NA_real_, k)
  b  <- as.matrix(XtXi %*% crossprod(X1, Y))
  r  <- Y - X1 %*% b
  s2 <- colSums(r^2) / df
  t_obs <- b[p1, ] / sqrt(s2 * XtXi[p1, p1])
  cnt <- integer(k)
  for (i in seq_len(B)) {
    yp  <- fitted0 + Res[sample.int(n), , drop = FALSE]         # 行重排 = 受试者标签置换
    bp  <- as.matrix(XtXi %*% crossprod(X1, yp))
    rp  <- yp - X1 %*% bp
    s2p <- colSums(rp^2) / df
    tp  <- bp[p1, ] / sqrt(s2p * XtXi[p1, p1])
    cnt <- cnt + (abs(tp) >= abs(t_obs))
  }
  data.frame(t_ext = t_obs, p_emp_ext = (1 + cnt) / (B + 1))
}

prog("Part A: Freedman-Lane 5000 x 86 集 (extended 模型) ...")
pred <- dat$n_asm
resA <- fl_perm(as.matrix(Sv), Xc_ext, pred, B = 5000)
resA$gs <- gsn
resA$p_emp_fdr_ext <- p.adjust(resA$p_emp_ext, "BH")

cmp <- merge(res_base[, .(gs, description, beta, p_emp, p_emp_fdr)],
             as.data.table(resA), by = "gs")
cmp[, dir_kept := sign(beta) > 0 & sign(t_ext) > 0 | sign(beta) < 0 & sign(t_ext) < 0]
cmp[, beta_retained := p_emp < 0.05 & p_emp_ext < 0.05]
fwrite(cmp, file.path(outdir, "asm_severity_decomposition.csv"))

sig <- cmp[p_emp_fdr < 0.05]
prog(sprintf("28 显著集: extended 后方向保持 %d/%d, p_emp_ext<0.05 %d/%d (最小 %.4f)",
             sum(sig$dir_kept), nrow(sig), sum(sig$p_emp_ext < 0.05), nrow(sig),
             min(sig$p_emp_ext)))
prog(sprintf("全 86 集: extended 后 p_emp_fdr<0.05 共 %d 集", sum(cmp$p_emp_fdr_ext < 0.05)))

## ---------- Part C: 家族特异性 (用 base 结果) ----------
fam <- function(gs, desc) {
  d <- tolower(paste(gs, desc))
  ifelse(grepl("gaba", d), "GABAergic",
  ifelse(grepl("^DGN", gs), "Epilepsy-DGN",
  ifelse(grepl("torc|mtor", d), "mTOR", "Other-GO/Reactome")))
}
cmp[, family := fam(gs, description)]
fsum <- cmp[, .(n = .N,
                n_sig_base   = sum(p_emp_fdr < 0.05),
                n_sig_ext    = sum(p_emp_fdr_ext < 0.05),
                med_beta     = round(median(beta), 6),
                med_t_ext    = round(median(t_ext), 3)), by = family]
fwrite(fsum, file.path(outdir, "family_specificity.csv"))
prog("Part C 家族汇总:"); print(fsum)
wt <- suppressWarnings(wilcox.test(abs(beta) ~ I(family == "GABAergic"), data = cmp))
prog(sprintf("Wilcoxon |beta| GABA vs 其他: p = %.4g", wt$p.value))
wt2 <- suppressWarnings(wilcox.test(abs(beta) ~ I(family == "Epilepsy-DGN"), data = cmp))
prog(sprintf("Wilcoxon |beta| DGN vs 其他: p = %.4g", wt2$p.value))

## ---------- Part B: 剂量反应 (代表集) ----------
reps <- sig[order(p_emp_fdr)][c(1, .N)]  # 最强 GABA + 最强 DGN 各 1
if (!any(grepl("1902710|0004890", reps$gs))) {
  g1 <- sig[grepl("1902710|0004890", gs)][1]
  if (!is.na(g1$gs)) reps <- rbind(reps, g1)
}
reps <- unique(reps, by = "gs")
dose <- dat[, .(n_grp = cut(n_asm, c(0, 1, 2, 3, Inf), labels = c("1", "2", "3", ">=4")))]
dose <- cbind(dose, as.data.table(Sv[, reps$gs, drop = FALSE]))
dm <- melt(dose, id.vars = "n_grp")
dsum <- dm[, .(mean = mean(value, na.rm = TRUE),
               se   = sd(value, na.rm = TRUE) / sqrt(.N),
               n    = .N), by = .(variable, n_grp)]
fwrite(dsum, file.path(outdir, "dose_response_representative.csv"))

## ---------- 汇总文本 ----------
lines <- c("================================================================",
           "节点3rd: vol x Number_ASMs 机制/混杂分解 — 汇总",
           sprintf("运行: %s | 注: IDEAS 仅发布 ASM 数量 (无药名), 药物机制分层不可行;",
                   format(Sys.time())),
           "本分析以 (A) 严重度全调整 + (B) 剂量反应 + (C) 基因家族特异性 回答同一问题。",
           "================================================================", "",
           sprintf("A. 严重度全调整 (cov: age_scan+sex+duration+fqFUS+fqFBTCS+SE+pathology; n=%d, FL 5000):", nrow(dat)),
           sprintf("   原 28 显著集: 方向保持 %d/%d | 调整后 p_emp<0.05: %d/%d | 全 86 集 FDR<0.05: %d",
                   sum(sig$dir_kept), nrow(sig), sum(sig$p_emp_ext < 0.05), nrow(sig),
                   sum(cmp$p_emp_fdr_ext < 0.05)), "",
           "C. 家族特异性 (base 模型 beta):")
for (i in seq_len(nrow(fsum))) lines <- c(lines, sprintf("   %-18s n=%3d 显著base=%2d 显著ext=%2d med_beta=%.2e med_t=%.2f",
                                                         fsum$family[i], fsum$n[i], fsum$n_sig_base[i], fsum$n_sig_ext[i], fsum$med_beta[i], fsum$med_t_ext[i]))
lines <- c(lines, "",
  sprintf("   Wilcoxon |beta| GABA vs 其他 p=%.3g | DGN vs 其他 p=%.3g", wt$p.value, wt2$p.value), "",
  "B. 剂量反应代表集见 dose_response_representative.csv 与 fig7 面板 b.",
  "================================================================")
writeLines(lines, file.path(outdir, "asm_mechanism_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
saveRDS(list(cmp = cmp, fsum = fsum, dsum = dsum), file.path(outdir, "asm_results_cache.rds"))
prog("34 DONE")
