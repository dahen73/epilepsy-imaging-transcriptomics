# =====================================================================
# 50_r1m1_resection_site.R -- R1-M1 resolution test (2026-10-23)
# 问题: 三个个体水平信号簇是否由"切除部位"(侧别/脑叶)混杂驱动?
# 设计: 在 24_robustness.R S1 (切除负担) 基础上, 逐级加入切除侧别与切除脑叶
#   M0: score ~ pred + age_scan + sex                       (原模型)
#   M1: M0 + resect_burden                                  (= S1)
#   M2: M1 + op_side                                        (切除侧别 L/R)
#   M3: M2 + resect_lobe                                    (切除脑叶 T/F/P/O/H/multi)
# 另加 (iii): 非 HS 患者按切除脑叶分层, 簇内方向一致性
# 输入: 本目录 inputs/ 下已复制的 rds/csv (源自 D:\DJ, 只读不改)
# 输出: 本目录 outputs/
# =====================================================================
options(stringsAsFactors = FALSE)
ws  <- "C:/Epilepsia_Submission/analysis_r1m1"
inp <- file.path(ws, "inputs")
out <- file.path(ws, "outputs"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(out, "r1m1_log.txt")
prog <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logfile, append = TRUE) }
set.seed(20261023)

## ---- 最小 lm 工具 (与 24_robustness.R 相同的矩阵化实现) ----
lmfit_t <- function(X, Y) {
  Y <- as.matrix(Y); fit <- lm.fit(X, Y); res <- fit$residuals
  if (is.null(dim(res))) res <- matrix(res, ncol = 1)
  co <- fit$coefficients
  if (is.null(dim(co))) co <- matrix(co, ncol = 1, dimnames = list(colnames(X), NULL))
  n <- nrow(X); k <- ncol(X); df_r <- n - k
  sig <- sqrt(colSums(res^2) / df_r)
  xtxinv <- chol2inv(qr.R(fit$qr))
  se <- sqrt(outer(diag(xtxinv), sig^2))
  tval <- co / se; pval <- 2 * pt(-abs(tval), df = df_r)
  list(t = tval, p = pval, beta = co)
}
run_lm <- function(X, Y) {
  tt <- lmfit_t(X, Y); ki <- grep("^pred", rownames(tt$t))[1]
  data.frame(beta = as.numeric(tt$beta[ki, ]), t = as.numeric(tt$t[ki, ]), p = as.numeric(tt$p[ki, ]))
}

## ---- 数据 ----
md  <- readRDS(file.path(inp, "phenotypes.rds")); md$ID <- as.character(md$ID)
bdf <- read.csv(file.path(inp, "resection_burden.csv")); bdf$ID <- as.character(bdf$ID)
Sobs_all <- list(area = readRDS(file.path(inp, "Sobs_area.rds")),
                 vol  = readRDS(file.path(inp, "Sobs_vol.rds")))

## 切除部位协变量: 侧别 + 脑叶 (由 Op_Type 派生)
parse_lobe <- function(ot) {
  ot <- trimws(ot)
  if (is.na(ot) || ot == "") return(NA_character_)
  if (ot == "Hx") return("Hemispheric")
  toks <- strsplit(ot, " ")[[1]]
  lobes <- unique(toks[toks %in% c("T", "F", "P", "O")])
  if (length(lobes) == 0) return(NA_character_)
  if (length(lobes) > 1) return("Multi")
  c(T = "Temporal", F = "Frontal", P = "Parietal", O = "Occipital")[lobes]
}
md$resect_lobe <- vapply(md$Op_Type, parse_lobe, character(1))
md$op_side     <- factor(md$Op_Side)
md$resect_lobe <- factor(md$resect_lobe, levels = c("Temporal", "Frontal", "Parietal", "Occipital", "Multi", "Hemispheric"))
prog(sprintf("lobe dist: %s", paste(names(table(md$resect_lobe)), as.integer(table(md$resect_lobe)), collapse = " ")))

## 组装数据框 (与 24_robustness.R build_dat 一致, 另加侧别/脑叶)
build_dat <- function(mod, pred) {
  S <- Sobs_all[[mod]]
  mdc <- md[match(rownames(S), md$ID), ]
  memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
  fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)
  pv <- switch(pred,
    Number_ASMs    = data.frame(pred = as.numeric(mdc$Number_ASMs)),
    age_onset      = data.frame(pred = mdc$age_onset),
    FCDII_vs_other = data.frame(pred = factor(ifelse(fcdii, "FCDII", "other"))))
  bur <- data.frame(resect_burden = bdf$burden_mean[match(rownames(S), bdf$ID)])
  dat <- cbind(pv,
               data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex)),
               bur,
               data.frame(op_side = mdc$op_side, resect_lobe = mdc$resect_lobe,
                          pathology = mdc$pathology))
  rownames(dat) <- rownames(S)
  dat
}

## ---- 全部显著集 (fdr<0.05), 与正文一致 ----
combos <- list(
  list(mod = "area", pred = "FCDII_vs_other", csv = file.path(inp, "individual_targeted_FCDII_vs_other.csv")),
  list(mod = "vol",  pred = "Number_ASMs",    csv = file.path(inp, "individual_targeted_Number_ASMs.csv")),
  list(mod = "vol",  pred = "age_onset",      csv = file.path(inp, "individual_targeted_age_onset.csv")))
sig_sets <- list()
for (cb in combos) {
  d <- read.csv(cb$csv)
  d <- d[!is.na(d$p_emp_fdr) & d$p_emp_fdr < 0.05, ]
  d <- d[order(d$p_emp_fdr), ]
  sig_sets[[length(sig_sets) + 1]] <- list(mod = cb$mod, pred = cb$pred, sets = d$gene_set)
  prog(sprintf("SIG %s x %s: %d sets (fdr<0.05)", cb$mod, cb$pred, nrow(d)))
}

## ---- Part A: 逐级协变量模型 ----
prog("PART A: hierarchical covariate models M0->M3")
rows <- list()
for (rp in sig_sets) {
  dat <- build_dat(rp$mod, rp$pred)
  ok  <- complete.cases(dat[, c("pred", "age_scan", "sex", "resect_burden", "op_side", "resect_lobe")])
  S   <- Sobs_all[[rp$mod]][ok, rp$sets, drop = FALSE]
  d   <- dat[ok, , drop = FALSE]
  X0 <- model.matrix(~ pred + age_scan + sex, data = d)
  X1 <- model.matrix(~ pred + age_scan + sex + resect_burden, data = d)
  X2 <- model.matrix(~ pred + age_scan + sex + resect_burden + op_side, data = d)
  X3 <- model.matrix(~ pred + age_scan + sex + resect_burden + op_side + resect_lobe, data = d)
  r0 <- run_lm(X0, S); r1 <- run_lm(X1, S); r2 <- run_lm(X2, S); r3 <- run_lm(X3, S)
  for (k in seq_along(rp$sets)) {
    rows[[length(rows) + 1]] <- data.frame(
      mod = rp$mod, pred = rp$pred, gene_set = rp$sets[k], N = sum(ok),
      beta_M0 = r0$beta[k], p_M0 = r0$p[k],
      beta_M1 = r1$beta[k], p_M1 = r1$p[k],
      beta_M2 = r2$beta[k], p_M2 = r2$p[k],
      beta_M3 = r3$beta[k], p_M3 = r3$p[k],
      dir_kept_M3 = sign(r0$beta[k]) == sign(r3$beta[k]),
      sig05_kept_M3 = (r0$p[k] < 0.05) == (r3$p[k] < 0.05),
      beta_change_M0toM3_pct = 100 * (r3$beta[k] - r0$beta[k]) / abs(r0$beta[k]))
  }
  prog(sprintf("A done: %s x %s (N=%d, %d sets)", rp$mod, rp$pred, sum(ok), length(rp$sets)))
}
resA <- do.call(rbind, rows)
write.csv(resA, file.path(out, "r1m1_covariate_models.csv"), row.names = FALSE)
prog(sprintf("PART A TOTAL: %d tests | direction kept M3 %d/%d | p<0.05 status kept %d/%d | median |beta change| %.1f%% (max %.1f%%)",
             nrow(resA), sum(resA$dir_kept_M3), nrow(resA), sum(resA$sig05_kept_M3), nrow(resA),
             median(abs(resA$beta_change_M0toM3_pct)), max(abs(resA$beta_change_M0toM3_pct))))

## FDR 水平稳健性: M3 模型下按簇内 BH 重算, 看 fdr<0.05 集保留率
for (rp in sig_sets) {
  sub <- resA[resA$mod == rp$mod & resA$pred == rp$pred, ]
  fdr3 <- p.adjust(sub$p_M3, method = "BH")
  prog(sprintf("FDR recheck %s x %s: base p<0.05 %d/%d | M3 p<0.05 %d/%d | M3 BH-FDR<0.05 %d/%d",
               rp$mod, rp$pred, sum(sub$p_M0 < 0.05), nrow(sub), sum(sub$p_M3 < 0.05), nrow(sub), sum(fdr3 < 0.05), nrow(sub)))
}

## ---- Part B: (iii) 非 HS 按切除脑叶分层 ----
prog("PART B: lobe-stratified replication in non-HS patients")
rowsB <- list()
for (rp in sig_sets) {
  dat <- build_dat(rp$mod, rp$pred)
  ok  <- complete.cases(dat[, c("pred", "age_scan", "sex", "resect_burden", "op_side", "resect_lobe")])
  dat <- dat[ok, , drop = FALSE]; S <- Sobs_all[[rp$mod]][ok, rp$sets, drop = FALSE]
  nonhs <- dat$pathology != "HS"
  for (lb in c("Temporal", "Frontal")) {
    sel <- nonhs & dat$resect_lobe == lb
    if (sum(sel) < 30) { prog(sprintf("B skip %s x %s lobe=%s: n=%d", rp$mod, rp$pred, lb, sum(sel))); next }
    Xs <- model.matrix(~ pred + age_scan + sex, data = dat[sel, , drop = FALSE])
    rs <- run_lm(Xs, S[sel, , drop = FALSE])
    ## 与全样本 M0 方向比较
    Xa <- model.matrix(~ pred + age_scan + sex, data = dat)
    ra <- run_lm(Xa, S)
    for (k in seq_along(rp$sets)) {
      rowsB[[length(rowsB) + 1]] <- data.frame(
        mod = rp$mod, pred = rp$pred, gene_set = rp$sets[k], lobe_stratum = lb, N = sum(sel),
        beta_full = ra$beta[k], p_full = ra$p[k],
        beta_stratum = rs$beta[k], p_stratum = rs$p[k],
        dir_kept = sign(ra$beta[k]) == sign(rs$beta[k]))
    }
  }
  prog(sprintf("B done: %s x %s", rp$mod, rp$pred))
}
resB <- do.call(rbind, rowsB)
write.csv(resB, file.path(out, "r1m1_lobe_stratified.csv"), row.names = FALSE)
if (!is.null(resB)) {
  agg <- aggregate(dir_kept ~ mod + pred + lobe_stratum, data = resB, FUN = function(x) c(mean = mean(x), n = length(x)))
  print(agg)
  prog(sprintf("PART B TOTAL: direction kept %.1f%% (%d/%d strata-tests)",
               100 * mean(resB$dir_kept), sum(resB$dir_kept), nrow(resB)))
}
prog("R1M1 COMPLETE")
