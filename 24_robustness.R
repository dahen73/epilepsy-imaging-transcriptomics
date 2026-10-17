# =====================================================================
# 24_robustness.R -- IDEAS 内部稳健性分析 (2026-10-15)
#
# 目标: 检验三大 FDR 信号簇的稳健性
#   簇1: area x FCDII_vs_other   (18 集癫痫基因, fdr=0.020)
#   簇2: vol   x Number_ASMs     (28 集 GABA+癫痫DGN, 全正向)
#   簇3: vol   x age_onset       (2 集 + region 级左额梯度)
#
# Part 0: 环境 + 四模态 Sobs (受试者 x 86 靶向集分数矩阵, 首次持久化)
# Part 1: 切除负担计算 (table_resected.csv -> 每患者皮层切除负担)
# Part 2: S1 敏感性 -- 三大簇代表集 lm 加入 resect_burden 协变量
# Part 3: S2 留一法 jackknife -- 逐受试者剔除, beta/p 稳定性
# Part 4: S3 随机亚样本重复 -- 80% 子样本 x 200, 方向/显著性稳定率
#
# 注: 组水平切除区剔除见 25_group_resectexcl.R (spin 5000, 慢, 单独跑)
# 输出: full_run/robustness/
# =====================================================================
suppressPackageStartupMessages({
  library(BrainEnrich); library(data.table); library(Matrix)
})
options(stringsAsFactors = FALSE)
ws      <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir  <- file.path(ws, "full_run/robustness"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(outdir, "robustness_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE); flush.console() }
prog("ROBUSTNESS: start")
set.seed(20261015)

## ---------------- Part 0: 环境 ----------------
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname     <- rownames(gene_data)
gn        <- as.matrix(gene_data); Zg <- scale(gn)
corr_fast <- function(brain_mat) atanh(crossprod(Zg, scale(brain_mat)) / (ncol(brain_mat) - 1))
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
  list(t = tval, p = pval, beta = co, df_r = df_r)
}

## 靶向集 (与 10 号完全一致: 直接读已存 rds)
tgt <- readRDS(file.path(ws, "full_run/targeted/targeted_sets.rds"))
targeted <- tgt$sets; tdesc <- tgt$desc
nset <- length(targeted); universe <- colnames(gene_data)
prog(sprintf("LOAD: targeted sets=%d", nset))

ng <- ncol(gn)
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
rownames(M) <- names(targeted); colnames(M) <- universe

## 四模态 indiv 矩阵 + 各自匹配后的元数据
mods <- list(
  thickness = list(dir = "ideas_prep",      rds = "brain_data_indiv_lh.rds"),
  area      = list(dir = "ideas_prep_area", rds = "brain_data_indiv_lh.rds"),
  vol       = list(dir = "ideas_prep_vol",  rds = "brain_data_indiv_lh.rds"),
  sulc      = list(dir = "ideas_prep_sulc", rds = "brain_data_indiv_lh.rds")
)
md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
Sobs_all <- list()
for (mn in names(mods)) {
  sv <- file.path(outdir, sprintf("Sobs_%s.rds", mn))
  if (file.exists(sv)) { Sobs_all[[mn]] <- readRDS(sv); next }
  indiv <- readRDS(file.path(ws, mods[[mn]]$dir, mods[[mn]]$rds))
  indiv <- indiv[gname, , drop = FALSE]
  idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
  indiv <- indiv[, keep, drop = FALSE]
  Gobs <- corr_fast(indiv)
  S <- t(as.matrix((M %*% Gobs) / Msum))
  dimnames(S) <- list(colnames(indiv), names(targeted))
  saveRDS(S, sv)
  Sobs_all[[mn]] <- S
  prog(sprintf("Sobs %s: %d subj x %d sets (saved)", mn, nrow(S), ncol(S)))
}
for (mn in names(Sobs_all)) prog(sprintf("Sobs %s: %d subj (existing)", mn, nrow(Sobs_all[[mn]])))
prog("Part 0 done")

## ---------------- Part 1: 切除负担 ----------------
resect_csv <- "/f/DJ/table_resected.csv"
resect_csv <- "F:/DJ/table_resected.csv"
burden_sv <- file.path(outdir, "resection_burden.csv")
if (!file.exists(burden_sv)) {
  rt <- fread(resect_csv, header = TRUE, check.names = FALSE)
  rnames <- as.character(rt[[1]]); rnames <- gsub("'", "", rnames)
  mat <- as.matrix(rt[, -1]); colnames(mat) <- colnames(rt)[-1]; rownames(mat) <- rnames
  ctx <- grep("^ctx-(lh|rh)-", rownames(mat), value = TRUE)      # 68 皮层区
  prog(sprintf("RESECT: table %d regions x %d patients; %d ctx rows", nrow(mat), ncol(mat), length(ctx)))
  ## 映射到 Sobs 受试者 ID
  ## 每患者: 皮层 68 区切除比例 -> mean / max / n>0.5
  bfun <- function(pid) {
    if (!pid %in% colnames(mat)) return(c(NA, NA, NA))
    v <- suppressWarnings(as.numeric(mat[ctx, pid]))
    v <- v[!is.na(v)]
    c(if (length(v)) mean(v) else NA, if (length(v)) max(v) else NA, sum(v > 0.5))
  }
  ids <- rownames(Sobs_all$thickness)   # 行 = 受试者 ID; 切勿用 colnames (那是基因集名)
  bb <- t(vapply(ids, bfun, numeric(3)))
  colnames(bb) <- c("burden_mean", "burden_max", "n_removed_gt50")
  bdf <- data.frame(ID = ids, bb, row.names = NULL)
  bdf$is_control <- as.numeric(bdf$ID) >= 4000
  bdf$burden_mean[bdf$is_control] <- 0; bdf$burden_max[bdf$is_control] <- 0; bdf$n_removed_gt50[bdf$is_control] <- 0
  fwrite(bdf, burden_sv)
  prog(sprintf("RESECT: burden computed for %d subjects; patient mean burden=%.4f, max=%.3f, n>50%% removed: %d",
               nrow(bdf), mean(bdf$burden_mean[!bdf$is_control], na.rm = TRUE),
               max(bdf$burden_max, na.rm = TRUE), sum(bdf$n_removed_gt50 > 0, na.rm = TRUE)))
}
bdf <- fread(burden_sv); bdf$ID <- as.character(bdf$ID)
prog("Part 1 done")

## ---------------- 代表集选取 (三大簇) ----------------
sel_file <- file.path(outdir, "cluster_rep_sets.csv")
combos <- list(
  list(mod = "area", pred = "FCDII_vs_other", csv = file.path(ws, "full_run/targeted_area/individual_targeted_FCDII_vs_other.csv")),
  list(mod = "vol",  pred = "Number_ASMs",    csv = file.path(ws, "full_run/targeted_vol/individual_targeted_Number_ASMs.csv")),
  list(mod = "vol",  pred = "age_onset",      csv = file.path(ws, "full_run/targeted_vol/individual_targeted_age_onset.csv"))
)
rep_sets <- list()
if (file.exists(sel_file)) {
  seldf <- fread(sel_file)
  for (i in seq_len(nrow(seldf))) rep_sets[[i]] <- list(mod = seldf$mod[i], pred = seldf$pred[i], sets = strsplit(seldf$sets[i], ";")[[1]])
} else {
  for (cb in combos) {
    d <- fread(cb$csv)
    sig <- d[order(d$p_emp_fdr), ]
    sig <- sig[p_emp_fdr < 0.05]
    if (nrow(sig) == 0) { prog(sprintf("REPSET WARN %s x %s: no fdr<0.05", cb$mod, cb$pred)); next }
    take <- head(sig$gene_set, 5)
    rep_sets[[length(rep_sets) + 1]] <- list(mod = cb$mod, pred = cb$pred, sets = take)
    prog(sprintf("REPSET %s x %s: %d sig -> top%d", cb$mod, cb$pred, nrow(sig), length(take)))
  }
  fwrite(data.frame(mod = sapply(rep_sets, `[[`, "mod"), pred = sapply(rep_sets, `[[`, "pred"),
                    sets = sapply(rep_sets, function(x) paste(x$sets, collapse = ";"))), sel_file)
}
prog(sprintf("REPRSETS: %d combos", length(rep_sets)))

## 通用: 组装 lm 数据框 (mdc 按 S 行名即时匹配, 避免跨模态受试者集不同导致的错位)
build_dat <- function(mod, pred) {
  S <- Sobs_all[[mod]]
  mdc <- md[match(rownames(S), md$ID), ]
  memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
  fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)
  pv <- switch(pred,
    Number_ASMs    = data.frame(pred = as.numeric(mdc$Number_ASMs)),
    age_onset      = data.frame(pred = mdc$age_onset),
    FCDII_vs_other = data.frame(pred = factor(ifelse(fcdii, "FCDII", "other")))
  )
  cov <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))
  bur <- data.frame(resect_burden = bdf$burden_mean[match(rownames(S), bdf$ID)])
  dat <- cbind(pv, cov, bur); rownames(dat) <- rownames(S)
  dat
}
run_lm <- function(X, Y) {
  tt <- lmfit_t(X, Y); ki <- grep("^pred", rownames(tt$t))[1]
  data.frame(beta = as.numeric(tt$beta[ki, ]), t_obs = as.numeric(tt$t[ki, ]), p_param = as.numeric(tt$p[ki, ]))
}

## ---------------- Part 2: S1 切除负担敏感性 ----------------
prog("PART 2: resection-burden sensitivity")
s1rows <- list()
for (rp in rep_sets) {
  dat <- build_dat(rp$mod, rp$pred)
  ok <- complete.cases(dat)
  S <- Sobs_all[[rp$mod]][ok, rp$sets, drop = FALSE]
  X0 <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
  X1 <- model.matrix(~ pred + age_scan + sex + resect_burden, data = dat[ok, , drop = FALSE])
  r0 <- run_lm(X0, S); r1 <- run_lm(X1, S)
  for (k in seq_along(rp$sets)) {
    s1rows[[length(s1rows) + 1]] <- data.frame(
      mod = rp$mod, pred = rp$pred, gene_set = rp$sets[k],
      N = sum(ok),
      beta_base = r0$beta[k],  p_base = r0$p_param[k],
      beta_adj  = r1$beta[k],  p_adj  = r1$p_param[k],
      beta_change_pct = 100 * (r1$beta[k] - r0$beta[k]) / abs(r0$beta[k]),
      direction_kept = sign(r0$beta[k]) == sign(r1$beta[k]),
      sig_kept_p05 = (r0$p_param[k] < 0.05) == (r1$p_param[k] < 0.05))
  }
}
s1 <- rbindlist(s1rows)
fwrite(s1, file.path(outdir, "s1_resection_sensitivity.csv"))
prog(sprintf("S1 done: %d tests | direction kept %d/%d | p<0.05 status kept %d/%d | median |beta change|=%.1f%%",
             nrow(s1), sum(s1$direction_kept), nrow(s1), sum(s1$sig_kept_p05), nrow(s1),
             median(abs(s1$beta_change_pct))))

## ---------------- Part 3: S2 留一法 jackknife ----------------
prog("PART 3: leave-one-out jackknife")
jrows <- list()
for (rp in rep_sets) {
  dat <- build_dat(rp$mod, rp$pred)
  ok <- complete.cases(dat)
  dat_ok <- dat[ok, , drop = FALSE]
  S <- Sobs_all[[rp$mod]][ok, rp$sets, drop = FALSE]
  Xf <- model.matrix(~ pred + age_scan + sex, data = dat_ok)
  base <- run_lm(Xf, S)
  n <- nrow(dat_ok)
  jk_beta <- matrix(NA_real_, n, length(rp$sets), dimnames = list(NULL, rp$sets))
  jk_p    <- jk_beta
  for (i in seq_len(n)) {
    Xi <- model.matrix(~ pred + age_scan + sex, data = dat_ok[-i, , drop = FALSE])
    ri <- run_lm(Xi, S[-i, , drop = FALSE])
    jk_beta[i, ] <- ri$beta; jk_p[i, ] <- ri$p_param
  }
  for (k in seq_along(rp$sets)) {
    jrows[[length(jrows) + 1]] <- data.frame(
      mod = rp$mod, pred = rp$pred, gene_set = rp$sets[k], N = n,
      beta_base = base$beta[k], p_base = base$p_param[k],
      jk_beta_mean = mean(jk_beta[, k]), jk_beta_sd = sd(jk_beta[, k]),
      jk_beta_min = min(jk_beta[, k]), jk_beta_max = max(jk_beta[, k]),
      jk_sign_flips = sum(sign(jk_beta[, k]) != sign(base$beta[k])),
      jk_p_lt05_prop = mean(jk_p[, k] < 0.05))
  }
  prog(sprintf("S2 %s x %s: done (%d subjects)", rp$mod, rp$pred, n))
}
jk <- rbindlist(jrows)
fwrite(jk, file.path(outdir, "s2_jackknife.csv"))
prog(sprintf("S2 done: sign flips total %d | median jk_p<0.05 prop=%.3f",
             sum(jk$jk_sign_flips), median(jk$jk_p_lt05_prop)))

## ---------------- Part 4: S3 随机亚样本重复 (80% x 200) ----------------
prog("PART 4: subsample stability (80% x 200)")
n_boot <- 200
brows <- list()
for (rp in rep_sets) {
  dat <- build_dat(rp$mod, rp$pred)
  ok <- complete.cases(dat)
  dat_ok <- dat[ok, , drop = FALSE]
  S <- Sobs_all[[rp$mod]][ok, rp$sets, drop = FALSE]
  n <- nrow(dat_ok); nsub <- floor(0.8 * n)
  bs_beta <- matrix(NA_real_, n_boot, length(rp$sets), dimnames = list(NULL, rp$sets))
  bs_p    <- bs_beta
  for (b in seq_len(n_boot)) {
    ii <- sample(seq_len(n), nsub)
    Xb <- model.matrix(~ pred + age_scan + sex, data = dat_ok[ii, , drop = FALSE])
    rb <- run_lm(Xb, S[ii, , drop = FALSE])
    bs_beta[b, ] <- rb$beta; bs_p[b, ] <- rb$p_param
  }
  for (k in seq_along(rp$sets)) {
    brows[[length(brows) + 1]] <- data.frame(
      mod = rp$mod, pred = rp$pred, gene_set = rp$sets[k],
      beta_base = run_lm(model.matrix(~ pred + age_scan + sex, data = dat_ok), S)$beta[k],
      boot_beta_mean = mean(bs_beta[, k]), boot_beta_sd = sd(bs_beta[, k]),
      boot_beta_ci_lo = quantile(bs_beta[, k], 0.025), boot_beta_ci_hi = quantile(bs_beta[, k], 0.975),
      boot_dir_consistency = mean(sign(bs_beta[, k]) == sign(mean(bs_beta[, k]))),
      boot_p_lt05_prop = mean(bs_p[, k] < 0.05), n_boot = n_boot, n_sub = nsub)
  }
  prog(sprintf("S3 %s x %s: done", rp$mod, rp$pred))
}
bs <- rbindlist(brows)
fwrite(bs, file.path(outdir, "s3_subsample_stability.csv"))
prog(sprintf("S3 done: median CI coverage of base beta=%d/%d | median boot_p<0.05 prop=%.3f",
             sum(bs$beta_base >= bs$boot_beta_ci_lo & bs$beta_base <= bs$boot_beta_ci_hi), nrow(bs),
             median(bs$boot_p_lt05_prop)))

prog("ROBUSTNESS COMPLETE (fast parts)")
