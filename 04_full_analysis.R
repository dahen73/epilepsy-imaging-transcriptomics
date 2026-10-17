# 04_full_analysis.R
# IDEAS x 转录组学: 全样本全参数正式分析 (幂等可续跑)
#
# 阶段 A  : 组水平 brainenrich() -- 10 基因集库 x spin_brain 5000
#           + DGN 敏感性: spearman / median / resample_gene
# 阶段 B1 : 个体水平参数化全景 -- 10 库 x 5 预测变量 (null = none)
# 阶段 B2 : 个体水平置换检验 -- DGN + SynGO, spin_brain 5000 x 5 预测变量
#           (增量分块实现: 每置换即聚合->回归->仅保留 t 值计数, 内存 O(50MB))
#
# 进度: full_run/progress_log.txt (追加) + full_run/status.txt (当前状态)
# 结果: full_run/group/ | full_run/individual_param/ | full_run/individual/

suppressPackageStartupMessages({
  library(BrainEnrich); library(data.table); library(Matrix)
})
options(stringsAsFactors = FALSE)

ws    <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
run   <- file.path(ws, "full_run")
dir.create(file.path(run, "group"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run, "individual_param"), showWarnings = FALSE)
dir.create(file.path(run, "individual"), showWarnings = FALSE)
LOG   <- file.path(run, "progress_log.txt")
ST    <- file.path(run, "status.txt")

prog <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n"); flush.console()
  con <- file(LOG, open = "a", encoding = "UTF-8")
  writeLines(line, con); close(con)
  con2 <- file(ST, open = "w", encoding = "UTF-8")
  writeLines(line, con2); close(con2)
}

## ---------------- 公共数据 ----------------
prog("LOAD: gene expression (desikan L) + inputs")
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname     <- rownames(gene_data)

group_df <- as.data.frame(readRDS(file.path(ws, "ideas_prep/brain_data_group_lh.rds")))
stopifnot(identical(rownames(group_df), gname))   # 34 区, 已按 gene_data 排序

indiv <- readRDS(file.path(ws, "ideas_prep/brain_data_indiv_lh.rds"))
indiv <- indiv[gname, , drop = FALSE]             # 行序对齐
md    <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)

data(perm_id_dk_lh_5000)                          # 34 x 5000 (包内置)
perm_id <- perm_id_dk_lh_5000
n_perm <- 5000
prog(sprintf("LOAD ok: %d regions x %d genes; %d subjects; perm_id %d",
             nrow(gene_data), ncol(gene_data), ncol(indiv), ncol(perm_id_dk_lh_5000)))

LIBS <- c("DGN", "SynGO", "GO_BP", "GO_MF", "GO_CC", "KEGG",
          "Reactome", "CellTypes_Seidlitz2020", "CellTypes_Martins2021", "WikiPathways")

## ============================================================
## 阶段 A: 组水平 brainenrich()
## ============================================================
ga <- function(lib, cor_m, agg_m, null_m) {
  sprintf("%s_%s_%s_%s", lib, cor_m, agg_m, null_m)
}
group_configs <- list()
for (lb in LIBS) group_configs[[length(group_configs) + 1]] <-
  list(lib = lb, cor = "pearson", agg = "mean", null = "spin_brain")
group_configs[[length(group_configs) + 1]] <- list(lib = "DGN", cor = "spearman", agg = "mean",    null = "spin_brain")
group_configs[[length(group_configs) + 1]] <- list(lib = "DGN", cor = "pearson",  agg = "median",  null = "spin_brain")
group_configs[[length(group_configs) + 1]] <- list(lib = "DGN", cor = "pearson",  agg = "mean",    null = "resample_gene")

prog(sprintf("STAGE A: group-level, %d configs", length(group_configs)))
for (ci in seq_along(group_configs)) {
  cfg <- group_configs[[ci]]
  tag <- ga(cfg$lib, cfg$cor, cfg$agg, cfg$null)
  csv <- file.path(run, "group", paste0(tag, ".csv"))
  rds <- file.path(run, "group", paste0(tag, "_gseaResult.rds"))
  if (file.exists(csv)) { prog(sprintf("A %d/%d SKIP (done): %s", ci, length(group_configs), tag)); next }
  prog(sprintf("A %d/%d START: %s", ci, length(group_configs), tag))
  anno <- get_annoData(cfg$lib)
  t0 <- Sys.time()
  res <- brainenrich(
    brain_data = group_df, gene_data = gene_data, annoData = anno,
    cor_method = cfg$cor, aggre_method = cfg$agg,
    null_model = cfg$null, perm_id = perm_id_dk_lh_5000,
    n_perm = n_perm, pvalueCutoff = 1, padjCutoff = NULL)
  df <- tryCatch(as.data.frame(res), error = function(e) as.data.frame(res@result))
  fwrite(df, csv)
  saveRDS(res, rds)
  nrow_sig <- sum(df$p.adjust < 0.05)
  prog(sprintf("A %d/%d DONE : %s | %d sets | p.adj<0.05: %d | %.1f min",
               ci, length(group_configs), tag, nrow(df), nrow_sig,
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(res, df, anno); gc(verbose = FALSE)
}

## ============================================================
## 阶段 B 公共: 个体水平设计矩阵 / 分数计算 / 批量 lm
## ============================================================
mdc_idx <- match(colnames(indiv), md$ID)
if (any(is.na(mdc_idx))) {
  prog(sprintf("WARN: %d subjects in z-table missing from metadata -> excluded",
               sum(is.na(mdc_idx))))
  keep <- !is.na(mdc_idx)
  indiv <- indiv[, keep, drop = FALSE]
  mdc <- md[mdc_idx[keep], ]
} else mdc <- md[mdc_idx, ]
stopifnot(!any(is.na(mdc$ID)), ncol(indiv) == nrow(mdc))

# 5 个预测变量: 名称 -> (数值向量/因子, 说明)
make_preds <- function(m) {
  p <- list()
  p$Number_ASMs <- data.frame(pred = as.numeric(m$Number_ASMs))
  p$age_onset   <- data.frame(pred = m$age_onset)
  p$HS_vs_other <- data.frame(pred = factor(ifelse(is.na(m$pathology), NA,
                        ifelse(as.character(m$pathology) == "HS", "HS", "other"))))
  p$SE          <- data.frame(pred = factor(m$n_SE, levels = c("0", "1")))
  p$ILAE_Y1     <- data.frame(pred = factor(ifelse(is.na(m$ILAE_Year1), NA,
                        ifelse(m$ILAE_Year1 == 1, "seizure_free", "not_free"))))
  p
}
PRED_NAMES <- c("Number_ASMs", "age_onset", "HS_vs_other", "SE", "ILAE_Y1")
cov_df <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))

# 每个 pred: complete cases -> 设计矩阵 X (含截距), 记录主效应列名/行索引
pred_setup <- list()
for (pn in PRED_NAMES) {
  pd <- make_preds(mdc)[[pn]]
  dat <- cbind(pd, cov_df)
  ok  <- complete.cases(dat)
  if (length(unique(na.omit(dat$pred[ok]))) < 2) { prog(sprintf("WARN: %s constant, skipped", pn)); next }
  X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
  pred_setup[[pn]] <- list(rows = which(ok), X = X,
                           term = grep("^pred", colnames(X), value = TRUE)[1])
  prog(sprintf("PRED %s: N=%d, term=%s", pn, sum(ok), pred_setup[[pn]]$term))
}

# 标准化基因矩阵 (一次): r2z pearson 相关 = t(Zg) %*% Zb / (n-1)
gn  <- as.matrix(gene_data)                 # 34 x 6513
Zg  <- scale(gn)                            # 34 x 6513
corr_fast <- function(brain_mat) {          # brain_mat: 34 x N -> 6513 x N (atanh r)
  Zb <- scale(brain_mat)
  atanh(crossprod(Zg, Zb) / (nrow(Zb) - 1))
}

# 批量 lm t 值: Y = N x GS, X = N x k -> list(t = k x GS, beta, p)
lmfit_t <- function(X, Y) {
  Y <- as.matrix(Y)
  fit  <- lm.fit(X, Y)
  res  <- fit$residuals
  if (is.null(dim(res))) res <- matrix(res, ncol = 1)   # 单列响应时 lm.fit 退化为 vector
  co   <- fit$coefficients
  if (is.null(dim(co))) co <- matrix(co, ncol = 1, dimnames = list(colnames(X), NULL))
  n    <- nrow(X); k <- ncol(X)
  df_r <- n - k
  sig  <- sqrt(colSums(res^2) / df_r)
  xtxinv <- chol2inv(qr.R(fit$qr))          # (X'X)^-1, k x k
  se   <- sqrt(outer(diag(xtxinv), sig^2))  # k x GS
  tval <- co / se
  pval <- 2 * pt(-abs(tval), df = df_r)
  list(t = tval, p = pval, beta = co, df_r = df_r)
}

# 基因集 membership 稀疏矩阵 (GS x 6513) + 等价性校验
build_membership <- function(selected.gs) {
  M <- Matrix(0, nrow = length(selected.gs), ncol = ncol(gn), sparse = TRUE)
  for (i in seq_along(selected.gs)) {
    j <- match(selected.gs[[i]], colnames(gn))
    M[i, j[!is.na(j)]] <- 1
  }
  M
}

## ---- 校验: corr_fast + membership(mean) + lmfit_t vs 包函数 ----
prog("VERIFY: equivalence of fast path vs package functions")
anno_v <- get_annoData("DGN")
gs_v   <- get_geneSetList(anno_v)
gs_v   <- filter_geneSetList(colnames(gn), gs_v, minGSSize = 10, maxGSSize = 200)
stopifnot(length(gs_v) > 0)
nchk <- 50   # 50 例保证 sex 两水平均在, X 满秩 (前 10 例可能单水平导致秩亏)
ri <- pred_setup$Number_ASMs$rows[1:nchk]
stopifnot(length(unique(mdc$sex[ri])) == 2)
G_pkg  <- corr_brain_gene(gene_data = gene_data, brain_data = indiv[, ri], method = "pearson")
G_fast <- corr_fast(indiv)          # 全体受试者, 6513 x N
stopifnot(isTRUE(all.equal(as.vector(G_pkg), as.vector(G_fast[, ri]), tolerance = 1e-8)))
Mv <- build_membership(gs_v)
Msum_v <- as.numeric(Mv %*% Matrix(1, ncol(gn), 1, sparse = TRUE))
s_pkg  <- sapply(gs_v[1:50], function(g) mean(G_pkg[match(g, rownames(G_pkg)), 1]))
s_fast <- as.numeric((Mv[1:50, ] %*% G_fast[, ri[1]]) / Msum_v[1:50])
stopifnot(isTRUE(all.equal(unname(s_pkg), unname(s_fast), tolerance = 1e-8)))
Xv <- pred_setup$Number_ASMs$X[1:nchk, , drop = FALSE]
Yv <- t(as.matrix(Mv[1:50, ] %*% G_fast[, ri, drop = FALSE]))
t_fast <- lmfit_t(Xv, Yv)$t[grep("^pred", colnames(Xv)), ]
# 与 stats::lm (黄金标准) 逐列对比
chk_dat <- data.frame(mdc[ri, ], check.names = FALSE)
t_ref <- sapply(1:ncol(Yv), function(j) {
  d <- cbind(Yv[, j], pred = mdc$Number_ASMs[ri], chk_dat[, c("age_scan", "sex")])
  summary(lm(Yv[, j] ~ pred + age_scan + sex, data = d))$coefficients["pred", "t value"]
})
stopifnot(isTRUE(all.equal(unname(t_fast), unname(t_ref), tolerance = 1e-6)))
prog("VERIFY OK: fast path == package (corr / aggregate / lm t)")

## ============================================================
## 阶段 B1: 个体水平参数化全景 (10 库 x 5 pred, null = none)
## ============================================================
prog(sprintf("STAGE B1: individual parametric, %d libs x %d preds", length(LIBS), length(PRED_NAMES)))
for (lb in LIBS) {
  csv <- file.path(run, "individual_param", paste0(lb, "_param_all.csv"))
  if (file.exists(csv)) { prog(sprintf("B1 SKIP (done): %s", lb)); next }
  prog(sprintf("B1 START: %s", lb))
  anno <- get_annoData(lb)
  gs   <- filter_geneSetList(colnames(gn), get_geneSetList(anno), minGSSize = 10, maxGSSize = 200)
  prog(sprintf("B1 %s: %d gene sets after filter", lb, length(gs)))
  M    <- build_membership(gs)
  Msum <- as.numeric(M %*% Matrix(1, ncol(gn), 1, sparse = TRUE))
  S    <- t(as.matrix((M %*% corr_fast(indiv)) / Msum))   # N x GS empirical scores
  colnames(S) <- names(gs)          # 稀疏乘法丢失 dimnames, 显式补回
  rownames(S) <- colnames(indiv)
  out_all <- list()
  for (pn in names(pred_setup)) {
    ps <- pred_setup[[pn]]
    Y  <- S[ps$rows, , drop = FALSE]
    tt <- lmfit_t(ps$X, Y)
    ki <- match(ps$term, rownames(tt$t))
    out <- data.frame(gene_set = colnames(S), stringsAsFactors = FALSE)
    out$beta     <- as.numeric(tt$beta[ki, ])
    out$t        <- as.numeric(tt$t[ki, ])
    out$p_param  <- as.numeric(tt$p[ki, ])
    out$p_fdr <- p.adjust(out$p_param, "fdr")
    desc <- tryCatch(get_termDescription(colnames(S), anno), error = function(e) NA)
    out$description <- unname(desc[colnames(S)])
    out$pred <- pn; out$N <- nrow(Y)
    out_all[[pn]] <- out
  }
  fwrite(rbindlist(out_all, use.names = TRUE, fill = TRUE), csv)
  prog(sprintf("B1 DONE : %s | %d sets x %d preds -> %s", lb, ncol(S), length(out_all), basename(csv)))
  rm(anno, gs, M, S, out_all); gc(verbose = FALSE)
}

## ============================================================
## 阶段 B2: 个体水平置换检验 (spin_brain 5000) DGN + SynGO x 5 pred
## ============================================================
B2_LIBS <- c("DGN", "SynGO")
CHUNK   <- 250
for (lb in B2_LIBS) {
  tag <- paste0(lb, "_spin5000")
  done <- file.path(run, "individual", paste0(tag, "_all.csv"))
  if (file.exists(done)) { prog(sprintf("B2 SKIP (done): %s", lb)); next }
  prog(sprintf("B2 START: %s (spin_brain %d, chunk %d)", lb, n_perm, CHUNK))
  anno <- get_annoData(lb)
  gs   <- filter_geneSetList(colnames(gn), get_geneSetList(anno), minGSSize = 10, maxGSSize = 200)
  M    <- build_membership(gs); Msum <- as.numeric(M %*% Matrix(1, ncol(gn), 1, sparse = TRUE))
  Gobs <- corr_fast(indiv)
  Sobs <- t(as.matrix((M %*% Gobs) / Msum))          # N_all x GS
  colnames(Sobs) <- names(gs)                        # 稀疏乘法丢失 dimnames, 显式补回
  rownames(Sobs) <- colnames(indiv)
  # 观测 t 与参数化 p (per pred)
  obs <- list()
  for (pn in names(pred_setup)) {
    ps <- pred_setup[[pn]]
    tt <- lmfit_t(ps$X, Sobs[ps$rows, , drop = FALSE])
    ki <- match(ps$term, rownames(tt$t))
    obs[[pn]] <- data.frame(gene_set = colnames(Sobs), stringsAsFactors = FALSE)
    obs[[pn]]$beta    <- as.numeric(tt$beta[ki, ])
    obs[[pn]]$t_obs   <- as.numeric(tt$t[ki, ])
    obs[[pn]]$p_param <- as.numeric(tt$p[ki, ])
    obs[[pn]]$N       <- nrow(ps$X)
  }
  # null 计数 (per pred): |t_null| >= |t_obs|
  cnt <- setNames(as.data.frame(matrix(0, ncol(Sobs), length(obs)), stringsAsFactors = FALSE), names(obs))
  rownames(cnt) <- colnames(Sobs)
  n_done <- 0; t0 <- Sys.time()
  for (b in seq(1, n_perm, by = CHUNK)) {
    idx <- b:min(b + CHUNK - 1, n_perm)
    for (i in idx) {
      Gnull <- corr_fast(indiv[perm_id[, i], , drop = FALSE])   # 行(区)重排, 与包内 spin_brain 一致
      Sn    <- t(as.matrix((M %*% Gnull) / Msum))               # N x GS
      rm(Gnull)
      for (pn in names(pred_setup)) {
        ps <- pred_setup[[pn]]
        tn <- lmfit_t(ps$X, Sn[ps$rows, , drop = FALSE])$t
        ki <- match(ps$term, rownames(tn))
        cnt[[pn]] <- cnt[[pn]] + as.numeric(abs(tn[ki, ]) >= abs(obs[[pn]]$t_obs) - 1e-12)
      }
    }
    n_done <- length(idx) + n_done - 1
    prog(sprintf("B2 %s: %d/%d permutations (%.0f%%) elapsed %.1f min",
                 lb, n_done, n_perm, 100 * n_done / n_perm,
                 as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  # 汇总: 经验 p + FDR + 正态性 (KS on empirical scores, 全体)
  out_all <- list()
  ks_p <- sapply(seq_len(ncol(Sobs)), function(j) tryCatch(ks.test(scale(Sobs[, j]), "pnorm")$p.value, error = function(e) NA))
  ks_fdr <- p.adjust(ks_p, "fdr")
  for (pn in names(obs)) {
    o <- obs[[pn]]
    o$p_emp <- (cnt[[pn]] + 1) / (n_perm + 1)
    o$p_emp_fdr <- p.adjust(o$p_emp, "fdr")
    o$score_nonnormal_fdr <- ks_fdr[match(o$gene_set, colnames(Sobs))]
    o$pred <- pn
    desc <- tryCatch(get_termDescription(o$gene_set, anno), error = function(e) NA)
    o$description <- unname(desc[o$gene_set])
    out_all[[pn]] <- o
  }
  fwrite(rbindlist(out_all, use.names = TRUE, fill = TRUE), done)
  saveRDS(list(gene_sets = gs, observed_scores = Sobs, null_counts = cnt, n_perm = n_perm),
          file.path(run, "individual", paste0(tag, "_support.rds")))
  prog(sprintf("B2 DONE : %s | %d sets x %d preds -> %s", lb, ncol(Sobs), length(obs), basename(done)))
  rm(anno, gs, M, Gobs, Sobs, obs, cnt, out_all); gc(verbose = FALSE)
}

prog("ALL STAGES COMPLETE")
