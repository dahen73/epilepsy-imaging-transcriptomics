# =====================================================================
# 25_individual_resectexcl.R -- 个体水平切除区敏感性 (2026-10-15)
#
# 直接检验簇1 (area x FCDII_vs_other, 个体水平 18 集 fdr<0.05) 是否由
# 术后切除区伪影驱动: 剔除"任一 FCDII 患者皮层切除比例 > 0.5"的 19 个 DK 区,
# 仅在保留的 15 区上重算 86 集个体分数 (Sobs15), 重跑 FCDII_vs_other
# lm + spin 5000, 与原始 (10 号) 结果对比。
#
# 注: 组水平 FCDII 萎缩图富集在原始分析中本无显著集 (min p.adj=0.192),
#     故组水平无需剔除区重跑; 个体水平才是簇1的信号所在。
# 输出: full_run/robustness/s4_individual_FCDII_exclany50.csv
# =====================================================================
suppressPackageStartupMessages({
  library(BrainEnrich); library(data.table); library(Matrix)
})
options(stringsAsFactors = FALSE)
ws      <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir  <- file.path(ws, "full_run/robustness"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(outdir, "robustness_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE); flush.console() }
prog("INDIV-RESECTEXCL: start")

gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn <- as.matrix(gene_data); Zg <- scale(gn)
gname <- rownames(gene_data)               # 34 区名
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
indiv <- readRDS(file.path(ws, "ideas_prep_area/brain_data_indiv_lh.rds"))
indiv <- indiv[gname, , drop = FALSE]
md    <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
indiv <- indiv[, keep, drop = FALSE]; mdc <- md[idx[keep], ]
data(perm_id_dk_lh_5000); perm_id <- perm_id_dk_lh_5000; n_perm <- 5000

tgt <- readRDS(file.path(ws, "full_run/targeted/targeted_sets.rds"))
targeted <- tgt$sets; tdesc <- tgt$desc
universe <- colnames(gene_data)
nset <- length(targeted); ng <- ncol(gn)
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
prog(sprintf("LOAD: %d regions x %d genes; %d subjects; %d sets", nrow(indiv), ng, ncol(indiv), nset))

memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)

## 剔除区: 任一 FCDII 患者皮层切除比例 > 0.5
rt <- fread("F:/DJ/table_resected.csv", header = TRUE, check.names = FALSE)
rnames <- gsub("'", "", as.character(rt[[1]]))
mat <- as.matrix(rt[, -1]); colnames(mat) <- colnames(rt)[-1]; rownames(mat) <- rnames
ids_fcdii <- colnames(indiv)[fcdii]
ctx_ids <- intersect(ids_fcdii, colnames(mat))
excl_mat <- matrix(NA_real_, nrow(indiv), length(ctx_ids), dimnames = list(rownames(indiv), ctx_ids))
for (r in rownames(indiv)) {
  reg <- sub("^[LR]_", "", r); hem <- substr(r, 1, 1)
  rn0 <- grep(paste0("^ctx-", ifelse(hem == "L", "lh", "rh"), "-", reg, "$"), rownames(mat), value = TRUE)
  if (length(rn0)) excl_mat[r, ] <- suppressWarnings(as.numeric(mat[rn0, ctx_ids]))
}
any50 <- rowSums(excl_mat > 0.5, na.rm = TRUE) > 0
keepmask <- !any50; nkeep <- sum(keepmask)
Zk <- Zg[keepmask, , drop = FALSE]
prog(sprintf("EXCLUSION: %d regions excluded, %d kept", sum(any50), nkeep))

## 保留区上的个体分数 (15 区下 |r| 可达 1, clamp 防 atanh 发散)
gkeep <- function(bm) {
  r <- crossprod(Zk, scale(bm[keepmask, , drop = FALSE])) / (nkeep - 1)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  atanh(r)
}
Gobs <- gkeep(indiv)                                    # 6513 x N
Sobs <- t(as.matrix((M %*% Gobs) / Msum))               # N x 86
colnames(Sobs) <- names(targeted); rownames(Sobs) <- colnames(indiv)
prog(sprintf("Sobs15 ready: %d x %d", nrow(Sobs), ncol(Sobs)))

## FCDII_vs_other lm + spin
pred <- factor(ifelse(fcdii, "FCDII", "other"))
cov  <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))
dat  <- cbind(data.frame(pred = pred), cov); ok <- complete.cases(dat)
X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
term <- grep("^pred", colnames(X), value = TRUE)[1]
Yo <- Sobs[ok, , drop = FALSE]
tt <- lmfit_t(X, Yo); ki <- match(term, rownames(tt$t))
t_obs <- as.numeric(tt$t[ki, ]); p_par <- as.numeric(tt$p[ki, ]); beta <- as.numeric(tt$beta[ki, ])
cnt <- integer(nset); t0 <- Sys.time()
for (i in seq_len(n_perm)) {
  Gp <- gkeep(indiv[perm_id[, i], , drop = FALSE])
  Sp <- t(as.matrix((M %*% Gp) / Msum))[ok, , drop = FALSE]
  tp <- lmfit_t(X, Sp)$t[ki, ]
  cnt <- cnt + (abs(tp) >= abs(t_obs))
  if (i %% 500 == 0) prog(sprintf("  perms %d/%d (%.0f%%) elapsed %.1f min",
                                  i, n_perm, 100 * i / n_perm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
p_emp <- (cnt + 1) / (n_perm + 1)
out <- data.frame(gene_set = names(targeted), description = tdesc,
                  beta_excl = beta, t_obs = t_obs, p_param = p_par,
                  p_emp = p_emp, p_emp_fdr = p.adjust(p_emp, "BH"), N = sum(ok))
fwrite(out, file.path(outdir, "s4_individual_FCDII_exclany50.csv"))

## 对比原始 (10 号 area individual FCDII_vs_other)
orig <- fread(file.path(ws, "full_run/targeted_area/individual_targeted_FCDII_vs_other.csv"))
m <- merge(out, orig[, .(gene_set, beta_orig = beta, p_emp_orig = p_emp, p_emp_fdr_orig = p_emp_fdr)],
           by = "gene_set")
m$dir_kept <- sign(m$beta_excl) == sign(m$beta_orig)
m$sig_orig <- m$p_emp_fdr_orig < 0.05
fwrite(m, file.path(outdir, "s4_individual_FCDII_exclany50_compare.csv"))
prog(sprintf("S4 DONE: %d orig-sig (fdr<0.05) | dir kept %d/%d | still fdr<0.05 %d/%d | %.1f min",
             sum(m$sig_orig), sum(m$dir_kept[m$sig_orig]), sum(m$sig_orig),
             sum(m$p_emp_fdr[m$sig_orig] < 0.05, na.rm = TRUE), sum(m$sig_orig),
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
prog("INDIV-RESECTEXCL COMPLETE")
