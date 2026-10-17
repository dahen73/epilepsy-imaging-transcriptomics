# 16_targeted_area_vol.R -- area / vol 两模态靶向管线 (复用 86 靶向集)
# A. region 级 age_onset 梯度 | B. 组水平 FCDII/HS/DNT brainenrich | C. 个体 6 预测变量 spin 5000
suppressPackageStartupMessages({ library(BrainEnrich); library(data.table); library(Matrix) })
options(stringsAsFactors = FALSE)
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
logfile <- file.path(ws, "full_run/progress_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE); flush.console() }
prog("TARGETED SULC: start")

## ---- 基因数据 + 靶向集 (复用) ----
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname <- rownames(gene_data); universe <- colnames(gene_data)
ts <- readRDS(file.path(ws, "full_run/targeted/targeted_sets.rds"))
targeted <- ts$sets; tdesc <- ts$desc; nset <- length(targeted); ng <- ncol(gene_data)
ext2 <- new.env(); for (g in unique(unlist(targeted))) ext2[[g]] <- names(targeted)[sapply(targeted, function(s) g %in% s)]
annoT <- new.env(); annoT$EXTID2PATHID <- ext2; annoT$PATHID2EXTID <- targeted; annoT$PATHID2NAME <- tdesc
data(perm_id_dk_lh_5000); perm_id <- perm_id_dk_lh_5000; n_perm <- 5000
gn <- as.matrix(gene_data); Zg <- scale(gn)
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
  list(t = co / se, p = 2 * pt(-abs(co / se), df = df_r), beta = co)
}
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
rownames(M) <- names(targeted); colnames(M) <- universe

md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)

for (mod in c("sulc")) {
  outdir <- file.path(ws, paste0("full_run/targeted_", mod))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  prog(sprintf("=== %s ===", mod))
  indiv <- readRDS(file.path(ws, sprintf("ideas_prep_%s/brain_data_indiv_lh.rds", mod)))
  indiv <- indiv[gname, , drop = FALSE]
  idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
  keep[keep] <- complete.cases(t(indiv[, keep, drop = FALSE]))   # 剔除含 NA 区的受试者 (sulc 有 1 例)
  indiv <- indiv[, keep, drop = FALSE]; mdc <- md[idx[keep], ]
  memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
  fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)
  covC <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))

  ## Part A: region 级 age_onset
  okA <- complete.cases(data.frame(pred = mdc$age_onset, covC))
  reg_out <- do.call(rbind, lapply(seq_len(nrow(indiv)), function(r) {
    dat <- data.frame(z = as.numeric(indiv[r, okA]), pred = mdc$age_onset[okA], covC[okA, , drop = FALSE])
    f <- summary(lm(z ~ pred + age_scan + sex, data = dat))$coefficients
    data.frame(region = rownames(indiv)[r], beta = f["pred", "Estimate"], t = f["pred", "t value"], p = f["pred", "Pr(>|t|)"])
  }))
  reg_out$p_fdr <- p.adjust(reg_out$p, "BH")
  fwrite(reg_out, file.path(outdir, "region_age_onset.csv"))
  gmA <- colMeans(indiv[, okA, drop = FALSE])
  datA <- data.frame(z = gmA[okA], pred = mdc$age_onset[okA], covC[okA, , drop = FALSE])
  fg <- summary(lm(z ~ pred + age_scan + sex, data = datA))$coefficients
  prog(sprintf("%s A: %d/34 fdr<0.05, min p=%.3g (%s) | global mean p=%.3g",
               mod, sum(reg_out$p_fdr < 0.05), min(reg_out$p), reg_out$region[which.min(reg_out$p)], fg["pred", "Pr(>|t|)"]))

  ## Part B: 组水平靶向 brainenrich
  subgroups <- list(FCDII = fcdii, HS = mdc$pathology == "HS", DNT = mdc$pathology == "DNT")
  for (sg in names(subgroups)) {
    sel <- subgroups[[sg]]; sel[is.na(sel)] <- FALSE
    csvB <- file.path(outdir, sprintf("group_%s_targeted.csv", sg))
    if (sum(sel) < 10) { prog(sprintf("%s B %s: n=%d skipped", mod, sg, sum(sel))); next }
    if (file.exists(csvB)) { prog(sprintf("%s B %s: exists", mod, sg)); next }
    mapdf <- data.frame(dev = rowMeans(indiv[, sel, drop = FALSE]), row.names = gname)
    res <- brainenrich(brain_data = mapdf, gene_data = gene_data, annoData = annoT,
                       cor_method = "pearson", aggre_method = "mean",
                       null_model = "spin_brain", perm_id = perm_id, n_perm = n_perm,
                       pvalueCutoff = 1, padjCutoff = NULL)
    df <- tryCatch(as.data.frame(res), error = function(e) as.data.frame(res@result))
    fwrite(df, csvB); saveRDS(res, file.path(outdir, sprintf("group_%s_targeted.rds", sg)))
    prog(sprintf("%s B %s DONE: n=%d, fdr<0.05: %d", mod, sg, sum(sel), sum(df$p.adjust < 0.05)))
  }

  ## Part C: 个体 6 预测变量
  make_preds <- function(m, fii) {
    p <- list()
    p$Number_ASMs <- data.frame(pred = as.numeric(m$Number_ASMs))
    p$age_onset   <- data.frame(pred = m$age_onset)
    p$HS_vs_other <- data.frame(pred = factor(ifelse(is.na(m$pathology), NA,
                          ifelse(as.character(m$pathology) == "HS", "HS", "other"))))
    p$SE          <- data.frame(pred = factor(m$n_SE, levels = c("0", "1")))
    p$ILAE_Y1     <- data.frame(pred = factor(ifelse(is.na(m$ILAE_Year1), NA,
                          ifelse(m$ILAE_Year1 == 1, "seizure_free", "not_free"))))
    p$FCDII_vs_other <- data.frame(pred = factor(ifelse(fii, "FCDII", "other")))
    p
  }
  Gobs <- corr_fast(indiv)
  Sobs <- t(as.matrix((M %*% Gobs) / Msum))
  colnames(Sobs) <- names(targeted); rownames(Sobs) <- colnames(indiv)
  for (pn in c("Number_ASMs", "age_onset", "HS_vs_other", "SE", "ILAE_Y1", "FCDII_vs_other")) {
    csvC <- file.path(outdir, sprintf("individual_targeted_%s.csv", pn))
    if (file.exists(csvC)) { prog(sprintf("%s C %s: exists", mod, pn)); next }
    pd <- make_preds(mdc, fcdii)[[pn]]
    dat <- cbind(pd, covC); ok <- complete.cases(dat)
    if (length(unique(na.omit(dat$pred[ok]))) < 2) { prog(sprintf("%s C WARN %s constant", mod, pn)); next }
    X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
    term <- grep("^pred", colnames(X), value = TRUE)[1]
    tt <- lmfit_t(X, Sobs[ok, , drop = FALSE])
    ki <- match(term, rownames(tt$t))
    t_obs <- as.numeric(tt$t[ki, ]); p_par <- as.numeric(tt$p[ki, ]); beta <- as.numeric(tt$beta[ki, ])
    cnt <- integer(nset); t0 <- Sys.time()
    for (i in seq_len(n_perm)) {
      Gp <- corr_fast(indiv[perm_id[, i], , drop = FALSE])
      Sp <- t(as.matrix((M %*% Gp) / Msum))
      tp <- lmfit_t(X, Sp[ok, , drop = FALSE])$t[ki, ]
      cnt <- cnt + (abs(tp) >= abs(t_obs))
      if (i %% 1000 == 0) prog(sprintf("%s C %s: %d/%d (%.0f%%) %.1f min",
                                       mod, pn, i, n_perm, 100 * i / n_perm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    p_emp <- (cnt + 1) / (n_perm + 1); p_emp_fdr <- p.adjust(p_emp, "BH")
    out <- data.frame(gene_set = names(targeted), description = tdesc, pred = pn,
                      beta = beta, t_obs = t_obs, p_param = p_par,
                      p_emp = p_emp, p_emp_fdr = p_emp_fdr,
                      p_param_fdr = p.adjust(p_par, "BH"), N = sum(ok))
    fwrite(out, csvC)
    prog(sprintf("%s C %s DONE: p_emp<0.05: %d | fdr<0.05: %d", mod, pn, sum(p_emp < 0.05), sum(p_emp_fdr < 0.05)))
  }
}
prog("TARGETED SULC ALL COMPLETE")
