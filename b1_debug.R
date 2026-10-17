suppressPackageStartupMessages({library(BrainEnrich); library(data.table); library(Matrix)})
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn <- as.matrix(gene_data); Zg <- scale(gn)
indiv <- readRDS(file.path(ws, "ideas_prep/brain_data_indiv_lh.rds"))
indiv <- indiv[rownames(gene_data), , drop = FALSE]
md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
indiv <- indiv[, keep, drop = FALSE]; md <- md[idx[keep], ]
cov_df <- data.frame(age_scan = md$age_scan, sex = factor(md$sex))

# pred_setup 完整复刻 (与主脚本一致)
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
pred_setup <- list()
for (pn in PRED_NAMES) {
  pd <- make_preds(md)[[pn]]
  dat <- cbind(pd, cov_df)
  ok  <- complete.cases(dat)
  X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
  pred_setup[[pn]] <- list(rows = which(ok), X = X,
                           term = grep("^pred", colnames(X), value = TRUE)[1])
}

corr_fast <- function(brain_mat) {
  Zb <- scale(brain_mat)
  atanh(crossprod(Zg, Zb) / (nrow(Zb) - 1))
}
lmfit_t <- function(X, Y) {
  fit  <- lm.fit(X, Y)
  res  <- fit$residuals
  n    <- nrow(X); k <- ncol(X); df_r <- n - k
  sig  <- sqrt(colSums(res^2) / df_r)
  xtxinv <- chol2inv(qr.R(fit$qr))
  se   <- sqrt(outer(diag(xtxinv), sig^2))
  tval <- fit$coefficients / se
  pval <- 2 * pt(-abs(tval), df = df_r)
  list(t = tval, p = pval, beta = fit$coefficients, df_r = df_r)
}

anno <- get_annoData("DGN")
gs   <- filter_geneSetList(colnames(gn), get_geneSetList(anno), minGSSize = 10, maxGSSize = 200)
M <- Matrix(0, nrow = length(gs), ncol = ncol(gn), sparse = TRUE)
for (i in seq_along(gs)) { j <- match(gs[[i]], colnames(gn)); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ncol(gn), 1, sparse = TRUE))
S <- t(as.matrix((M %*% corr_fast(indiv)) / Msum))

con <- file(file.path(ws, "b1_debug.txt"), open = "wt")
for (pn in names(pred_setup)) {
  ps <- pred_setup[[pn]]
  Y  <- S[ps$rows, , drop = FALSE]
  tt <- lmfit_t(ps$X, Y)
  ki <- match(ps$term, rownames(tt$t))
  writeLines(c(
    sprintf("pred=%s term=%s N=%d dim(Y)=%dx%d dim(tt$t)=%dx%d",
            pn, ps$term, nrow(ps$X), nrow(Y), ncol(Y), nrow(tt$t), ncol(tt$t)),
    sprintf("  rownames(tt$t): %s", paste(rownames(tt$t), collapse = ",")),
    sprintf("  ki: %s", ki)
  ), con)
}
close(con)
