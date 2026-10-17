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
dat <- data.frame(pred = as.numeric(md$Number_ASMs), cov_df)
ok <- complete.cases(dat)
ri <- which(ok)[1:50]
X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])[1:50, ]

corr_fast <- function(brain_mat) {
  Zb <- scale(brain_mat)
  atanh(crossprod(Zg, Zb) / (nrow(Zb) - 1))
}
anno_v <- get_annoData("DGN")
gs_v <- get_geneSetList(anno_v)
gs_v <- filter_geneSetList(colnames(gn), gs_v, minGSSize = 10, maxGSSize = 200)
M <- Matrix(0, nrow = length(gs_v), ncol = ncol(gn), sparse = TRUE)
for (i in seq_along(gs_v)) { j <- match(gs_v[[i]], colnames(gn)); M[i, j[!is.na(j)]] <- 1 }
G_fast <- corr_fast(indiv)
Msum <- as.numeric(M %*% Matrix(1, ncol(gn), 1, sparse = TRUE))
Yv <- t(as.matrix(M[1:50, ] %*% G_fast[, ri, drop = FALSE]))
fit <- lm.fit(X, Yv)
res <- fit$residuals; df_r <- nrow(X) - ncol(X)
sig <- sqrt(colSums(res^2) / df_r)
xtxinv <- chol2inv(qr.R(fit$qr))
se <- sqrt(outer(diag(xtxinv), sig^2))
t_fast <- (fit$coefficients / se)[2, ]
sl <- simple_lm(dependent_df = as.data.frame(Yv), pred_df = data.frame(pred = md$Number_ASMs[ri]),
                cov_df = cov_df[ri, , drop = FALSE], stat2return = "tval")
t_pkg <- sl$tval[match(colnames(Yv), sl$Dependent_vars)]
con <- file(file.path(ws, "t_diff2.txt"), open = "wt")
writeLines(c(
  sprintf("max abs diff: %g", max(abs(t_fast - t_pkg), na.rm = TRUE)),
  sprintf("n diff>1e-6: %d / 50", sum(abs(t_fast - t_pkg) > 1e-6, na.rm = TRUE)),
  sprintf("t_fast[1:3]: %s", paste(sprintf("%.6f", t_fast[1:3]), collapse = " | ")),
  sprintf("t_pkg [1:3]: %s", paste(sprintf("%.6f", t_pkg[1:3]), collapse = " | ")),
  sprintf("match NA: %d", sum(is.na(match(colnames(Yv), sl$Dependent_vars))))
), con); close(con)
