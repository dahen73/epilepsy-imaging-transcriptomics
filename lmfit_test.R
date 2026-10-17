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
X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])[1:50, ]
Y <- t(indiv[, 1:50])   # 50 x 34
fit <- lm.fit(X, Y)
res <- fit$residuals; df_r <- nrow(X) - ncol(X)
sig <- sqrt(colSums(res^2) / df_r)
xtxinv <- chol2inv(qr.R(fit$qr))
se <- sqrt(outer(diag(xtxinv), sig^2))
tval <- fit$coefficients / se
con <- file(file.path(ws, "lmfit_test.txt"), open = "wt")
writeLines(c(
  sprintf("rownames(coefficients): %s", paste(rownames(fit$coefficients), collapse = ",")),
  sprintf("rownames(tval): %s", paste(rownames(tval), collapse = ",")),
  sprintf("class tval: %s ; dim: %s", class(tval)[1], paste(dim(tval), collapse = "x")),
  sprintf("match('pred', rownames(tval)): %s", match("pred", rownames(tval)))
), con); close(con)
