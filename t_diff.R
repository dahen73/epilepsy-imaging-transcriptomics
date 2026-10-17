suppressPackageStartupMessages({library(BrainEnrich); library(data.table); library(Matrix)})
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn <- as.matrix(gene_data)
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
stopifnot(length(unique(md$sex[ri])) == 2)

Y <- t(indiv[, ri])   # 50 x 34 直接用脑区分数代替基因集(等价测试 lm 机制)
fit <- lm.fit(X, Y)
res <- fit$residuals; n <- nrow(X); k <- ncol(X); df_r <- n - k
sig <- sqrt(colSums(res^2) / df_r)
xtxinv <- chol2inv(qr.R(fit$qr))
se <- sqrt(outer(diag(xtxinv), sig^2))
t_fast <- fit$coefficients / se
sl <- simple_lm(dependent_df = as.data.frame(Y), pred_df = data.frame(pred = md$Number_ASMs[ri]),
                cov_df = cov_df[ri, , drop = FALSE], stat2return = "tval")
t_pkg <- sl$tval[match(colnames(Y), sl$Dependent_vars)]
con <- file(file.path(ws, "t_diff.txt"), open = "wt")
writeLines(c(
  sprintf("N=%d, df_r=%d", n, df_r),
  sprintf("max abs diff (pred t): %g", max(abs(t_fast[2, ] - t_pkg), na.rm = TRUE)),
  sprintf("t_fast[2,1:3]: %s", paste(sprintf("%.8f", t_fast[2, 1:3]), collapse = " | ")),
  sprintf("t_pkg [1:3]  : %s", paste(sprintf("%.8f", t_pkg[1:3]), collapse = " | ")),
  sprintf("X colnames: %s", paste(colnames(X), collapse = ",")),
  sprintf("any NA t_fast: %d ; t_pkg: %d", sum(is.na(t_fast[2, ])), sum(is.na(t_pkg)))
), con); close(con)
