# =====================================================================
# 30_individual_rh.R -- 节点2nd: 右半球个体水平敏感性 (2026-10-15)
#
# 4 模态 (thick/area/vol/sulc, rh 32 区) x 6 预测变量 x 86 靶向集
# brainscore 分数 -> lm(pred + age_scan + sex) -> 参数 t + spin 5000 (rh 置换)
# 优化: 每次置换同时算 6 个预测变量的 lm, 置换次数 5000/模态 (lh 管线为 6x5000)。
# 输出: full_run/rh/individual_rh_<mod>_<pred>.csv
# =====================================================================
suppressPackageStartupMessages({ library(BrainEnrich); library(data.table); library(Matrix) })
options(stringsAsFactors = FALSE)
ws      <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir  <- file.path(ws, "full_run/rh"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(outdir, "rh_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE); flush.console() }
prog("INDIV-RH start")

gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "R"))
gn <- as.matrix(gene_data); Zg <- scale(gn)          # 32 区 x 6513
gname <- rownames(gene_data)
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
tgt <- readRDS(file.path(ws, "full_run/targeted/targeted_sets.rds"))
targeted <- tgt$sets; tdesc <- tgt$desc
universe <- colnames(gene_data); nset <- length(targeted); ng <- ncol(gn)
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
perm_id <- readRDS(file.path(outdir, "perm_id_dk_rh_5000.rds"))
n_perm <- ncol(perm_id)
md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
prog(sprintf("LOAD: %d regions x %d genes; %d sets; %d perms", nrow(gn), ng, nset, n_perm))

mods <- list(
  thick = list(dir = "ideas_prep_rh",     rds = "brain_data_indiv_thick_rh.rds"),
  area  = list(dir = "ideas_prep_rh",     rds = "brain_data_indiv_area_rh.rds"),
  vol   = list(dir = "ideas_prep_rh",     rds = "brain_data_indiv_vol_rh.rds"),
  sulc  = list(dir = "ideas_prep_rh",     rds = "brain_data_indiv_sulc_rh.rds")
)

for (mn in names(mods)) {
  indiv <- readRDS(file.path(ws, mods[[mn]]$dir, mods[[mn]]$rds))
  indiv <- indiv[gname, , drop = FALSE]
  na_subj <- apply(indiv, 2, function(x) any(is.na(x)))
  if (any(na_subj)) {
    prog(sprintf("%s: dropping %d subjects with NA: %s", mn, sum(na_subj), paste(colnames(indiv)[na_subj], collapse = ",")))
    indiv <- indiv[, !na_subj, drop = FALSE]
  }
  idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
  indiv <- indiv[, keep, drop = FALSE]; mdc <- md[idx[keep], ]
  Gobs <- corr_fast(indiv)
  Sobs <- t(as.matrix((M %*% Gobs) / Msum))
  rownames(Sobs) <- colnames(indiv); colnames(Sobs) <- names(targeted)

  memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
  fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)
  preds <- list(
    Number_ASMs    = data.frame(pred = as.numeric(mdc$Number_ASMs)),
    age_onset      = data.frame(pred = mdc$age_onset),
    HS_vs_other    = data.frame(pred = factor(ifelse(is.na(mdc$pathology), NA, ifelse(as.character(mdc$pathology) == "HS", "HS", "other")))),
    SE             = data.frame(pred = factor(mdc$n_SE, levels = c("0", "1"))),
    ILAE_Y1        = data.frame(pred = factor(ifelse(is.na(mdc$ILAE_Year1), NA, ifelse(mdc$ILAE_Year1 == 1, "seizure_free", "not_free")))),
    FCDII_vs_other = data.frame(pred = factor(ifelse(fcdii, "FCDII", "other")))
  )
  covC <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))
  ## 每个 pred 的模型与观测统计
  P <- list()
  for (pn in names(preds)) {
    dat <- cbind(preds[[pn]], covC); ok <- complete.cases(dat)
    if (length(unique(na.omit(dat$pred[ok]))) < 2) { prog(sprintf("%s WARN %s constant", mn, pn)); next }
    X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
    term <- grep("^pred", colnames(X), value = TRUE)[1]
    tt <- lmfit_t(X, Sobs[ok, , drop = FALSE]); ki <- match(term, rownames(tt$t))
    P[[pn]] <- list(ok = ok, X = X, ki = ki,
                    t_obs = as.numeric(tt$t[ki, ]), beta = as.numeric(tt$beta[ki, ]),
                    p_param = as.numeric(tt$p[ki, ]))
  }
  cnt <- setNames(rep(0, length(P)), names(P)); t0 <- Sys.time()
  for (i in seq_len(n_perm)) {
    Gp <- corr_fast(indiv[perm_id[, i], , drop = FALSE])
    Sp <- t(as.matrix((M %*% Gp) / Msum))
    for (pn in names(P)) {
      tp <- lmfit_t(P[[pn]]$X, Sp[P[[pn]]$ok, , drop = FALSE])$t[P[[pn]]$ki, ]
      cnt[pn] <- cnt[pn] + sum(abs(tp) >= abs(P[[pn]]$t_obs))
    }
    if (i %% 1000 == 0) prog(sprintf("  %s: %d/%d perms (%.0f%%) elapsed %.1f min",
                                     mn, i, n_perm, 100 * i / n_perm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  for (pn in names(P)) {
    p_emp <- (cnt[pn] / nset + 1) / (n_perm + 1)
    out <- data.frame(gene_set = names(targeted), description = tdesc, pred = pn,
                      beta = P[[pn]]$beta, t_obs = P[[pn]]$t_obs, p_param = P[[pn]]$p_param,
                      p_emp = p_emp, p_emp_fdr = p.adjust(p_emp, "BH"),
                      p_param_fdr = p.adjust(P[[pn]]$p_param, "BH"), N = sum(P[[pn]]$ok))
    fwrite(out, file.path(outdir, sprintf("individual_rh_%s_%s.csv", mn, pn)))
    prog(sprintf("%s x %s DONE: %d sets | p_emp<0.05: %d | fdr<0.05: %d", mn, pn, nset, sum(p_emp < 0.05), sum(p_emp < 0.05 & p.adjust(p_emp, "BH") < 0.05)))
  }
  prog(sprintf("MODALITY %s COMPLETE (%.1f min)", mn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
prog("INDIV-RH ALL COMPLETE")
