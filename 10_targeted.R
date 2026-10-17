# =====================================================================
# 10_targeted.R -- 靶向假设驱动分析 (2026-10-14)
#
# A. Region 级 age_onset 梯度: 34 区逐一 lm(z ~ age_onset + age_scan + sex)
#    + 全局混淆检查 (平均厚度 vs age_onset)
# B. 组水平靶向富集: FCDII / HS / DNT 亚组萎缩图 x 靶向基因集, brainenrich spin 5000
# C. 个体水平靶向检验: 靶向集分数 x 6 预测变量 (参数化 + spin 5000 置换)
#
# 靶向基因集: 从已装库 (KEGG/Reactome/GO_BP/GO_CC/GO_MF/SynGO/DGN) 按先验挖掘
#   - mTOR 通路 (FCDII = mTORopathy)
#   - GABA 能突触/受体 (HS = GABA 能神经元丢失)
#   - 癫痫/发作疾病基因 (DGN: epilep|seizure)
# 输出: full_run/targeted/
# =====================================================================
suppressPackageStartupMessages({
  library(BrainEnrich); library(data.table); library(Matrix)
})
options(stringsAsFactors = FALSE)
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir <- file.path(ws, "full_run/targeted"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(ws, "full_run/progress_log.txt")
prog <- function(...) { msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE); flush.console() }
prog("TARGETED: start")

## ---- 载入 (与 04 相同) ----
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname     <- rownames(gene_data)
indiv <- readRDS(file.path(ws, "ideas_prep/brain_data_indiv_lh.rds"))
indiv <- indiv[gname, , drop = FALSE]
md    <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
indiv <- indiv[, keep, drop = FALSE]; mdc <- md[idx[keep], ]
stopifnot(ncol(indiv) == nrow(mdc))
data(perm_id_dk_lh_5000)
perm_id <- perm_id_dk_lh_5000; n_perm <- 5000
prog(sprintf("LOAD: %d regions x %d genes; %d subjects", nrow(indiv), ncol(gene_data), ncol(indiv)))

## ---- FCDII 提取 (OP MEMO: type IIA/IIB/II; 排除 IIID) ----
memo <- ifelse(is.na(mdc[["OP MEMO"]]), "", as.character(mdc[["OP MEMO"]]))
fcdii <- mdc$pathology == "FCD" & grepl("^type II($|[AB])", memo)
prog(sprintf("FCDII: n=%d (of %d FCD, %d total)", sum(fcdii), sum(mdc$pathology == "FCD", na.rm = TRUE), nrow(mdc)))

## ---- 靶向基因集挖掘 ----
libs <- c("KEGG", "Reactome", "GO_BP", "GO_CC", "GO_MF", "SynGO", "DGN")
pat_mtor <- "mTOR|mechanistic target|rapamycin|TORC1|TORC2|TOR signaling"
pat_gaba <- "GABA"
pat_epi  <- "epilep|seizure"
universe <- colnames(gene_data)
targeted <- list(); tdesc <- character(0)
for (lb in libs) {
  anno <- get_annoData(lb)
  sets <- anno$PATHID2EXTID
  nms  <- anno$PATHID2NAME
  for (tm in names(sets)) {
    lab <- if (!is.null(nms) && tm %in% names(nms)) nms[[tm]] else ""
    hay <- paste(tm, lab)
    pat <- if (lb == "DGN") pat_epi else paste(pat_mtor, "|", pat_gaba)
    if (!grepl(pat, hay, ignore.case = (lb == "DGN"))) next
    g <- intersect(sets[[tm]], universe)
    if (length(g) < 10 || length(g) > 200) next
    key <- paste0(lb, ":", tm)
    targeted[[key]] <- g
    tdesc[key] <- paste0(lab, " (", lb, ")")
  }
}
prog(sprintf("TARGETED SETS: %d mined", length(targeted)))
for (k in names(targeted)) prog(sprintf("  %s | n=%d | %s", k, length(targeted[[k]]), tdesc[k]))
stopifnot(length(targeted) > 0)

## 自定义 annoData 环境 (DOSE 结构)
ext2 <- new.env(); for (g in unique(unlist(targeted))) ext2[[g]] <- names(targeted)[sapply(targeted, function(s) g %in% s)]
annoT <- new.env()
annoT$EXTID2PATHID <- ext2
annoT$PATHID2EXTID <- targeted
annoT$PATHID2NAME  <- tdesc
saveRDS(list(sets = targeted, desc = tdesc), file.path(outdir, "targeted_sets.rds"))

## ---- 通用工具 (与 04 等价) ----
gn  <- as.matrix(gene_data); Zg <- scale(gn)
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
ng <- ncol(gn); nset <- length(targeted)
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
rownames(M) <- names(targeted); colnames(M) <- universe

## ============================================================
## Part A: region 级 age_onset 梯度
## ============================================================
prog("PART A: region-level age_onset gradient")
covA <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))
okA  <- complete.cases(data.frame(pred = mdc$age_onset, covA))
reg_out <- do.call(rbind, lapply(seq_len(nrow(indiv)), function(r) {
  dat <- data.frame(z = as.numeric(indiv[r, okA]), pred = mdc$age_onset[okA], covA[okA, , drop = FALSE])
  f <- summary(lm(z ~ pred + age_scan + sex, data = dat))$coefficients
  data.frame(region = rownames(indiv)[r], beta = f["pred", "Estimate"], t = f["pred", "t value"], p = f["pred", "Pr(>|t|)"])
}))
reg_out$p_fdr <- p.adjust(reg_out$p, "BH")
fwrite(reg_out, file.path(outdir, "region_age_onset.csv"))
gmA <- colMeans(indiv[, okA, drop = FALSE])
datA <- data.frame(z = gmA[okA], pred = mdc$age_onset[okA], covA[okA, , drop = FALSE])
fg <- summary(lm(z ~ pred + age_scan + sex, data = datA))$coefficients
prog(sprintf("A global-confound: mean-thickness ~ age_onset beta=%.4f p=%.3g", fg["pred", "Estimate"], fg["pred", "Pr(>|t|)"]))
prog(sprintf("A done: %d/34 regions p_fdr<0.05; min p=%.3g (%s)",
             sum(reg_out$p_fdr < 0.05), min(reg_out$p), reg_out$region[which.min(reg_out$p)]))

## ============================================================
## Part B: 组水平靶向富集 (FCDII / HS / DNT 萎缩图)
## ============================================================
prog("PART B: group-level targeted brainenrich")
subgroups <- list(FCDII = fcdii, HS = mdc$pathology == "HS", DNT = mdc$pathology == "DNT")
for (sg in names(subgroups)) {
  sel <- subgroups[[sg]]; sel[is.na(sel)] <- FALSE
  if (sum(sel) < 10) { prog(sprintf("B %s: n=%d skipped", sg, sum(sel))); next }
  csvB <- file.path(outdir, sprintf("group_%s_targeted.csv", sg))
  if (file.exists(csvB)) { prog(sprintf("B %s: exists, skipped", sg)); next }
  mapdf <- data.frame(dev = rowMeans(indiv[, sel, drop = FALSE]), row.names = gname)
  t0 <- Sys.time()
  res <- brainenrich(brain_data = mapdf, gene_data = gene_data, annoData = annoT,
                     cor_method = "pearson", aggre_method = "mean",
                     null_model = "spin_brain", perm_id = perm_id, n_perm = n_perm,
                     pvalueCutoff = 1, padjCutoff = NULL)
  df <- tryCatch(as.data.frame(res), error = function(e) as.data.frame(res@result))
  fwrite(df, file.path(outdir, sprintf("group_%s_targeted.csv", sg)))
  saveRDS(res, file.path(outdir, sprintf("group_%s_targeted.rds", sg)))
  prog(sprintf("B %s DONE: n=%d, %d sets, p.adjust<0.05: %d | %.1f min",
               sg, sum(sel), nrow(df), sum(df$p.adjust < 0.05), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ============================================================
## Part C: 个体水平靶向检验 (参数化 + spin 5000)
## ============================================================
prog("PART C: individual-level targeted tests")
make_preds <- function(m) {
  p <- list()
  p$Number_ASMs <- data.frame(pred = as.numeric(m$Number_ASMs))
  p$age_onset   <- data.frame(pred = m$age_onset)
  p$HS_vs_other <- data.frame(pred = factor(ifelse(is.na(m$pathology), NA,
                        ifelse(as.character(m$pathology) == "HS", "HS", "other"))))
  p$SE          <- data.frame(pred = factor(m$n_SE, levels = c("0", "1")))
  p$ILAE_Y1     <- data.frame(pred = factor(ifelse(is.na(m$ILAE_Year1), NA,
                        ifelse(m$ILAE_Year1 == 1, "seizure_free", "not_free"))))
  p$FCDII_vs_other <- data.frame(pred = factor(ifelse(fcdii, "FCDII", "other")))
  p
}
covC <- data.frame(age_scan = mdc$age_scan, sex = factor(mdc$sex))
Gobs <- corr_fast(indiv)                                # 6513 x N
Sobs <- t(as.matrix((M %*% Gobs) / Msum))               # N x nset
colnames(Sobs) <- names(targeted)
rownames(Sobs) <- colnames(indiv)
prog(sprintf("C scores ready: %d x %d", nrow(Sobs), ncol(Sobs)))

for (pn in c("Number_ASMs", "age_onset", "HS_vs_other", "SE", "ILAE_Y1", "FCDII_vs_other")) {
  pd <- make_preds(mdc)[[pn]]
  dat <- cbind(pd, covC); ok <- complete.cases(dat)
  if (length(unique(na.omit(dat$pred[ok]))) < 2) { prog(sprintf("C WARN %s constant", pn)); next }
  X <- model.matrix(~ pred + age_scan + sex, data = dat[ok, , drop = FALSE])
  term <- grep("^pred", colnames(X), value = TRUE)[1]
  Yo <- Sobs[ok, , drop = FALSE]
  tt <- lmfit_t(X, Yo)
  ki <- match(term, rownames(tt$t))
  t_obs <- as.numeric(tt$t[ki, ]); p_par <- as.numeric(tt$p[ki, ]); beta <- as.numeric(tt$beta[ki, ])
  csvC <- file.path(outdir, sprintf("individual_targeted_%s.csv", pn))
  if (file.exists(csvC)) { prog(sprintf("C %s: exists, skipped", pn)); next }
  cnt <- integer(nset)
  t0 <- Sys.time()
  for (i in seq_len(n_perm)) {
    Gp <- corr_fast(indiv[perm_id[, i], , drop = FALSE])  # 重排脑区行 (spin_brain)
    Sp <- t(as.matrix((M %*% Gp) / Msum))
    Sp <- Sp[ok, , drop = FALSE]
    tp <- lmfit_t(X, Sp)$t[ki, ]
    cnt <- cnt + (abs(tp) >= abs(t_obs))
    if (i %% 250 == 0) prog(sprintf("C %s: %d/%d perms (%.0f%%) elapsed %.1f min",
                                    pn, i, n_perm, 100 * i / n_perm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  p_emp <- (cnt + 1) / (n_perm + 1)
  p_emp_fdr <- p.adjust(p_emp, "BH")
  out <- data.frame(gene_set = names(targeted), description = tdesc, pred = pn,
                    beta = beta, t_obs = t_obs, p_param = p_par,
                    p_emp = p_emp, p_emp_fdr = p.adjust(p_emp, "BH"), p_param_fdr = p.adjust(p_par, "BH"), N = sum(ok))
  fwrite(out, file.path(outdir, sprintf("individual_targeted_%s.csv", pn)))
  prog(sprintf("C %s DONE: %d sets | p_emp<0.05: %d | p_emp_fdr<0.05: %d | %.1f min",
               pn, nset, sum(p_emp < 0.05), sum(p_emp_fdr < 0.05), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
prog("TARGETED ALL COMPLETE")
