# =====================================================================
# 41_stratified_pathology.R
# Pathology-spectrum stratified individual-level analysis (node: user request)
#   - Strata: HS / non-HS (primary, FDR within cell) + DNT/FCD/CAV/DUAL/OTHER
#     (small, nominal p only, power-limited; GL/TBC/TREBLE n<10 skipped)
#   - Predictors: Number_ASMs, age_onset, SE, ILAE_Y1, FCDII_vs_other
#     (auto-skip if constant or any level < 10 within stratum)
#   - 4 modalities (thickness/area/vol/sulc), spin 5000, BH-FDR within
#     stratum x modality x predictor cell across 86 targeted sets
#   - Speed: restrict permutation score computation to genes in the 86 sets
# Output: full_run/stratified_pathology/
# =====================================================================
suppressPackageStartupMessages({
  library(BrainEnrich); library(data.table); library(Matrix)
})
options(stringsAsFactors = FALSE)
ws     <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir <- file.path(ws, "full_run/stratified_pathology")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(outdir, "log.txt")
prog <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n"); cat(line, "\n", file = logfile, append = TRUE)
}

## ---- load gene expression + targeted sets (same as 10/24) ----
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn  <- as.matrix(gene_data); universe <- colnames(gene_data)
Zg  <- scale(gn)                                   # 34 x 6513, gene-standardized across regions
corr_fast <- function(brain_mat) atanh(crossprod(Zg, scale(brain_mat)) / (ncol(brain_mat) - 1))

tgt <- readRDS(file.path(ws, "full_run/targeted/targeted_sets.rds"))
targeted <- tgt$sets; tdesc <- tgt$desc
nset <- length(targeted); ng <- ncol(gn)
M <- Matrix(0, nrow = nset, ncol = ng, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], universe); M[i, j[!is.na(j)]] <- 1 }
Msum <- as.numeric(M %*% Matrix(1, ng, 1, sparse = TRUE))
rownames(M) <- names(targeted); colnames(M) <- universe

## ---- speed-up: restrict to genes present in the 86 sets ----
needed <- sort(unique(unlist(targeted)))
needed <- needed[needed %in% universe]
nG <- length(needed)
ZgN <- Zg[, needed, drop = FALSE]                  # 34 x nG
Mn <- Matrix(0, nrow = nset, ncol = nG, sparse = TRUE)
for (i in seq_len(nset)) { j <- match(targeted[[i]], needed); Mn[i, j[!is.na(j)]] <- 1 }
MsumN <- as.numeric(Mn %*% Matrix(1, nG, 1, sparse = TRUE))
prog(sprintf("LOAD: %d sets, %d/%d genes needed", nset, nG, ng))

## ---- per-stratum observed scores via fast path (identical math to corr_fast) ----
fast_scores <- function(indiv, ok_cols) {
  S <- scale(indiv[, ok_cols, drop = FALSE])       # 34 x n (column = subject, as in base pipeline)
  GpN <- atanh(pmin(pmax(crossprod(ZgN, S) / (ncol(S) - 1), -0.999999), 0.999999))  # nG x n
  Y <- t(as.matrix((Mn %*% GpN) / MsumN))          # n x 86
  colnames(Y) <- names(targeted); rownames(Y) <- ok_cols
  Y
}
## verification: fast path must equal corr_fast route on full thickness sample
indiv0 <- readRDS(file.path(ws, "ideas_prep", "brain_data_indiv_lh.rds"))
indiv0 <- indiv0[rownames(gene_data), , drop = FALSE]
Sobs_ref <- readRDS(file.path(ws, "full_run/robustness/Sobs_thickness.rds"))
ids0 <- intersect(colnames(indiv0), rownames(Sobs_ref))
Ycheck <- fast_scores(indiv0, ids0)
dmax <- max(abs(Ycheck - Sobs_ref[rownames(Ycheck), ]), na.rm = TRUE)
prog(sprintf("VERIFY fast vs corr_fast route: max abs diff = %.2e", dmax))
if (dmax > 1e-6) { prog("VERIFY FAILED - aborting"); quit(status = 1) }

## ---- per-stratum observed scores via fast path (identical math to corr_fast) ----
lmfit_t <- function(X, Y) {
  Y <- as.matrix(Y); fit <- lm.fit(X, Y); res <- fit$residuals
  if (is.null(dim(res))) res <- matrix(res, ncol = 1)
  co <- fit$coefficients
  if (is.null(dim(co))) co <- matrix(co, ncol = 1, dimnames = list(colnames(X), NULL))
  n <- nrow(X); k <- ncol(X); df_r <- n - k
  sig <- sqrt(colSums(res^2) / df_r)
  xtxinv <- chol2inv(qr.R(fit$qr))
  se <- sqrt(outer(diag(xtxinv), sig^2))
  tval <- co / se
  list(t = tval, beta = co, df_r = df_r)
}

## ---- metadata + strata ----
md <- readRDS(file.path(ws, "ideas_prep/phenotypes.rds")); md$ID <- as.character(md$ID)
fcdii_id <- md$ID[!is.na(md$`OP MEMO`) & grepl("^type II($|[AB])", md$`OP MEMO`)]
md$path4 <- ifelse(is.na(md$pathology), NA, as.character(md$pathology))
strata <- list(
  HS    = "HS", nonHS = c("DNT", "FCD", "CAV", "DUAL", "OTHER", "GL", "TBC", "TREBLE"),
  DNT   = "DNT", FCD = "FCD", CAV = "CAV", DUAL = "DUAL", OTHER = "OTHER"
)
mods <- list(
  thickness = list(dir = "ideas_prep",      rds = "brain_data_indiv_lh.rds"),
  area      = list(dir = "ideas_prep_area", rds = "brain_data_indiv_lh.rds"),
  vol       = list(dir = "ideas_prep_vol",  rds = "brain_data_indiv_lh.rds"),
  sulc      = list(dir = "ideas_prep_sulc", rds = "brain_data_indiv_lh.rds")
)
data(perm_id_dk_lh_5000); perm_id <- perm_id_dk_lh_5000; n_perm <- 5000
primary <- c("HS", "nonHS")
min_level_n <- 10

make_pred_list <- function(m, ids) {
  p <- list()
  p$Number_ASMs    <- data.frame(pred = as.numeric(m$Number_ASMs))
  p$age_onset      <- data.frame(pred = as.numeric(m$age_onset))
  p$SE             <- data.frame(pred = factor(m$n_SE, levels = c("0", "1")))
  p$ILAE_Y1        <- data.frame(pred = factor(ifelse(is.na(m$ILAE_Year1), NA,
                        ifelse(m$ILAE_Year1 == 1, "seizure_free", "not_free"))))
  p$FCDII_vs_other <- data.frame(pred = factor(ifelse(ids %in% fcdii_id, "FCDII", "other")))
  p
}

results <- list(); t_start <- Sys.time()
for (sname in names(strata)) {
  smd <- md[!is.na(md$path4) & md$path4 %in% strata[[sname]], ]
  for (mn in names(mods)) {
    indiv <- readRDS(file.path(ws, mods[[mn]]$dir, mods[[mn]]$rds))
    indiv <- indiv[rownames(gene_data), , drop = FALSE]
    ids_ok <- intersect(smd$ID, colnames(indiv))
    ids_ok <- ids_ok[complete.cases(t(indiv[, ids_ok, drop = FALSE]))]   # drop sulc NA case
    if (length(ids_ok) < 30) { prog(sprintf("SKIP %s x %s: n=%d", sname, mn, length(ids_ok))); next }
    msub <- md[match(ids_ok, md$ID), ]
    Yobs <- fast_scores(indiv, ids_ok)
    plist <- make_pred_list(msub, ids_ok)
    ## keep valid predictors only
    valid <- list()
    for (pn in names(plist)) {
      pd <- plist[[pn]]
      dat <- cbind(pd, data.frame(age_scan = msub$age_scan, sex = factor(msub$sex)))
      ok <- complete.cases(dat)
      lv <- table(droplevels(factor(as.character(dat$pred[ok]))))
      const <- length(unique(na.omit(dat$pred[ok]))) < 2
      if (const || (is.factor(dat$pred) && any(lv < min_level_n)) || sum(ok) < 30) next
      valid[[pn]] <- list(dat = dat, ok = ok)
    }
    if (!length(valid)) { prog(sprintf("SKIP %s x %s: no valid predictor", sname, mn)); next }
    Xs <- lapply(valid, function(v) {
      X <- model.matrix(~ pred + age_scan + sex, data = v$dat[v$ok, , drop = FALSE])
      list(X = X, rows = which(v$ok), term = grep("^pred", colnames(X), value = TRUE)[1])
    })
    t_obs <- lapply(Xs, function(v) {
      tt <- lmfit_t(v$X, Yobs[v$rows, , drop = FALSE])
      list(t = as.numeric(tt$t[v$term, ]), beta = as.numeric(tt$beta[v$term, ]),
           df = tt$df_r)
    })
    cnt <- setNames(replicate(length(Xs), integer(nset), simplify = FALSE), names(Xs))
    t0 <- Sys.time()
    Ssub <- scale(indiv[, ids_ok, drop = FALSE])     # 34 x n (column = subject)
    nS <- ncol(Ssub)
    for (i in seq_len(n_perm)) {
      GpN <- atanh(pmin(pmax(crossprod(ZgN, Ssub[perm_id[, i], , drop = FALSE]) / (nS - 1),
                             -0.999999), 0.999999))
      Sp <- t(as.matrix((Mn %*% GpN) / MsumN))       # n x 86
      for (pn in names(Xs)) {
        tp <- lmfit_t(Xs[[pn]]$X, Sp[Xs[[pn]]$rows, , drop = FALSE])$t[Xs[[pn]]$term, ]
        cnt[[pn]] <- cnt[[pn]] + (abs(tp) >= abs(t_obs[[pn]]$t))
      }
      if (i %% 500 == 0) prog(sprintf("RUN %s x %s: %d/%d perms (%.0f%%)",
                                      sname, mn, i, n_perm, 100 * i / n_perm))
    }
    for (pn in names(Xs)) {
      p_emp <- (cnt[[pn]] + 1) / (n_perm + 1)
      out <- data.frame(gene_set = names(targeted), description = tdesc,
                        stratum = sname, modality = mn, pred = pn, N = length(Xs[[pn]]$rows),
                        beta = t_obs[[pn]]$beta, t_obs = t_obs[[pn]]$t, p_emp = p_emp,
                        p_emp_fdr = p.adjust(p_emp, "BH"),
                        p_param = 2 * pt(-abs(t_obs[[pn]]$t), df = t_obs[[pn]]$df))
      results[[length(results) + 1]] <- out
    }
    prog(sprintf("DONE %s x %s: %d preds | %s | elapsed %.1f min",
                 sname, mn, length(Xs), paste(names(Xs), collapse = "+"),
                 as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  ## incremental save per stratum
  if (any(sapply(results, function(o) unique(o$stratum)) == sname)) {
    fwrite(rbindlist(results[sapply(results, function(o) unique(o$stratum)[1]) == sname]),
           file.path(outdir, sprintf("strat_%s_all.csv", sname)))
  }
}
all_res <- rbindlist(results)
fwrite(all_res, file.path(outdir, "stratified_all_results.csv"))

## ---- summary ----
lines <- c("================================================================",
           "Pathology-spectrum stratified individual-level analysis (spin 5000)",
           sprintf("Generated: %s | total runtime %.1f min", format(Sys.time()),
                   as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
           "Primary strata: HS vs non-HS (FDR within cell). Small strata: exploratory (power-limited).",
           "================================================================", "")
for (sname in names(strata)) {
  d <- all_res[stratum == sname]
  if (!nrow(d)) next
  tag <- if (sname %in% primary) "PRIMARY" else "exploratory"
  lines <- c(lines, sprintf("---- %s (%s) ----", sname, tag))
  for (mn in unique(d$modality)) for (pn in unique(d$pred[d$modality == mn])) {
    dd <- d[modality == mn & pred == pn]
    hitf <- dd[p_emp_fdr < 0.05]; hitn <- dd[p_emp < 0.05]
    lines <- c(lines, sprintf("  %s x %s (N=%s): FDR<0.05 %d/%d | nominal<0.05 %d | min p_emp=%.4g | min q=%.4g",
              mn, pn, dd$N[1], nrow(hitf), nrow(dd), nrow(hitn), min(dd$p_emp), min(dd$p_emp_fdr)))
    if (nrow(hitf)) lines <- c(lines, paste0("      FDR hits: ",
              paste(sprintf("%s(b=%.3f,q=%.3f)", hitf$gene_set, hitf$beta, hitf$p_emp_fdr), collapse = "; ")))
  }
  lines <- c(lines, "")
}
writeLines(lines, file.path(outdir, "stratified_summary.txt"))
prog("ALL STRATA COMPLETE")
