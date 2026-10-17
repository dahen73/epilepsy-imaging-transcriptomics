## 37_make_tables.R - single-paper package: Tables 1-3 + supplements
suppressPackageStartupMessages({library(data.table)})
ws   <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
outdir <- file.path(ws, "tables"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
prog <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## ---------- Table 1: cohort demographics ----------
mp <- fread("F:/DJ/tables_metadata/Metadata_Release_Anon.csv")
mc <- fread("F:/DJ/tables_metadata/Metadata_Controls_Release.csv")
nrow_p <- nrow(mp); nrow_c <- nrow(mc)
sexp <- table(mp$Sex); sexc <- table(mc$Sex)
path_tab <- sort(table(mp$Pathology), decreasing = TRUE)
asm_tab <- table(mp$Number_ASMs, useNA = "ifany")
age_scan <- table(mp$Binned_Age_at_Scan)
onset <- table(mp$Binned_Onset_Age)
op_side <- table(mp$Op_Side); op_type <- table(mp$Op_Type)
se_n <- sum(mp$SE == TRUE); fus_n <- sum(mp$FUS == TRUE); fbtcs_n <- sum(mp$FBTCS == TRUE)

t1 <- data.frame(
  Variable = c("Participants, n", "Sex M/F",
               "Age at scan (binned, mode)", "Age at onset (binned, mode)",
               "Pathology, n (top)", "Number of ASMs, median [IQR]",
               "Focal to bilateral tonic-clonic seizures, n",
               "Focal impaired-awareness seizures (FUS), n",
               "Status epilepticus, n",
               "Surgery side L/R", "Surgery type (top)"),
  Patients = c(nrow_p,
               paste0(sexp["M"], "/", sexp["F"]),
               names(which.max(age_scan)), names(which.max(onset)),
               paste(names(path_tab)[1:4], path_tab[1:4], collapse = "; "),
               sprintf("%s [%s-%s]", median(mp$Number_ASMs, na.rm = TRUE),
                       quantile(mp$Number_ASMs, .25, na.rm = TRUE),
                       quantile(mp$Number_ASMs, .75, na.rm = TRUE)),
               fus_n * 0 + fbtcs_n, fus_n, se_n,
               paste0(op_side["L"], "/", op_side["R"]),
               names(which.max(op_type))),
  Controls = c(nrow_c, paste0(sexc["M"], "/", sexc["F"]),
               names(which.max(table(mc$Binned_Age_at_Scan))), "-", "-", "-", "-", "-", "-", "-", "-"))
fwrite(t1, file.path(outdir, "table1_demographics.csv"))
prog("Table 1 done"); print(t1)

## ---------- Table 2: three clusters ----------
res_vol <- fread(file.path(ws, "full_run/targeted_vol/individual_targeted_Number_ASMs.csv"))
res_are <- fread(file.path(ws, "full_run/targeted_area/individual_targeted_FCDII_vs_other.csv"))
res_voa <- fread(file.path(ws, "full_run/targeted_vol/individual_targeted_age_onset.csv"))
res_vof <- fread(file.path(ws, "full_run/targeted_vol/individual_targeted_FCDII_vs_other.csv"))
setnames(res_vol, "gene_set", "gs"); setnames(res_are, "gene_set", "gs")
setnames(res_voa, "gene_set", "gs"); setnames(res_vof, "gene_set", "gs")

fam <- function(gs, desc) {
  d <- tolower(paste(gs, desc))
  ifelse(grepl("gaba", d), "GABAergic",
  ifelse(grepl("^DGN", gs), "Epilepsy-DGN",
  ifelse(grepl("torc|mtor", d), "mTOR", "Other")))
}
c1 <- res_are[p_emp_fdr < 0.05][order(p_emp_fdr),
       .(Cluster = "1: area x FCD II", gs, description,
         beta = round(beta, 6), p_emp = round(p_emp, 5),
         q_fdr = round(p_emp_fdr, 4), Family = fam(gs, description))]
c2 <- res_vol[p_emp_fdr < 0.05][order(p_emp_fdr),
       .(Cluster = "2: vol x Number_ASMs", gs, description,
         beta = round(beta, 6), p_emp = round(p_emp, 5),
         q_fdr = round(p_emp_fdr, 4), Family = fam(gs, description))]
c3 <- res_voa[p_emp_fdr < 0.05][order(p_emp_fdr),
       .(Cluster = "3: vol x age at onset", gs, description,
         beta = round(beta, 6), p_emp = round(p_emp, 5),
         q_fdr = round(p_emp_fdr, 4), Family = fam(gs, description))]
c4 <- res_vof[p_emp_fdr < 0.05][order(p_emp_fdr),
       .(Cluster = "extra: vol x FCD II", gs, description,
         beta = round(beta, 6), p_emp = round(p_emp, 5),
         q_fdr = round(p_emp_fdr, 4), Family = fam(gs, description))]
t2 <- rbind(c1, c2, c3, c4)
fwrite(t2, file.path(outdir, "table2_clusters.csv"))
prog(sprintf("Table 2 done: %d rows (%s)", nrow(t2),
             paste(table(t2$Cluster), collapse = "/")))

## ---------- Table 3: robustness S1-S6 ----------
t3 <- data.frame(
  Test = c("S1 Resection-burden covariate", "S2 Leave-one-out jackknife",
           "S3 80% subsample x 200", "S4 High-resection-region exclusion (FCD II cluster)",
           "S5 Right-hemisphere replication (all clusters)",
           "S6 Severity adjustment (ASM cluster)"),
  Design = c("lm + per-patient cortical resection burden", "remove each patient in turn",
             "random subsamples, bootstrap", "exclude 19 high-resection regions, rescore + spin 5000",
             "full pipeline on 32 RH regions, RH spin 5000",
             "Freedman-Lane 5000, cov: age/sex/duration/SE/pathology"),
  Result = c("12/12 direction+significance retained; |beta change| median 1.5%",
             "0/4992 sign flips; P(p<0.05)=1.00",
             "direction consistency 1.00; 12/12 base beta inside 95% CI",
             "direction 18/18; r(beta)=0.913; p_emp>=0.13 (power-limited)",
             "48/49 sign retained; r_s(beta) 0.78-0.94 by cluster; no RH survivors",
             "28/28 p_emp<0.05 retained (min 0.010); no FDR survivors under reduced power"),
  Verdict = c("robust", "robust", "robust", "direction robust, power-limited",
              "effect-size concordant, RH-driven significance", "robust (nominal), power-limited"))
fwrite(t3, file.path(outdir, "table3_robustness.csv"))
prog("Table 3 done")

## ---------- Supplements ----------
fs <- c(Sys.glob(file.path(ws, "full_run/targeted",       "individual_targeted_*.csv")),
        Sys.glob(file.path(ws, "full_run/targeted_area", "individual_targeted_*.csv")),
        Sys.glob(file.path(ws, "full_run/targeted_vol",  "individual_targeted_*.csv")),
        Sys.glob(file.path(ws, "full_run/targeted_sulc", "individual_targeted_*.csv")))
sup1 <- rbindlist(lapply(fs, function(f) {
  d <- fread(f); setnames(d, "gene_set", "gs")
  d[, mod := basename(dirname(f))]
  d[, pred := sub("individual_targeted_|\\.csv$", "", basename(f))]
  d
}), use.names = TRUE, fill = TRUE)
fwrite(sup1, file.path(outdir, "suptable1_all_individual_targeted.csv"))
prog(sprintf("Supp 1: %d rows from %d files", nrow(sup1), length(fs)))

sup2 <- fread(file.path(ws, "full_run/rh/rh_vs_lh.csv"))
fwrite(sup2, file.path(outdir, "suptable2_rh_vs_lh.csv"))

sup3 <- fread(file.path(ws, "full_run/targeted_vol/region_age_onset.csv"))
fwrite(sup3, file.path(outdir, "suptable3_region_age_onset_vol.csv"))
prog(sprintf("Supp 2: %d rows | Supp 3: %d rows", nrow(sup2), nrow(sup3)))
