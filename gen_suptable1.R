suppressPackageStartupMessages(library(data.table))
tgt <- readRDS("D:/DJ/Epilepsia_work/04_full_run_results/targeted/targeted_sets.rds")
sets <- tgt$sets; desc <- tgt$desc
lib <- sub(":.*", "", names(sets))
fam <- ifelse(grepl("mTOR|TORC", desc, ignore.case = TRUE) | grepl("mTOR|TORC", names(sets), ignore.case = TRUE), "mTOR",
       ifelse(grepl("GABA", desc, ignore.case = TRUE), "GABAergic",
       ifelse(lib == "DGN", "Epilepsy-DGN", "other")))
sz <- sapply(sets, length)
dt <- data.table(gene_set = names(sets), description = desc, library = lib, family = fam, n_genes = sz)
fwrite(dt, "C:/Epilepsia_Submission/tables_new/suptable1_targeted_panel_86sets.csv")
print(table(fam)); cat("total:", nrow(dt), "\n")
