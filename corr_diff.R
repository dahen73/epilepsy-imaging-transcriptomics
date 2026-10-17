suppressPackageStartupMessages({library(BrainEnrich); library(data.table)})
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn  <- as.matrix(gene_data)
indiv <- readRDS(file.path(ws, "ideas_prep/brain_data_indiv_lh.rds"))
indiv <- indiv[rownames(gene_data), , drop = FALSE]
sub10 <- 1:10
G_pkg  <- corr_brain_gene(gene_data = gene_data, brain_data = indiv[, sub10], method = "pearson")
Zg <- scale(gn); Zb <- scale(indiv[, sub10])
G_fast <- atanh(crossprod(Zg, Zb) / (nrow(Zb) - 1))
con <- file(file.path(ws, "corr_diff.txt"), open = "wt")
writeLines(c(
  sprintf("dim G_pkg: %d x %d ; dim G_fast: %d x %d", nrow(G_pkg), ncol(G_pkg), nrow(G_fast), ncol(G_fast)),
  sprintf("rownames equal: %s", identical(rownames(G_pkg), rownames(G_fast))),
  sprintf("colnames equal: %s", identical(colnames(G_pkg), colnames(G_fast))),
  sprintf("max abs diff: %g", max(abs(G_pkg - G_fast), na.rm = TRUE)),
  sprintf("G_pkg[1:3,1] : %s", paste(sprintf("%.10f", G_pkg[1:3,1]), collapse=" | ")),
  sprintf("G_fast[1:3,1]: %s", paste(sprintf("%.10f", G_fast[1:3,1]), collapse=" | ")),
  sprintf("class gene_data: %s ; indiv: %s", class(gene_data)[1], class(indiv)[1]),
  sprintf("any NA in G_pkg: %d ; G_fast: %d", sum(is.na(G_pkg)), sum(is.na(G_fast)))
), con); close(con)
