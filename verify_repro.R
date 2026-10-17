suppressPackageStartupMessages({library(BrainEnrich); library(data.table); library(Matrix)})
ws <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58"
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gn <- as.matrix(gene_data); Zg <- scale(gn)
indiv <- readRDS(file.path(ws, "ideas_prep/brain_data_indiv_lh.rds"))
indiv <- indiv[rownames(gene_data), , drop = FALSE]

anno_v <- get_annoData("DGN")
gs_v   <- get_geneSetList(anno_v)
gs_v   <- filter_geneSetList(rownames(gn), gs_v, minGSSize = 10, maxGSSize = 200)
cat("filter with rownames(gn) returns NULL:", is.null(gs_v), "\n")
gs_ok  <- filter_geneSetList(colnames(gn), gs_v, minGSSize = 10, maxGSSize = 200)
cat("filter with colnames(gn) sets:", length(gs_ok), "\n")

G_pkg  <- corr_brain_gene(gene_data = gene_data, brain_data = indiv[, 1:10], method = "pearson")
Zb <- scale(indiv[, 1:10])
G_fast <- atanh(crossprod(Zg, Zb) / (nrow(Zb) - 1))
ae <- all.equal(unname(G_pkg), unname(G_fast), tolerance = 1e-8)
cat("all.equal result:", isTRUE(ae), "\n")
if (!isTRUE(ae)) { cat("diff detail:\n"); print(ae) }
