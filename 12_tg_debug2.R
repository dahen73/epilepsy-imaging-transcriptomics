suppressPackageStartupMessages(library(BrainEnrich))
md <- readRDS("ideas_prep/phenotypes.rds"); md$ID <- as.character(md$ID)
indiv <- readRDS("ideas_prep/brain_data_indiv_lh.rds")
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname <- rownames(gene_data)
idx <- match(colnames(indiv), md$ID); keep <- !is.na(idx)
indiv2 <- indiv[, keep, drop = FALSE]; mdc <- md[idx[keep], ]
con <- file("tg_debug2.txt", "wt")
writeLines(c(paste("indiv2 dim:", paste(dim(indiv2), collapse = "x")),
             paste("mdc rows:", nrow(mdc))), con)
sel <- mdc$pathology == "HS"; sel[is.na(sel)] <- FALSE
writeLines(paste("sel len:", length(sel), " TRUE:", sum(sel)), con)
cm <- colMeans(indiv2[, sel, drop = FALSE])
writeLines(c(paste("colMeans len:", length(cm)),
             paste("dev len test:", length(data.frame(dev = cm)$dev))), con)
df <- tryCatch(data.frame(dev = cm, row.names = gname), error = function(e) conditionMessage(e))
writeLines(paste("df build:", if (is.data.frame(df)) "OK" else paste("ERR:", df)), con)
close(con)
cat("DONE\n")
