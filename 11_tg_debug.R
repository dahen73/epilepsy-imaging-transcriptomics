md <- readRDS("ideas_prep/phenotypes.rds"); md$ID <- as.character(md$ID)
indiv <- readRDS("ideas_prep/brain_data_indiv_lh.rds")
suppressPackageStartupMessages(library(BrainEnrich))
gene_data <- as.data.frame(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
gname <- rownames(gene_data)
con <- file("tg_debug.txt", "wt")
writeLines(c(
  paste("gname len:", length(gname), "| dup:", anyDuplicated(gname)),
  paste("indiv dim:", paste(dim(indiv), collapse = "x")),
  paste("md cols:", paste(names(md), collapse = ",")),
  paste("pathology class:", class(md$pathology)[1], "| vals:", paste(head(unique(as.character(md$pathology)), 12), collapse = ",")),
  paste("OP.MEMO exists:", "OP.MEMO" %in% names(md))
), con)
if ("OP.MEMO" %in% names(md)) {
  memo <- as.character(md$OP.MEMO)
  writeLines(paste("memo sample FCD:", paste(head(unique(memo[md$Pathology == "FCD" & !is.na(memo)]), 8), collapse = " ;; ")), con)
  fcdii <- md$Pathology == "FCD" & grepl("^type II($|[AB])", memo)
  writeLines(paste("fcdii count (Pathology capital):", sum(fcdii, na.rm = TRUE)), con)
}
close(con)
cat("DONE\n")
