suppressPackageStartupMessages({library(BrainEnrich); library(data.table)})
g <- rownames(get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L"))
strip <- function(x){x<-sub("^(lh|rh)_","",x); x<-sub("_(thickness|thickavg)$","",x); x}
zp <- fread("F:/DJ/zscores/thick_zscores_patients.csv", nrows = 1)
r68 <- setdiff(names(zp), "ID")
lh <- grep("^lh_", r68, value = TRUE)
writeLines(c(
  paste("gene_head5:", paste(head(g, 5), collapse = ",")),
  paste("ideas_head5:", paste(head(lh, 5), collapse = ",")),
  paste("unmatched_from_ideas:", paste(setdiff(strip(lh), strip(g)), collapse = ",")),
  paste("extra_in_gene:", paste(setdiff(strip(g), strip(lh)), collapse = ","))
), "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/align_check.txt")
