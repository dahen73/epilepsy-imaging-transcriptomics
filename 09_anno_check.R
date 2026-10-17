suppressPackageStartupMessages(library(BrainEnrich))
e <- get_annoData("KEGG")
con <- file("anno_struct.txt", "wt")
writeLines(paste("objects:", paste(ls(e), collapse = ",")), con)
if ("annoData" %in% ls(e)) {
  a <- e$annoData
  writeLines(paste("annoData class:", class(a), " dim:", paste(dim(a), collapse = "x")), con)
  writeLines(paste("cols:", paste(names(a), collapse = ",")), con)
  writeLines(paste("head:", paste(apply(head(a, 3), 1, paste, collapse = "|"), collapse = " ;; ")), con)
}
for (nm in setdiff(ls(e), "annoData")) {
  o <- e[[nm]]
  writeLines(sprintf("%s : %s %s", nm, class(o)[1],
    if (is.data.frame(o)) paste(dim(o), collapse = "x") else
    if (is.list(o)) paste("list len", length(o)) else
    paste("len", length(o))), con)
}
close(con)
cat("DONE\n")
