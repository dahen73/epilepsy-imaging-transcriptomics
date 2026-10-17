d <- read.csv("ideas_prep/phenotypes.csv")
m <- d$OP.MEMO[!is.na(d$OP.MEMO) & d$Pathology == "FCD"]
con <- file("fcd_memo.txt", "wt")
writeLines(paste(m, collapse = " | "), con)
# 同时检查所有病理的 MEMO 中含 II / type 的比例
memo_all <- d$OP.MEMO[!is.na(d$OP.MEMO)]
writeLines(paste0("n_memo_nonNA=", length(memo_all), " / nrow=", nrow(d)), con)
close(con)
cat("DONE\n")
