indiv <- readRDS("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/ideas_prep/brain_data_indiv_lh.rds")
n_bad <- sum(colSums(is.na(indiv)) > 0)
n_cell <- sum(is.na(indiv))
con <- file("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/na_check.txt", open = "wt")
writeLines(c(sprintf("subjects with any NA: %d / %d", n_bad, ncol(indiv)),
             sprintf("NA cells total: %d / %d", n_cell, length(indiv)),
             sprintf("clean subjects: %d", ncol(indiv) - n_bad)), con); close(con)
