suppressPackageStartupMessages(library(BrainEnrich))
d <- read.csv("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/smoke_test_results.csv")
d <- d[order(d$p.val), ]
anno <- get_annoData("DGN")
top <- head(d, 15)
desc <- get_termDescription(top$Dependent_vars, anno)
out <- data.frame(
  CUI = top$Dependent_vars,
  Disease = ifelse(is.na(desc[top$Dependent_vars]), "?", desc[top$Dependent_vars]),
  Std_Coef = round(top$Standardized_Coefficient, 3),
  t = round(top$`t_Value`, 2),
  p = signif(top$p.val, 3)
)
write.csv(out, "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/smoke_top_named.csv", row.names = FALSE)
print(out, row.names = FALSE)
