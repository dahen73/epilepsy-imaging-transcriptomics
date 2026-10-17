suppressPackageStartupMessages(library(data.table))
zp <- fread("F:/DJ/zscores/thick_zscores_patients.csv", select = "ID", nrows = 5)
md <- fread("F:/DJ/tables_metadata/Metadata_Release_Anon.csv", select = "ID")
zp_all <- fread("F:/DJ/zscores/thick_zscores_patients.csv", select = "ID")
con <- file("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/id_check.txt", open = "wt")
writeLines(c(
  sprintf("z表 n=%d, 元数据 n=%d", nrow(zp_all), nrow(md)),
  sprintf("z表 ID 样例: %s | %s | %s", zp$ID[1], zp$ID[2], zp$ID[3]),
  sprintf("元数据 ID 样例: %s | %s | %s", md$ID[1], md$ID[2], md$ID[3]),
  sprintf("z表 ID 在元数据中未找到: %d", sum(!(as.character(zp_all$ID) %in% as.character(md$ID)))),
  sprintf("元数据 ID 在 z 表中未找到: %d", sum(!(as.character(md$ID) %in% as.character(zp_all$ID))))
), con); close(con)
