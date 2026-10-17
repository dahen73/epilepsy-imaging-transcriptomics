# 02_smoke_test.R
# 10 例患者小样本冒烟测试：验证 brainscore.lm_test 全流程
# null_model = "none"（官方说明即用于 quick testing），DGN(DisGeNet) 基因集
suppressPackageStartupMessages({
  library(BrainEnrich)
  library(data.table)
})

out <- "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/ideas_prep"

## 1. 输入数据
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
brain_all <- readRDS(file.path(out, "brain_data_indiv_lh.rds"))
stopifnot(all(rownames(brain_all) %in% rownames(gene_data)))
brain_all <- brain_all[rownames(gene_data), , drop = FALSE]  # 行序对齐 gene_data

n_sub <- 10
brain10 <- brain_all[, 1:n_sub, drop = FALSE]
cat("[1] 冒烟样本:", n_sub, "人; ID:",
    paste(colnames(brain10), collapse = ", "), "\n")

## 2. 表型对齐
md <- readRDS(file.path(out, "phenotypes.rds"))
md$ID <- as.character(md$ID)
idx <- match(colnames(brain10), md$ID)
stopifnot(!any(is.na(idx)))
sub <- md[idx, ]

## 预测变量: Number_ASMs(连续); 协变量: age_scan + sex
## 小样本中如有 NA / 单一水平, 自动降级
pred_df <- data.frame(Number_ASMs = as.numeric(sub$Number_ASMs))
cov_df  <- data.frame(age_scan = sub$age_scan, sex = factor(sub$sex))
use_cov <- c()
if (!any(is.na(cov_df$age_scan))) use_cov <- c(use_cov, "age_scan")
if (length(levels(droplevels(cov_df$sex))) == 2) use_cov <- c(use_cov, "sex")
cov_df <- cov_df[, use_cov, drop = FALSE]
ok <- complete.cases(pred_df)
pred_df <- pred_df[ok, , drop = FALSE]; cov_df <- cov_df[ok, , drop = FALSE]
brain10 <- brain10[, ok, drop = FALSE]
cat("[2] 分析 N =", sum(ok), "; 预测变量=Number_ASMs; 协变量:",
    paste(use_cov, collapse = "+"), "\n")

## 3. 基因集 (DGN = DisGeNet 疾病基因集)
annoData <- get_annoData("DGN")

## 4. brainscore.lm_test (null_model="none" 仅参数化结果, 用于冒烟)
res <- brainscore.lm_test(
  pred_df     = pred_df,
  cov_df      = cov_df,
  brain_data  = brain10,
  gene_data   = gene_data,
  annoData    = annoData,
  cor_method  = "pearson",
  aggre_method = "mean",
  null_model  = "none",
  threshold_type = "none",
  gsea_obj    = FALSE,
  pvalueCutoff = 1
)

## 5. 输出
res_df <- as.data.frame(res)
cat("\n[3] 结果行数:", nrow(res_df), "\n")
fwrite(res_df, "C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58/smoke_test_results.csv")
cat("\n=== 按 p 值排序的前 15 行 ===\n")
print(head(res_df[order(res_df$p.val), ], 15))
cat("\n=== SMOKE_TEST_OK ===\n")
