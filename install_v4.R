# install_v4.R - 续装：BiocManager -> DOSE/enrichplot -> BrainEnrich -> 验证
lib <- "C:/Users/dahen/Rlibs/4.3"
.libPaths(c(lib, .libPaths()))
options(timeout = 3600)

# Step 0: BiocManager (官方 CRAN, 允许二进制)
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("BiocManager", lib = lib, type = "binary")
library(BiocManager)

# Step 1: DOSE / enrichplot (官方 Bioconductor)
options(repos = BiocManager::repositories())
install.packages(c("DOSE", "enrichplot"), lib = lib, ask = FALSE)

# Step 2: BrainEnrich 本地源码安装（纯 R 包）
setwd("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58")
if (!dir.exists("BrainEnrich-master")) untar("BrainEnrich-master.tar.gz")
install.packages("BrainEnrich-master", repos = NULL, type = "source", lib = lib)

# Step 3: 验证
suppressPackageStartupMessages(library(BrainEnrich))
cat("=== BrainEnrich version:", as.character(packageVersion("BrainEnrich")), "===\n")
gene <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
cat("gene_data dim:", dim(gene), "\n")
cat("前10个脑区名:\n")
print(rownames(gene)[1:10])
cat("=== INSTALL_V4_OK ===\n")
