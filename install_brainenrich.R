options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(Ncpus = 4)
options(timeout = 3600)

cat("== Step 1/4: remotes & BiocManager ==\n")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")

cat("== Step 2/4: DOSE (Bioconductor) ==\n")
if (!requireNamespace("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE", ask = FALSE, update = FALSE)
}

cat("== Step 3/4: BrainEnrich (GitHub) ==\n")
remotes::install_github("zh1peng/BrainEnrich",
                        dependencies = TRUE,
                        upgrade = "never",
                        build_vignettes = FALSE,
                        force = TRUE)

cat("== Step 4/4: verify installation ==\n")
library(BrainEnrich)
cat("BrainEnrich OK, version:", as.character(packageVersion("BrainEnrich")), "\n")
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
cat("gene_exp dim (regions x genes):", dim(gene_data)[1], "x", dim(gene_data)[2], "\n")
cat("region names (first 10):\n")
print(head(rownames(gene_data), 10))
cat("column name example (genes):\n")
print(head(colnames(gene_data), 5))
cat("\nALL DONE\n")
