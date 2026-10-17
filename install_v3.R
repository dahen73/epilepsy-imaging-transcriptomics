## install_v3.R —— 修复版: Rtools43 工具链 + 新库 + 二进制优先
rt <- "C:/rtools43/usr/bin"
if (dir.exists(rt)) Sys.setenv(PATH = paste(rt, Sys.getenv("PATH"), sep = ";"))
cat("make ->", Sys.which("make")[["make"]], "\n")
cat("sh   ->", Sys.which("sh")[["sh"]], "\n")
cat("cat  ->", Sys.which("cat")[["cat"]], "\n")
cat("lib  ->", .libPaths(), sep = "\n")

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 3600)

cat("\n== Step 1/4: CRAN dependencies ==\n")
cran_pkgs <- c("dplyr", "tibble", "tidyr", "purrr", "stringr", "broom",
               "parameters", "pls", "data.table", "matrixStats", "pbapply",
               "clue", "rlang", "cli", "ggplot2", "ggseg", "kableExtra",
               "generics", "cpp11", "testthat")
miss <- cran_pkgs[!vapply(cran_pkgs, function(p)
  requireNamespace(p, quietly = TRUE), logical(1))]
cat("missing:", paste(miss, collapse = ", "), "\n")
if (length(miss)) install.packages(miss)

cat("\n== Step 2/4: Bioconductor (official) ==\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cat("Bioc:", as.character(BiocManager::version()), "\n")
for (p in c("DOSE", "enrichplot")) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = FALSE)
}

cat("\n== Step 3/4: BrainEnrich (local source) ==\n")
setwd("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58")
if (!dir.exists("BrainEnrich-master")) untar("BrainEnrich-master.tar.gz", exdir = ".")
if (file.exists("BrainEnrich_1.0.0.tar.gz")) file.remove("BrainEnrich_1.0.0.tar.gz")
system2(file.path(R.home("bin"), "R.exe"),
        c("CMD", "build", "--no-build-vignettes", "BrainEnrich-master"))
stopifnot(file.exists("BrainEnrich_1.0.0.tar.gz"))
install.packages("BrainEnrich_1.0.0.tar.gz", repos = NULL, type = "source")

cat("\n== Step 4/4: verify ==\n")
library(BrainEnrich)
cat("BrainEnrich", as.character(packageVersion("BrainEnrich")), "OK\n")
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
cat("regions:", nrow(gene_data), " genes:", ncol(gene_data), "\n")
cat("REGION NAMES:\n")
print(rownames(gene_data))
cat("\nALL DONE\n")
