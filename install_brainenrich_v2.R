options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://bioconductor.org")
options(Ncpus = 4)
options(timeout = 3600)

cat("== Step 1/4: CRAN dependencies ==\n")
cran_pkgs <- c("dplyr", "tibble", "tidyr", "purrr", "broom", "parameters",
               "pls", "data.table", "matrixStats", "pbapply", "clue", "rlang",
               "ggplot2", "ggseg", "kableExtra")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

cat("== Step 2/4: Bioconductor packages (official mirror) ==\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cat("Bioc version:", as.character(BiocManager::version()), "\n")
for (p in c("DOSE", "enrichplot")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

cat("== Step 3/4: BrainEnrich (local source) ==\n")
setwd("C:/Users/dahen/WorkBuddy/2026-10-09-14-08-58")
if (!dir.exists("BrainEnrich-master")) {
  untar("BrainEnrich-master.tar.gz", exdir = ".")
}
stopifnot(file.exists(file.path("BrainEnrich-master", "DESCRIPTION")))
remotes::install_local(path = "BrainEnrich-master",
                       dependencies = TRUE,
                       upgrade = "never",
                       build_vignettes = FALSE,
                       force = TRUE)

cat("== Step 4/4: verify ==\n")
library(BrainEnrich)
cat("BrainEnrich OK, version:", as.character(packageVersion("BrainEnrich")), "\n")
gene_data <- get_geneExp(atlas = "desikan", rdonor = "r0.4", hem = "L")
cat("gene_exp dim:", nrow(gene_data), "regions x", ncol(gene_data), "genes\n")
cat("REGION NAMES (all):\n")
print(rownames(gene_data))
cat("\nALL DONE\n")
