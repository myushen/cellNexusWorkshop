# Versions from sessionInfo.txt: R 4.6.1, Bioconductor 3.24.

install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

cran_versions <- c(
  BiocManager = "1.30.27",
  remotes = "2.5.0",
  dplyr = "1.2.1",
  tidyr = "1.3.2",
  tibble = "3.3.1",
  ggplot2 = "4.0.3",
  patchwork = "1.3.2",
  stringr = "1.6.0",
  Matrix = "1.7-6",
  harmony = "2.0.5",
  Seurat = "5.5.1",
  dbplyr = "2.6.0",
  duckdb = "1.5.5",
  DBI = "1.3.0",
  arrow = "25.0.1"
)

for (pkg in names(cran_versions)) {
  remotes::install_version(
    pkg,
    version = cran_versions[[pkg]],
    upgrade = "never"
  )
}

# cmdstanr 0.9.0, CmdStan 2.37.0
remotes::install_github(
  "stan-dev/cmdstanr",
  ref = "v0.9.0",
  upgrade = "never"
)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan(version = "2.37.0")

# cellxgene.census 1.16.1 is git tag v1.18.0
remotes::install_github(
  "chanzuckerberg/cellxgene-census",
  ref = "v1.18.0",
  subdir = "api/r/cellxgene.census",
  upgrade = "never",
  repos = c(
    "https://tiledb-inc.r-universe.dev",
    "https://chanzuckerberg.r-universe.dev",
    BiocManager::repositories()
  )
)

BiocManager::install(
  c(
    "SummarizedExperiment", "SingleCellExperiment",
    "scater", "scran", "scuttle", "scrapper", "bluster",
    "DESeq2", "edgeR", "limma",
    "org.Hs.eg.db", "AnnotationDbi", "singscore",
    "tidySingleCellExperiment", "tidySummarizedExperiment",
    "tidybulk", "tidyseurat", "tidyomics",
    "HDF5Array"
  ),
  version = "3.24",
  ask = FALSE,
  update = FALSE
)

bioc_repos <- BiocManager::repositories()

remotes::install_github(
  "MangiolaLaboratory/cellNexus",
  upgrade = "never",
  repos = bioc_repos
)

remotes::install_github(
  "MangiolaLaboratory/sccomp",
  upgrade = "never",
  repos = bioc_repos
)

remotes::install_github(
  "myushen/cellNexusWorkshop",
  upgrade = "never",
  repos = bioc_repos
)
