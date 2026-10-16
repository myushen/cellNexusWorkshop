# cellNexusWorkshop

Hands-on BiocAsia workshop. From a curated Human Cell Atlas query to donor-level composition and expression analysis, in R.

<p>
The page to follow during the session is this site. Read this page, then the articles in order.
</p>

## Instructor names and contact information

* Mengyuan Shen <shen.m at wehi.edu.au>
* Stefano Mangiola <mangiolastefano at gmail.com>

## Syllabus

Follow these articles in one R session.

1. [Query the atlas](articles/query_the_atlas.html)
2. [Explore and visualise](articles/explore_and_visualise.html)
3. [Compositional analysis](articles/compositional_analysis.html)
4. [Gene expression](articles/gene_expression.html)
5. [Solutions](articles/solutions.html)

The long-form protocol behind the case study, with the design notes for each step, is [book_chapter_Mangiola_et_al_cellNexus](https://github.com/stemangiola/book_chapter_Mangiola_et_al_cellNexus). This workshop follows the R stream of that chapter.

### Session outline

| Block | Minutes | Page |
| --- | --- | --- |
| Welcome and environment | 15 | this page |
| Query, quality control, cohort balance | 40 | Query the atlas |
| SingleCellExperiment, Harmony, UMAP | 35 | Explore and visualise |
| Immune composition with sccomp | 40 | Compositional analysis |
| Epithelial pseudobulk and DESeq2 | 40 | Gene expression |
| Close | 10 | the installed workflow script |

The signature-scoring section at the end of article 4 is the piece to skip if the room is short of time.

## Workshop package installation

The session machines should already have this stack. The commands below are for installing it afterwards, or for preparing a machine. The workshop targets R 4.5 or later. `cellNexus` is installed from GitHub. Install CmdStan before the first `sccomp` fit so the session does not stop on a compiler.

```r
install.packages("BiocManager")
install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()

BiocManager::install(c(
  "dplyr", "tidyr", "ggplot2", "patchwork", "stringr", "Matrix",
  "SummarizedExperiment", "SingleCellExperiment",
  "scater", "scran", "scuttle", "scrapper", "bluster", "harmony",
  "DESeq2", "org.Hs.eg.db", "AnnotationDbi", "singscore",
  "tidySingleCellExperiment", "tidySummarizedExperiment",
  "tidybulk", "tidyseurat", "tidyomics", "Seurat",
  "dbplyr", "duckdb"
))

remotes::install_github("MangiolaLaboratory/cellNexus")
remotes::install_github("MangiolaLaboratory/sccomp")
remotes::install_github("myushen/cellNexusWorkshop")
```

From a clone of this repository, install the workshop package with `remotes::install_local(".", dependencies = FALSE)` after the packages above.

A Docker image that runs the same install is `Dockerfile`, based on `bioconductor/bioconductor_docker:devel`.

Open the reusable script after installation:

```r
file.edit(system.file(
  "scripts", "workshop_respiratory_ageing.R",
  package = "cellNexusWorkshop"
))
```

## Workshop description

The Human Cell Atlas and CELLxGENE now host hundreds of millions of harmonised single cells spanning healthy and diseased tissues. This data availability created enormous opportunities to unveil biological mechanisms at high resolution across cell types, tissues and population niches. However, identifying a collection of atlases is only a small part of the effort. Making large-scale data analysis ready requires weeks or months of work.

This hands-on workshop introduces cellNexus, an API for R and Python that delivers quality control, annotation, and analysis layers for CELLxGENE-scale single-cell data. It allows you to conveniently create ready-to-use SingleCellExperiment and/or Seurat objects of arbitrary selections of cells. Beyond expression data, cellNexus also standardises cell type labels across studies through multi-reference consensus annotation and enriches samples with standardised metadata, including age, ethnicity, and tissue groups. This enhanced metadata supports large-scale modelling, for example of cell-type composition patterns. For large-scale compositional analyses of single-cell data, we introduce sccomp, a robust Bayesian framework for detecting differences in cell type composition and variability between conditions.

In our hands-on session, participants will learn to use cellNexus to query and download a curated subset of atlas data by tissue, cell type, or disease status, then explore and visualise the resulting data using the tidyomics ecosystem. They will subsequently perform gene expression analyses and use sccomp to perform a compositional analysis. By the end of the session, participants will have made a reusable R workflow that takes them from the Human Cell Atlas to analyses.

The live session is R. The same case study has a Python stream in the chapter repository linked above.

### Pre-requisites

* Basic familiarity with single-cell RNA-seq analysis
* Basic familiarity with the tidyverse
* An R 4.5 session with the packages above, and CmdStan installed

### Workshop participation

The workshop is a half-day hands-on session of demos, exercises, and questions. Run the article code in order, in one R session. Exercises are marked in the articles. Worked code is on the Solutions page.

The website shows the code you run. It does not execute the atlas query. Reference figures are taken from the full chapter analysis (12,000-cell cap, Hamiltonian Monte Carlo). The session defaults are a 4,000-cell cap and `sccomp` Pathfinder, so your figures will differ in detail. Set `CELLNEXUS_WORKSHOP_MAX_CELLS=12000` and `CELLNEXUS_WORKSHOP_INFERENCE=hmc` to match the chapter.

### _R_ / _Bioconductor_ packages used

* cellNexus
* sccomp
* tidySingleCellExperiment, tidySummarizedExperiment, tidybulk, tidyseurat
* scater, scran, scrapper, bluster, harmony
* DESeq2, org.Hs.eg.db, singscore

### Workshop goals and objectives

Participants leave with a script that defines a cohort from cellNexus metadata, downloads the matching counts, and fits two donor-level models: cell-type composition, and gene expression within one epithelial cluster.

#### Learning goals

* Query curated atlas metadata and turn a filter into a SingleCellExperiment or a Seurat object
* Keep the donor as the replicate for both composition and expression
* Detect composition changes with sccomp
* Test expression on sample-by-cluster pseudobulk with tidybulk

#### What you will learn

* How to filter cellNexus metadata by cell type, by a pinned set of studies, and by Census disease and tissue labels
* How tidySingleCellExperiment connects a SingleCellExperiment to dplyr and ggplot2
* A respiratory ageing case study: immune composition across nine studies, then epithelial differential expression

#### What you will not learn

* The molecular technology of single-cell sequencing
* How to train a cell-type classifier
* The Python interface, which is in the chapter repository

## Publish this site

The articles are the package vignettes. GitHub Actions (`.github/workflows/rworkflows.yml`) builds the pkgdown site and deploys the `gh-pages` branch. In the repository settings, set GitHub Pages to that branch. The address is the `url` field in `_pkgdown.yml`. Change that URL, and the GitHub links, if the repository is published under a different account.

Build the site locally with `pkgdown::build_site()`.
