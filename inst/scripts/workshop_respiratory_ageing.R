# Respiratory ageing workflow
#
# BiocAsia workshop script. Run top to bottom in one R session.
# Articles that explain each block:
#   vignette query_the_atlas
#   vignette explore_and_visualise
#   vignette compositional_analysis
#   vignette gene_expression
#
# Workshop pace (defaults):
#   CELLNEXUS_WORKSHOP_MAX_CELLS=4000
#   CELLNEXUS_WORKSHOP_INFERENCE=pathfinder
# Larger example run:
#   CELLNEXUS_WORKSHOP_MAX_CELLS=12000
#   CELLNEXUS_WORKSHOP_INFERENCE=hmc

library(cellNexus)
library(cellNexusWorkshop)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scater)
library(scuttle)
library(scran)
library(scrapper)
library(bluster)
library(harmony)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(tidySingleCellExperiment)
library(tidySummarizedExperiment)
library(tidybulk)
library(sccomp)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)

max_cells <- as.integer(Sys.getenv("CELLNEXUS_WORKSHOP_MAX_CELLS", "4000"))
inference_method <- Sys.getenv("CELLNEXUS_WORKSHOP_INFERENCE", "pathfinder")

# ---- Cohort ----------------------------------------------------------------

cell_types_query <- c(immune_cell_types(), epithelial_cell_types())

obs_lung <- get_metadata() |>
  keep_quality_cells() |>
  dplyr::filter(
    dataset_id %in% respiratory_dataset_ids(),
    cell_type_unified_ensemble %in% cell_types_query,
    age_days >= 365,
    age_days < 80 * 365
  ) |>
  annotate_cell_nexus(min_cells_per_sample = 50)

# ---- Composition ------------------------------------------------------------

obs_composition <- obs_lung |>
  dplyr::filter(cell_type_unified_ensemble %in% immune_cell_types()) |>
  dplyr::mutate(age_binary = as.character(age_binary))

sccomp_res <- obs_composition |>
  sccomp_estimate(
    formula_composition = ~ age_binary + dataset_id,
    sample = "sample_id",
    cell_group = "cell_type_unified_ensemble",
    bimodal_mean_variability_association = TRUE,
    inference_method = inference_method
  ) |>
  sccomp_remove_outliers() |>
  sccomp_test()

# ---- Epithelial integration ------------------------------------------------

set.seed(42)
epi_cells <- obs_lung |>
  dplyr::filter(cell_type_unified_ensemble == "epithelial")
epi_cells <- dplyr::slice_sample(epi_cells, n = min(max_cells, nrow(epi_cells)))

sce_epi <- epi_cells |>
  get_single_cell_experiment(assays = "counts") |>
  logNormCounts()

hvgs <- sce_epi |> modelGeneVar() |> getTopHVGs(n = 2000)

sce_epi <- sce_epi |>
  runPCA(ncomponents = 30, subset_row = hvgs) |>
  RunHarmony(
    group.by.vars = "dataset_id",
    dims.use = 1:30,
    verbose = FALSE,
    reduction.save = "HARMONY"
  ) |>
  runUMAP(dimred = "HARMONY", n_neighbors = 15)

sce_epi$epi_cluster <- factor(paste0(
  "E",
  clusterCells(
    sce_epi,
    use.dimred = "HARMONY",
    BLUSPARAM = KmeansParam(centers = 6L)
  )
))

pb_se <- scrapper::aggregateAcrossCells.se(
  sce_epi,
  factors = list(
    sample_id = sce_epi$sample_id,
    epi_cluster = sce_epi$epi_cluster
  ),
  counts.name = "ncells"
)
SummarizedExperiment::assayNames(pb_se)[
  SummarizedExperiment::assayNames(pb_se) == "sums"
] <- "counts"

pb_se <- pb_se |>
  dplyr::filter(ncells >= 10) |>
  dplyr::mutate(log_cell_count = as.numeric(scale(log1p(ncells))))

# ---- Differential expression ------------------------------------------------

de_design <- ~ age_binary + log_cell_count + dataset_id

epi_cluster_screen <- select_de_cluster(
  pb_se,
  de_design,
  min_cells_per_class = 100L,
  min_samples_per_class = 4L
)

pb_se <- pb_se[, pb_se$epi_cluster == epi_cluster_screen$cluster]
pb_se <- pb_se[, Matrix::colSums(SummarizedExperiment::assay(pb_se, "counts")) > 0]
SummarizedExperiment::assay(pb_se, "counts") <- round(
  SummarizedExperiment::assay(pb_se, "counts")
)
pb_se$sizeFactor <- NULL
pb_se$age_binary <- factor(
  ifelse(pb_se$age_binary, "old", "young"),
  levels = c("young", "old")
)

pb_se <- pb_se |>
  keep_abundant(
    formula_design = ~ age_binary,
    minimum_proportion = 0.5,
    minimum_count_per_million = 5
  ) |>
  test_differential_expression(
    de_design,
    contrasts = list(c("age_binary", "old", "young")),
    method = "DESeq2",
    omit_contrast_in_colnames = TRUE,
    test_above_log2_fold_change = 1,
    sfType = "poscounts"
  )

de_results <- SummarizedExperiment::rowData(pb_se) |>
  tibble::as_tibble(rownames = ".feature")
