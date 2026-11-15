# Respiratory ageing workflow
#
# BiocAsia workshop script. Run top to bottom in one R session.
# Articles that explain each block:
#   vignette query_the_atlas
#   vignette explore_and_visualise
#   vignette compositional_analysis
#   vignette gene_expression

library(cellNexus)
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
library(tidyseurat)
library(Seurat)
library(sccomp)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(singscore)
library(DESeq2)

# ---- 1. Query the atlas -----------------------------------------------------

metadata <- get_metadata() |>
  keep_quality_cells()

dataset_ids <- c(
  "9f222629-9e39-47d0-b83f-e08d610c7479",
  "093d3bfe-6f0f-4ac0-a7a1-829f94d0a49f",
  "8c42cfd0-0b0a-46d5-910c-fc833d83c45e",
  "53d208b0-2cfd-4366-9866-c3c6114081bc",
  "3de0ad6d-4378-4f62-b37b-ec0b75a50d94",
  "d8da613f-e681-4c69-b463-e94f5e66847f",
  "e04daea4-4412-45b5-989e-76a9be070a89",
  "2adb1f8a-a6b1-4909-8ee8-484814e2d4bf",
  "2672b679-8048-4f5e-9786-f1b196ccfd08"
)

# Built by inst/scripts/census_respiratory_parquet.R. cellxgene.census
# does not run on Windows, so the workshop reads this file.
census_path <- system.file(
  "extdata",
  "census_respiratory_2024-07-01.parquet",
  package = "cellNexusWorkshop"
)
if (!nzchar(census_path)) {
  census_path <- "inst/extdata/census_respiratory_2024-07-01.parquet"
}
census_metadata <- arrow::read_parquet(census_path)

metadata_respiratory <- metadata |>
  dplyr::filter(dataset_id %in% dataset_ids) |>
  dplyr::left_join(
    census_metadata |>
      dplyr::select(
        observation_joinid, dataset_id, tissue,
        self_reported_ethnicity, cell_type, assay,
        disease, sex
      ),
    by = c("observation_joinid", "dataset_id"),
    copy = TRUE
  )

respiratory_tissues <- c(
  "lung", "lung parenchyma", "alveolus of lung", "bronchus",
  "respiratory airway", "pleura", "pleural effusion",
  "middle lobe of right lung", "upper lobe of left lung",
  "lower lobe of left lung", "upper lobe of right lung",
  "lower lobe of right lung", "lingula of left lung",
  "right lung", "left lung"
)

immune_cell_types <- c(
  "cd16 mono", "nk", "mait", "cd14 mono", "cdc", "b memory", "pdc",
  "nkt", "cd4 th1/th17 em", "plasma", "cd4 th17 em", "b naive", "macrophage",
  "cd4 th2 em", "cd8 tcm", "cd8 naive", "cd8 tem", "tgd", "cd4 fh em",
  "cd4 th1 em", "cd4 tem", "cd4 naive", "mast", "treg", "cd4 tcm", "granulocyte"
)
epithelial_cell_types <- c("epithelial", "pneumocyte", "secretory")
cell_types_query <- c(immune_cell_types, epithelial_cell_types)

obs_lung <- metadata_respiratory |>
  dplyr::filter(
    cell_type_unified_ensemble %in% c(immune_cell_types, epithelial_cell_types),
    tissue %in% respiratory_tissues,
    age_days >= 365,
    age_days < 80 * 365
  ) |>
  dplyr::collect() |>
  dplyr::mutate(
    sample_id = as.character(sample_id),
    dataset_id = gsub("-", "_", as.character(dataset_id), fixed = TRUE),
    age_years = as.numeric(age_days / 365.25),
    age_binary = age_years > 50,
    age_bin = cut(
      age_years,
      breaks = c(0, 50, 80),
      right = FALSE,
      labels = c("20-49", "50-79")
    )
  )

obs_lung <- obs_lung |>
  dplyr::filter(dplyr::n() >= 50, .by = sample_id) |>
  dplyr::arrange(dataset_id, sample_id, observation_joinid) |>
  dplyr::distinct(observation_joinid, .keep_all = TRUE)

age_scaled <- obs_lung |>
  dplyr::distinct(sample_id, age_years) |>
  dplyr::mutate(
    age_years_scaled_centered = as.numeric(scale(age_years)[, 1])
  ) |>
  dplyr::select(sample_id, age_years_scaled_centered)

obs_lung <- obs_lung |>
  dplyr::left_join(age_scaled, by = "sample_id") |>
  dplyr::select(cell_id, dataset_id, sample_id, donor_id, atlas_id,
                dplyr::matches("age|file_id|ethnicity|cell_type|disease|sex|tissue|assay") )

if (file.exists("DESCRIPTION") && any(grepl("^Package: cellNexusWorkshop", readLines("DESCRIPTION", n = 5)))) {
  dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
  saveRDS(obs_lung, "inst/extdata/obs_lung.rds", compress = "xz")
}

# ---- 2. Explore and visualise -----------------------------------------------

max_cells_qc <- 12000L
cells_per_sample_qc <- 100L
n_samples_qc <- 120L

set.seed(42)
samples_qc <- dplyr::distinct(obs_lung, sample_id)
samples_qc <- dplyr::slice_sample(samples_qc, n = min(n_samples_qc, nrow(samples_qc)))

sampled_qc <- obs_lung |>
  dplyr::semi_join(samples_qc, by = "sample_id") |>
  dplyr::group_by(sample_id) |>
  dplyr::slice(sample.int(dplyr::n(), min(cells_per_sample_qc, dplyr::n()))) |>
  dplyr::ungroup()

# This takes ~5 minutes
sce_qc <- sampled_qc |>
  get_single_cell_experiment(assays = "counts")

sce_qc <- DelayedArray::realize(sce_qc)
saveRDS(sce_qc, "inst/extdata/sce_qc.rds", compress = "xz")

sce_qc

# This takes ~5 minutes
sce_qc <- logNormCounts(sce_qc)

hvgs_qc <- sce_qc |>
  modelGeneVar(block = sce_qc$sample_id) |>
  getTopHVGs(n = 2000)

sce_qc <- sce_qc |>
  runPCA(ncomponents = 30, subset_row = hvgs_qc) |>
  RunHarmony(
    group.by.vars = "dataset_id",
    dims.use = 1:30,
    verbose = FALSE,
    reduction.save = "HARMONY"
  ) |>
  runUMAP(dimred = "HARMONY", n_neighbors = 15)

cd45_id <- "ENSG00000081237"
sce_cd45 <- sampled_qc |>
  get_single_cell_experiment(assays = "counts", features = cd45_id)

sce_qc$cd45 <- NA_real_
i <- match(colnames(sce_qc), colnames(sce_cd45))
sce_qc$cd45[!is.na(i)] <- log1p(
  as.numeric(SummarizedExperiment::assay(sce_cd45, "counts")[1, i[!is.na(i)]])
)
rm(sce_cd45)

saveRDS(sce_qc, "inst/extdata/sce_qc_batch_corrected.rds", compress = "xz")

sce_qc |>
  ggplot(aes(UMAP1, UMAP2, colour = cd45)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_colour_viridis_c(option = "magma") +
  coord_equal() +
  theme_bw() +
  labs(colour = "CD45 (PTPRC)")

# Demonstration of get seurat object from one study
one_study <- obs_lung |>
  dplyr::filter(dataset_id == dplyr::first(dataset_id)) |>
  dplyr::slice_head(n = 100)

seurat_demo <- one_study |>
  get_seurat(assays = "counts")

# ---- 3. Compositional analysis ----------------------------------------------

obs_composition <- obs_lung |>
  dplyr::filter(cell_type_unified_ensemble %in% immune_cell_types) |>
  dplyr::mutate(age_binary = as.character(age_binary))

sccomp_res <- obs_composition |>
  sccomp_estimate(
    formula_composition = ~ age_binary + dataset_id,
    sample = "sample_id",
    cell_group = "cell_type_unified_ensemble",
    bimodal_mean_variability_association = TRUE,
    inference_method = "hmc",
    output_directory = "/vast/scratch/users/shen.m/sccomp_draws_files",
    cache_stan_model = "/vast/scratch/users/shen.m/"
  ) |>
  sccomp_remove_outliers(
    output_directory = "/vast/scratch/users/shen.m/sccomp_draws_files",
    cache_stan_model = "/vast/scratch/users/shen.m/"
  ) |>
  sccomp_test()

saveRDS(sccomp_res, "inst/extdata/sccomp_res.rds", compress = "xz")

sccomp_boxplot(
  sccomp_res,
  factor = "age_binary",
  significance_threshold = 0.05,
  remove_unwanted_effects = TRUE,
  cache_stan_model = "/vast/scratch/users/shen.m/"
)

# ---- 4. Gene expression -----------------------------------------------------

## Subset epithelial cells
epi_cell_cap <- 12000L
min_cells_pb <- 10L

set.seed(42)
epi_cells <- obs_lung |>
  dplyr::filter(cell_type_unified_ensemble == "epithelial")
epi_cells <- dplyr::slice_sample(
  epi_cells,
  n = min(epi_cell_cap, nrow(epi_cells))
)

sce_epi <- epi_cells |>
  get_single_cell_experiment(assays = "counts") |>
  logNormCounts()

hvgs <- sce_epi |>
  modelGeneVar() |>
  getTopHVGs(n = 2000)

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
  dplyr::filter(ncells >= min_cells_pb) |>
  dplyr::mutate(log_cell_count = as.numeric(scale(log1p(ncells))))

de_design <- ~ age_binary + log_cell_count + dataset_id

min_cells_per_class <- 150L
min_samples_per_class <- 5L
min_cell_fraction <- 0.25
min_datasets_spanning_age <- 2L

epi_cluster_screen <- as.data.frame(SummarizedExperiment::colData(pb_se)) |>
  dplyr::mutate(
    epi_cluster = as.character(epi_cluster),
    is_old = as.logical(age_binary)
  ) |>
  dplyr::summarise(
    n_pseudobulk = dplyr::n(),
    n_cells = sum(ncells),
    n_cells_young = sum(ncells[!is_old]),
    n_cells_old = sum(ncells[is_old]),
    n_samples_young = dplyr::n_distinct(sample_id[!is_old]),
    n_samples_old = dplyr::n_distinct(sample_id[is_old]),
    n_datasets = dplyr::n_distinct(dataset_id),
    n_datasets_spanning_age = sum(
      tapply(is_old, factor(as.character(dataset_id)), dplyr::n_distinct) == 2L,
      na.rm = TRUE
    ),
    .by = epi_cluster
  )

epi_cluster_screen$full_rank <- vapply(
  epi_cluster_screen$epi_cluster,
  function(cl) {
    se <- pb_se[, pb_se$epi_cluster == cl]
    se <- se[, Matrix::colSums(SummarizedExperiment::assay(se, "counts")) > 0]
    cd <- droplevels(as.data.frame(SummarizedExperiment::colData(se)))
    mm <- try(stats::model.matrix(de_design, data = cd), silent = TRUE)
    !inherits(mm, "try-error") && qr(mm)$rank == ncol(mm)
  },
  logical(1)
)

epi_cluster_screen <- epi_cluster_screen |>
  dplyr::mutate(
    frac_young = n_cells_young / pmax(n_cells, 1),
    frac_old = n_cells_old / pmax(n_cells, 1),
    eligible =
      n_cells_young >= min_cells_per_class &
      n_cells_old >= min_cells_per_class &
      n_samples_young >= min_samples_per_class &
      n_samples_old >= min_samples_per_class &
      frac_young >= min_cell_fraction &
      frac_old >= min_cell_fraction &
      n_datasets_spanning_age >= min_datasets_spanning_age &
      full_rank
  ) |>
  dplyr::arrange(
    dplyr::desc(eligible),
    dplyr::desc(n_datasets_spanning_age),
    dplyr::desc(n_pseudobulk),
    epi_cluster
  )

if (!any(epi_cluster_screen$eligible)) {
  stop(
    "No epithelial cluster can support this design. ",
    "Clusters are either too small, too unbalanced, or their studies each ",
    "supply a single age class."
  )
}

epi_cluster_de_top <- epi_cluster_screen$epi_cluster[[1]]

pb_se <- pb_se[, pb_se$epi_cluster == epi_cluster_de_top]
pb_se <- pb_se[, Matrix::colSums(assay(pb_se, "counts")) > 0]

assay(pb_se, "counts") <- round(assay(pb_se, "counts"))
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
  ) |>
  adjust_abundance(
    .factor_unwanted = c(log_cell_count, dataset_id),
    .factor_of_interest = age_binary,
    method = "limma_remove_batch_effect"
  )

saveRDS(pb_se, "inst/extdata/pb_se.rds", compress = "xz")

de_results <- rowData(pb_se) |>
  tibble::as_tibble(rownames = ".feature") |>
  dplyr::mutate(
    gene_symbol = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = stringr::str_replace(.feature, "\\.[0-9]+$", ""),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    ) |>
      unname() |>
      dplyr::coalesce(.feature)
  )

signature_up <- de_results |>
  dplyr::filter(pvalue < 0.05, log2FoldChange > 0) |>
  dplyr::slice_min(pvalue, n = 25)

signature_down <- de_results |>
  dplyr::filter(pvalue < 0.05, log2FoldChange < 0) |>
  dplyr::slice_min(pvalue, n = 25)

sce <- sce_epi |>
  dplyr::filter(epi_cluster == epi_cluster_de_top)


rank_mat <- singscore::rankGenes(as.matrix(assay(sce, "logcounts")))
sce$aging_signature_score <- singscore::simpleScore(
  rank_mat,
  upSet = signature_up$.feature,
  downSet = if (nrow(signature_down) > 0) signature_down$.feature else NULL
)$TotalScore

sce |>
  ggplot(aes(UMAP1, UMAP2, colour = aging_signature_score)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_colour_viridis_c(option = "magma") +
  facet_wrap(~ age_bin, ncol = 2) +
  coord_equal() +
  theme_bw() +
  labs(colour = "Ageing score", title = paste("Cluster", epi_cluster_de_top))

sce |>
  ggplot(aes(aging_signature_score, fill = age_bin, colour = age_bin)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  theme_bw() +
  labs(x = "Ageing signature score", y = "Density", fill = "Age", colour = "Age")
