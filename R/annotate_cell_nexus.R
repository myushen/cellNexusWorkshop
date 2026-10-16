#' Add age variables and drop sparse samples from a cellNexus query.
#'
#' Intended usage, after a `get_metadata()` query:
#'
#' ```
#' obs_lung <- get_metadata() |>
#'   keep_quality_cells() |>
#'   dplyr::filter(dataset_id %in% respiratory_dataset_ids()) |>
#'   annotate_cell_nexus(min_cells_per_sample = 50)
#' ```
#'
#' Adds `age_years`, `age_binary` (donors older than 50 years), `age_bin`
#' (`20-49` / `50-79`), and a sample-scaled `age_years_scaled_centered`.
#' Rewrites `dataset_id` with underscores so the identifier can be used in a
#' model formula. Rows are returned in a fixed order, so a seeded subsample is
#' reproducible. When the metadata still contains `tissue_groups`, rows are
#' limited to `tissue_groups`. Current cellNexus metadata omits that column;
#' the caller then defines tissue by `dataset_id`.
#'
#' @param x A metadata table from `get_metadata()`, or a data frame
#'   with the same columns.
#' @param min_cells_per_sample Samples with fewer cells than this are dropped.
#'   Every cell of every retained sample is kept.
#' @param tissue_groups Tissue-group label. Used only when `x` has a
#'   `tissue_groups` column. `NULL` skips the filter.
#' @return A data frame with one row per cell.
#' @export
annotate_cell_nexus <- function(
    x,
    min_cells_per_sample = 50L,
    tissue_groups = "respiratory system"
) {
  tissue <- tissue_groups
  min_n <- as.integer(min_cells_per_sample)

  if (!is.null(tissue)) {
    if ("tissue_groups" %in% colnames(x)) {
      x <- dplyr::filter(x, tissue_groups == !!tissue)
    } else {
      message(
        "cellNexusWorkshop: this metadata table has no tissue_groups column. ",
        "Define the tissue cohort with dataset_id before calling ",
        "annotate_cell_nexus()."
      )
    }
  }

  out <- x |>
    dplyr::collect() |>
    dplyr::mutate(
      sample_id = as.character(sample_id),
      # Atlas dataset_ids are UUIDs. Hyphens make model-matrix names
      # non-syntactic, which limma refuses. Underscores keep one identifier
      # usable in every formula; hyphens are recoverable.
      dataset_id = gsub("-", "_", as.character(dataset_id), fixed = TRUE),
      age_years = as.numeric(age_days / 365.25),
      age_binary = age_years > 50,
      age_bin = cut(
        age_years,
        breaks = c(0, 50, 80),
        right = FALSE,
        labels = c("20-49", "50-79")
      )
    ) |>
    dplyr::filter(dplyr::n() >= min_n, .by = sample_id) |>
    dplyr::arrange(dataset_id, sample_id, observation_joinid) |>
    dplyr::distinct(observation_joinid, .keep_all = TRUE)

  # Scale age once per sample so the value is constant within a sample and
  # survives aggregation to pseudobulk.
  age_scaled <- out |>
    dplyr::distinct(sample_id, age_years) |>
    dplyr::mutate(
      age_years_scaled_centered = as.numeric(scale(age_years)[, 1])
    ) |>
    dplyr::select(sample_id, age_years_scaled_centered)

  out |>
    dplyr::select(-dplyr::any_of("age_years_scaled_centered")) |>
    dplyr::left_join(age_scaled, by = "sample_id")
}
