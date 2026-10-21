# Census `tissue` labels that belong to one tissue group. Only the
# respiratory system is mapped; every other tissue stays unmatched.
.tissue_grouped <- list(
  "respiratory system" = c(
    "lung", "lung parenchyma", "alveolus of lung", "bronchus",
    "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
    "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
    "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
  )
)

#' Attach Census tissue to a cellNexus metadata query.
#'
#' `get_metadata()` does not carry `tissue` or disease. This registers the
#' Census observation table for `census_version` on the same DuckDB connection
#' and left-joins `tissue`, `disease`, `sex`, and the other Census columns used
#' downstream, on `observation_joinid` and `dataset_id`.
#'
#' Call this after the study and cell-type filters, and before
#' [annotate_cell_nexus()].
#'
#' @param metadata A lazy metadata table from `get_metadata()`.
#' @param census_version Census snapshot. cellNexus uses `2024-07-01`;
#'   match this to the atlas you queried.
#' @return `metadata`, with the Census columns attached.
#' @export
join_census_tissue <- function(metadata, census_version = "2024-07-01") {
  census_metadata <- cellNexus:::get_census_metadata(census_version)
  con <- dbplyr::remote_con(metadata)
  try(duckdb::duckdb_unregister(con, "census_metadata"), silent = TRUE)
  duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)
  metadata |>
    dplyr::left_join(
      dplyr::tbl(con, "census_metadata") |>
        dplyr::select(
          observation_joinid, dataset_id, tissue,
          self_reported_ethnicity, cell_type, assay,
          disease, sex
        ),
      by = c("observation_joinid", "dataset_id")
    )
}

#' Map Census `tissue` onto `tissue_groups`.
#'
#' Limited to the respiratory system. Tissues outside that list get
#' `tissue_groups = NA`.
#' @noRd
get_tissue_grouped <- function(tbl) {
  map <- data.frame(
    tissue = unlist(.tissue_grouped, use.names = FALSE),
    tissue_groups = rep(names(.tissue_grouped), lengths(.tissue_grouped)),
    stringsAsFactors = FALSE
  )
  tbl |>
    dplyr::select(-dplyr::any_of("tissue_groups")) |>
    dplyr::left_join(map, by = "tissue", copy = TRUE)
}

#' Add age variables and keep respiratory cells from a cellNexus query.
#'
#' Intended usage, after a `get_metadata()` query:
#'
#' ```
#' obs_lung <- get_metadata() |>
#'   keep_quality_cells() |>
#'   dplyr::filter(dataset_id %in% respiratory_dataset_ids()) |>
#'   join_census_tissue() |>
#'   annotate_cell_nexus(min_cells_per_sample = 50)
#' ```
#'
#' `get_metadata()` does not carry a tissue group. [join_census_tissue()]
#' attaches the CELLxGENE `tissue` label, and this function keeps the labels
#' that belong to the respiratory system. It then adds `age_years`,
#' `age_binary` (donors older than 50 years), `age_bin` (`20-49` / `50-79`),
#' and a sample-scaled `age_years_scaled_centered`. `dataset_id` is rewritten
#' with underscores so the identifier can be used in a model formula. Rows are
#' returned in a fixed order, so a seeded subsample is reproducible.
#'
#' @param x A metadata table from `get_metadata()` after
#'   [join_census_tissue()], or a data frame with a `tissue` column.
#' @param min_cells_per_sample Samples with fewer cells than this are dropped.
#'   Every cell of every retained sample is kept.
#' @param tissue_groups Tissue-group label to keep. `NULL` skips the filter.
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
    if (!"tissue" %in% colnames(x)) {
      stop(
        "annotate_cell_nexus() needs a tissue column. ",
        "Call join_census_tissue() before annotate_cell_nexus().",
        call. = FALSE
      )
    }
    x <- x |>
      get_tissue_grouped() |>
      dplyr::filter(tissue_groups == !!tissue)
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
