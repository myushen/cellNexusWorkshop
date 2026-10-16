#' Choose the epithelial cluster for the donor-level age model.
#'
#' k-means labels are arbitrary: they change with the seed, the cell cap, and
#' the gene selection. The cluster is chosen by a screen.
#'
#' A cluster qualifies when both age classes are well represented and the age
#' effect is separable from the study effect. If every study in a cluster
#' supplies a single age class, age is a linear combination of the study
#' indicators and no age coefficient can be separated from a study effect.
#' Among clusters that qualify, the one supported by the most studies that
#' contribute both age classes is chosen, then the one with the most pseudobulk
#' columns. Ties break on the label.
#'
#' @param pb_se A SummarizedExperiment of sample-by-cluster pseudobulks, with
#'   `epi_cluster`, `sample_id`, `dataset_id`, `age_binary`, and `ncells` in
#'   the column data, and a `counts` assay.
#' @param de_formula The design used for the expression model, for example
#'   `~ age_binary + log_cell_count + dataset_id`.
#' @param min_cells_per_class Minimum cells in each age class.
#' @param min_samples_per_class Minimum samples in each age class.
#' @param min_cell_fraction Minimum share of cells in each age class.
#' @param min_datasets_spanning_age Minimum number of studies that contribute
#'   donors of both age classes.
#' @return A list with `cluster`, the chosen label, and `screen`, the table the
#'   choice was made from.
#' @export
select_de_cluster <- function(
    pb_se,
    de_formula,
    min_cells_per_class = 150L,
    min_samples_per_class = 5L,
    min_cell_fraction = 0.25,
    min_datasets_spanning_age = 2L
) {
  screen <- as.data.frame(SummarizedExperiment::colData(pb_se)) |>
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
    ) |>
    dplyr::mutate(
      frac_young = n_cells_young / pmax(n_cells, 1),
      frac_old = n_cells_old / pmax(n_cells, 1),
      full_rank = vapply(
        epi_cluster,
        function(cl) design_is_full_rank(pb_se, cl, de_formula),
        logical(1)
      ),
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

  if (!any(screen$eligible)) {
    stop(
      "No epithelial cluster can support ", deparse1(de_formula),
      ". Clusters are either too small, too unbalanced, or their studies each ",
      "supply a single age class. Screen:\n",
      paste(utils::capture.output(print(screen)), collapse = "\n")
    )
  }

  list(cluster = screen$epi_cluster[[1]], screen = screen)
}

#' Is `de_formula` estimable on the pseudobulks of one cluster?
#' @noRd
design_is_full_rank <- function(pb_se, epi_cluster_label, de_formula) {
  se <- pb_se[, pb_se$epi_cluster == epi_cluster_label]
  se <- se[, Matrix::colSums(SummarizedExperiment::assay(se, "counts")) > 0]
  cd <- droplevels(as.data.frame(SummarizedExperiment::colData(se)))
  if (ncol(se) < 3L || dplyr::n_distinct(cd$age_binary) < 2L) return(FALSE)
  terms_used <- all.vars(de_formula)
  usable_terms <- vapply(
    terms_used,
    function(v) {
      column <- cd[[v]]
      is.numeric(column) || dplyr::n_distinct(column) >= 2L
    },
    logical(1)
  )
  if (!all(usable_terms)) return(FALSE)
  mm <- tryCatch(
    stats::model.matrix(de_formula, data = cd),
    error = function(e) NULL
  )
  if (is.null(mm)) return(FALSE)
  qr(mm)$rank == ncol(mm)
}
