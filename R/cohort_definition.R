#' Dataset identifiers for the respiratory reference cohort.
#'
#' Nine curated respiratory studies from the Human Cell Atlas. Pinning the
#' identifiers keeps the workshop cohort stable when the atlas gains studies.
#' Identifiers use the hyphenated form stored by `get_metadata()`.
#' `annotate_cell_nexus()` rewrites hyphens to underscores so the same
#' identifier can be placed in a model formula.
#'
#' @return A character vector of `dataset_id` values.
#' @export
#' @examples
#' respiratory_dataset_ids()
respiratory_dataset_ids <- function() {
  c(
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
}

#' Immune labels used as the composition denominator.
#'
#' Composition is modelled inside the immune compartment. Epithelial yield in
#' dissociated lung depends on the digestion protocol, so epithelial labels are
#' kept out of this list and analysed in the expression article.
#'
#' @return A character vector of `cell_type_unified_ensemble` labels.
#' @export
#' @examples
#' immune_cell_types()
immune_cell_types <- function() {
  c(
    "cd16 mono", "nk", "mait", "cd14 mono", "cdc", "b memory", "pdc",
    "nkt", "cd4 th1/th17 em", "plasma", "cd4 th17 em", "b naive", "macrophage",
    "cd4 th2 em", "cd8 tcm", "cd8 naive", "cd8 tem", "tgd", "cd4 fh em",
    "cd4 th1 em", "cd4 tem", "cd4 naive", "mast", "treg", "cd4 tcm", "granulocyte"
  )
}

#' Epithelial labels used for the expression analysis.
#'
#' @return A character vector of `cell_type_unified_ensemble` labels.
#' @export
#' @examples
#' epithelial_cell_types()
epithelial_cell_types <- function() {
  c("epithelial", "pneumocyte", "secretory")
}
