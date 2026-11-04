# Add age variables and keep respiratory cells from a cellNexus query.

Intended usage, after a `get_metadata()` query:

## Usage

``` r
annotate_cell_nexus(
  x,
  min_cells_per_sample = 50L,
  tissue_groups = "respiratory system"
)
```

## Arguments

- x:

  A metadata table from `get_metadata()` after
  [`join_census_tissue()`](https://myushen.github.io/cellNexusWorkshop/reference/join_census_tissue.md),
  or a data frame with a `tissue` column.

- min_cells_per_sample:

  Samples with fewer cells than this are dropped. Every cell of every
  retained sample is kept.

- tissue_groups:

  Tissue-group label to keep. `NULL` skips the filter.

## Value

A data frame with one row per cell.

## Details

    dataset_ids <- respiratory_dataset_ids()
    obs_lung <- get_metadata() |>
      keep_quality_cells() |>
      dplyr::filter(dataset_id %in% dataset_ids) |>
      join_census_tissue() |>
      annotate_cell_nexus(min_cells_per_sample = 50)

`get_metadata()` does not carry a tissue group.
[`join_census_tissue()`](https://myushen.github.io/cellNexusWorkshop/reference/join_census_tissue.md)
attaches the CELLxGENE `tissue` label, and this function keeps the
labels that belong to the respiratory system. It then adds `age_years`,
`age_binary` (donors older than 50 years), `age_bin` (`20-49` /
`50-79`), and a sample-scaled `age_years_scaled_centered`. `dataset_id`
is rewritten with underscores so the identifier can be used in a model
formula. Rows are returned in a fixed order, so a seeded subsample is
reproducible.
