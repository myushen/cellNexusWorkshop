# Add age variables and drop sparse samples from a cellNexus query.

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

  A metadata table from `get_metadata()`, or a data frame with the same
  columns.

- min_cells_per_sample:

  Samples with fewer cells than this are dropped. Every cell of every
  retained sample is kept.

- tissue_groups:

  Tissue-group label. Used only when `x` has a `tissue_groups` column.
  `NULL` skips the filter.

## Value

A data frame with one row per cell.

## Details

    obs_lung <- get_metadata() |>
      keep_quality_cells() |>
      dplyr::filter(dataset_id %in% respiratory_dataset_ids()) |>
      annotate_cell_nexus(min_cells_per_sample = 50)

Adds `age_years`, `age_binary` (donors older than 50 years), `age_bin`
(`20-49` / `50-79`), and a sample-scaled `age_years_scaled_centered`.
Rewrites `dataset_id` with underscores so the identifier can be used in
a model formula. Rows are returned in a fixed order, so a seeded
subsample is reproducible. When the metadata still contains
`tissue_groups`, rows are limited to `tissue_groups`. Current cellNexus
metadata omits that column; the caller then defines tissue by
`dataset_id`.
