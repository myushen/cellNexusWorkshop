# Attach Census tissue to a cellNexus metadata query.

`get_metadata()` does not carry `tissue` or disease. This registers the
Census observation table for `census_version` on the same DuckDB
connection and left-joins `tissue`, `disease`, `sex`, and the other
Census columns used downstream, on `observation_joinid` and
`dataset_id`.

## Usage

``` r
join_census_tissue(metadata, census_version = "2024-07-01")
```

## Arguments

- metadata:

  A lazy metadata table from `get_metadata()`.

- census_version:

  Census snapshot. cellNexus uses `2024-07-01`; match this to the atlas
  you queried.

## Value

`metadata`, with the Census columns attached.

## Details

Call this after the study and cell-type filters, and before
[`annotate_cell_nexus()`](https://myushen.github.io/cellNexusWorkshop/reference/annotate_cell_nexus.md).
