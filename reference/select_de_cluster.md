# Choose the epithelial cluster for the donor-level age model.

k-means labels are arbitrary: they change with the seed, the cell cap,
and the gene selection. The cluster is chosen by a screen.

## Usage

``` r
select_de_cluster(
  pb_se,
  de_formula,
  min_cells_per_class = 150L,
  min_samples_per_class = 5L,
  min_cell_fraction = 0.25,
  min_datasets_spanning_age = 2L
)
```

## Arguments

- pb_se:

  A SummarizedExperiment of sample-by-cluster pseudobulks, with
  `epi_cluster`, `sample_id`, `dataset_id`, `age_binary`, and `ncells`
  in the column data, and a `counts` assay.

- de_formula:

  The design used for the expression model, for example
  `~ age_binary + log_cell_count + dataset_id`.

- min_cells_per_class:

  Minimum cells in each age class.

- min_samples_per_class:

  Minimum samples in each age class.

- min_cell_fraction:

  Minimum share of cells in each age class.

- min_datasets_spanning_age:

  Minimum number of studies that contribute donors of both age classes.

## Value

A list with `cluster`, the chosen label, and `screen`, the table the
choice was made from.

## Details

A cluster qualifies when both age classes are well represented and the
age effect is separable from the study effect. If every study in a
cluster supplies a single age class, age is a linear combination of the
study indicators and no age coefficient can be separated from a study
effect. Among clusters that qualify, the one supported by the most
studies that contribute both age classes is chosen, then the one with
the most pseudobulk columns. Ties break on the label.
