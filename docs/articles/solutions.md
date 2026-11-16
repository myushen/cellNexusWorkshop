# Solutions

Source:
[`vignettes/solutions.Rmd`](https://github.com/myushen/cellNexusWorkshop/blob/HEAD/vignettes/solutions.Rmd)

The same answers are hidden under each exercise. This page collects them
in one place.

Run these in the session that produced `obs_lung`, `sce_qc`,
`sccomp_res`, and `epi_cluster_screen`. Counts depend on the atlas
release and on the cell cap, so the answers are the tables these
expressions print.

## Exercise 1

Age span of donors with European imputed ethnicity and healthy tissue.
Healthy tissue is `disease == "normal"`. Count each donor once.

``` r
metadata_respiratory |>
  dplyr::filter(
    imputed_ethnicity == "European",
    disease == "normal"
  ) |>
  dplyr::distinct(donor_id, age_days) |>
  dplyr::summarise(
    n_donors = dplyr::n(),
    min_age_years = min(age_days) / 365.25,
    max_age_years = max(age_days) / 365.25
  )
```

`distinct(donor_id)` keeps one row per donor, so a donor with many cells
does not pull the span. The span is `max_age_years - min_age_years`.

## Exercise 2

``` r
obs_lung |>
  dplyr::distinct(sample_id, age_bin) |>
  dplyr::count(age_bin)
```

Both bins should contain donors. A bin with a handful of samples cannot
support the age term once `dataset_id` is in the model.

## Exercise 3

Cross-tabulate distinct samples by `dataset_id` and `age_binary`, then
keep the studies with samples in both age classes.

``` r
age_by_study <- obs_lung |>
  dplyr::distinct(sample_id, dataset_id, age_binary) |>
  dplyr::mutate(age_class = ifelse(age_binary, "older", "younger")) |>
  dplyr::count(dataset_id, age_class) |>
  tidyr::pivot_wider(names_from = age_class, values_from = n, values_fill = 0)

age_by_study

age_by_study |>
  dplyr::filter(younger > 0, older > 0)
```

## Exercise 4

``` r
sce_qc |>
  ggplot2::ggplot(ggplot2::aes(UMAP1, UMAP2, colour = dataset_id)) +
  ggplot2::geom_point(size = 0.3, alpha = 0.6) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw()
```

## Exercise 5

``` r
sce_qc |>
  dplyr::count(cell_type_unified_ensemble, sort = TRUE)

setdiff(immune_cell_types, unique(sce_qc$cell_type_unified_ensemble))
```

## Exercise 6

``` r
sccomp_res |>
  dplyr::filter(factor == "age_binary", c_FDR < 0.05) |>
  dplyr::select(cell_type_unified_ensemble, dplyr::starts_with("c_"), c_FDR)
```

Read the sign of the effect column together with the boxplot. The
contrast is donors over 50 against the reference level `"FALSE"` (50 and
under), because `age_binary` was stored as character with `"FALSE"`
first.

## Exercise 7

``` r
epi_cluster_de_top
epi_cluster_screen
```

The chosen row is the first row of the screen after eligible clusters
are sorted by `n_datasets_spanning_age`, then by the number of
pseudobulk profiles. Quote `n_datasets_spanning_age` for that row. A
value of 1 means the age term rests on a single study.

## Exercise 8

`c("age_binary", "old", "young")` tests old against young. Positive
`log2FoldChange` means higher abundance in donors over 50. The boxplot
for a gene in `signature_up` should sit higher in the `50-79` bin. If it
sits lower, the factor levels and the contrast disagree.
