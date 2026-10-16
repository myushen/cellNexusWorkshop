# Solutions

The same answers are hidden under each exercise. This page collects them
in one place.

Run these in the session that produced `obs_lung`, `sce_qc`,
`sccomp_res`, and `epi_cluster_screen`. Counts depend on the atlas
release and on the cell cap, so the answers are the tables these
expressions print.

## Exercise 1

Cells per harmonised label for donors older than 50 years.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellNexus``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
\
`get_metadata``(``)`` ``|>`\
`  ``keep_quality_cells``(``)`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``age_days`` ``>`` ``50`` ``*`` ``365``)`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``cell_type_unified_ensemble``, sort ``=`` ``TRUE``)`

A label with a few dozen cells across the whole atlas will have a very
wide composition interval. Drop it from the cell-type list, or accept
that its estimate is descriptive.

## Exercise 2

\
`obs_lung`` ``|>`\
`  ``dplyr``::`[`distinct`](https://dplyr.tidyverse.org/reference/distinct.html)`(``sample_id``, ``age_bin``)`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``age_bin``)`

Both bins should contain donors. A bin with a handful of samples cannot
support the age term once `dataset_id` is in the model.

## Exercise 3

\
`obs_lung`` ``|>`\
`  ``dplyr``::`[`distinct`](https://dplyr.tidyverse.org/reference/distinct.html)`(``sample_id``, ``dataset_id``, ``age_binary``)`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``dataset_id``, ``age_binary``)`` ``|>`\
`  ``tidyr``::``pivot_wider``(``names_from ``=`` ``age_binary``, values_from ``=`` ``n``, values_fill ``=`` ``0``)`

A study with a zero in either age class does not identify an age effect
on its own. The fixed `dataset_id` term uses the studies that contain
both classes.

## Exercise 4

\
`sce_qc`` ``|>`\
`  ``ggplot2``::``ggplot``(``ggplot2``::``aes``(``UMAP1``, ``UMAP2``, colour ``=`` ``dataset_id``)``)`` ``+`\
`  ``ggplot2``::``geom_point``(``size ``=`` ``0.3``, alpha ``=`` ``0.6``)`` ``+`\
`  ``ggplot2``::``coord_equal``(``)`` ``+`\
`  ``ggplot2``::``theme_bw``(``)`

Mixing of colours inside a cell-type region is the pattern you want
after Harmony. A region painted by a single study is under-correction.
The expression model still includes `dataset_id`, because the counts
were not adjusted.

## Exercise 5

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`cellNexusWorkshop`](https://myushen.github.io/cellNexusWorkshop)`)`\
\
`sce_qc`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``cell_type_unified_ensemble``, sort ``=`` ``TRUE``)`\
\
[`setdiff`](https://generics.r-lib.org/reference/setops.html)`(`[`immune_cell_types`](https://myushen.github.io/cellNexusWorkshop/reference/immune_cell_types.md)`(``)``, `[`unique`](https://rdrr.io/r/base/unique.html)`(``sce_qc``$``cell_type_unified_ensemble``)``)`

Labels missing from this subsample are absent because of the cap, not
because they are absent from `obs_lung`. Composition in the next article
is fit on `obs_lung`, which keeps every cell of every retained sample.
Fitting `sccomp` on `sce_qc` would estimate the subsample.

## Exercise 6

\
`sccomp_res`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``factor`` ``==`` ``"age_binary"``, ``c_FDR`` ``<`` ``0.05``)`` ``|>`\
`  ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``cell_type_unified_ensemble``, ``dplyr``::`[`starts_with`](https://tidyselect.r-lib.org/reference/starts_with.html)`(``"c_"``)``, ``c_FDR``)`

Read the sign of the effect column together with the boxplot. The
contrast is donors over 50 against the reference level `"FALSE"` (50 and
under), because `age_binary` was stored as character with `"FALSE"`
first.

## Exercise 7

Epithelial, pneumocyte, and secretory labels are left out so the
proportions are shares of the immune compartment. Dissociation
efficiency changes epithelial yield from study to study, and that
technical shift would move every immune proportion. Adding the
epithelial labels back answers a different question: change in share
among immune and epithelial cells together. Report the label set with
the effects.

## Exercise 8

\
`epi_cluster_de_top`\
`epi_cluster_screen``$``screen`

The chosen row is the first row of the screen after eligible clusters
are sorted by `n_datasets_spanning_age`, then by the number of
pseudobulk profiles. Quote `n_datasets_spanning_age` for that row. A
value of 1 means the age term rests on a single study.

## Exercise 9

`c("age_binary", "old", "young")` tests old against young. Positive
`log2FoldChange` means higher abundance in donors over 50. The boxplot
for a gene in `signature_up` should sit higher in the `50-79` bin. If it
sits lower, the factor levels and the contrast disagree.
