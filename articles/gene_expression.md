# 4. Gene expression

Continue in the same session. This block integrates epithelial cells,
aggregates sample-by-cluster pseudobulk, and tests donor-level age
effects with DESeq2 through tidybulk (Love et al. 2014). The last
section scores the resulting signature back onto single cells (Foroutan
et al. 2018).

The reusable script for everything you have run is installed with the
package:

\
[`file.edit`](https://rdrr.io/r/utils/file.edit.html)`(`[`system.file`](https://rdrr.io/r/base/system.file.html)`(`\
`  ``"scripts"``, ``"workshop_respiratory_ageing.R"``,`\
`  package ``=`` ``"cellNexusWorkshop"`\
`)``)`

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellNexus``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`cellNexusWorkshop`](https://myushen.github.io/cellNexusWorkshop)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ggplot2`](https://ggplot2.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scater``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scuttle``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scran``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scrapper``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`bluster`](https://bioconductor.org/packages/bluster)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`harmony`](https://github.com/immunogenomics/harmony)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``SingleCellExperiment``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`SummarizedExperiment`](https://bioconductor.org/packages/SummarizedExperiment)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`Matrix`](https://codeberg.org/jaganmn/Matrix)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``tidySingleCellExperiment``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``tidySummarizedExperiment``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``tidybulk``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``AnnotationDbi``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`org.Hs.eg.db`](https://bioconductor.org/packages/org.Hs.eg.db)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`stringr`](https://stringr.tidyverse.org)`)`\
\
`if`` ``(``!`[`exists`](https://rdrr.io/r/base/exists.html)`(``"obs_lung"``)``)`` ``{`\
`  ``obs_lung`` ``<-`` ``get_metadata``(``)`` ``|>`\
`    ``keep_quality_cells``(``)`` ``|>`\
`    ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(`\
`      ``dataset_id`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`respiratory_dataset_ids`](https://myushen.github.io/cellNexusWorkshop/reference/respiratory_dataset_ids.md)`(``)``,`\
`      ``cell_type_unified_ensemble`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`c`](https://rdrr.io/r/base/c.html)`(`[`immune_cell_types`](https://myushen.github.io/cellNexusWorkshop/reference/immune_cell_types.md)`(``)``, `[`epithelial_cell_types`](https://myushen.github.io/cellNexusWorkshop/reference/epithelial_cell_types.md)`(``)``)``,`\
`      ``age_days`` ``>=`` ``365``,`\
`      ``age_days`` ``<`` ``80`` ``*`` ``365`\
`    ``)`` ``|>`\
`    `[`annotate_cell_nexus`](https://myushen.github.io/cellNexusWorkshop/reference/annotate_cell_nexus.md)`(``min_cells_per_sample ``=`` ``50``)`\
`}`

## Integrate and cluster the epithelial compartment

The broad label `epithelial` mixes ciliated, secretory, basal, and
alveolar states. We integrate cells, cut the embedding into six
clusters, and test age inside the cluster that can support the design.

The download is capped by `CELLNEXUS_WORKSHOP_MAX_CELLS` (default 4000).
The example figures use 12000 cells. Harmony here is for clustering and
the UMAP. The counts still carry study effects, and the model below
still includes `dataset_id`.

\
`epi_cell_cap`` ``<-`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(`[`Sys.getenv`](https://rdrr.io/r/base/Sys.getenv.html)`(``"CELLNEXUS_WORKSHOP_MAX_CELLS"``, ``"4000"``)``)`\
`min_cells_pb`` ``<-`` ``10L`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`epi_cells`` ``<-`` ``obs_lung`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``cell_type_unified_ensemble`` ``==`` ``"epithelial"``)`\
`epi_cells`` ``<-`` ``dplyr``::`[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)`(`\
`  ``epi_cells``,`\
`  n ``=`` `[`min`](https://rdrr.io/r/base/Extremes.html)`(``epi_cell_cap``, `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``epi_cells``)``)`\
`)`\
\
`sce_epi`` ``<-`` ``epi_cells`` ``|>`\
`  ``get_single_cell_experiment``(``assays ``=`` ``"counts"``)`` ``|>`\
`  ``logNormCounts``(``)`\
\
`hvgs`` ``<-`` ``sce_epi`` ``|>`\
`  ``modelGeneVar``(``)`` ``|>`\
`  ``getTopHVGs``(``n ``=`` ``2000``)`\
\
`sce_epi`` ``<-`` ``sce_epi`` ``|>`\
`  ``runPCA``(``ncomponents ``=`` ``30``, subset_row ``=`` ``hvgs``)`` ``|>`\
`  ``RunHarmony``(`\
`    group.by.vars ``=`` ``"dataset_id"``,`\
`    dims.use ``=`` ``1``:``30``,`\
`    verbose ``=`` ``FALSE``,`\
`    reduction.save ``=`` ``"HARMONY"`\
`  ``)`` ``|>`\
`  ``runUMAP``(``dimred ``=`` ``"HARMONY"``, n_neighbors ``=`` ``15``)`\
\
`sce_epi``$``epi_cluster`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(`\
`  ``"E"``,`\
`  ``clusterCells``(`\
`    ``sce_epi``,`\
`    use.dimred ``=`` ``"HARMONY"``,`\
`    BLUSPARAM ``=`` ``KmeansParam``(``centers ``=`` ``6L``)`\
`  ``)`\
`)``)`

\
`sce_epi`` ``|>`\
`  ``ggplot``(``aes``(``UMAP1``, ``UMAP2``, colour ``=`` ``epi_cluster``)``)`` ``+`\
`  ``geom_point``(``size ``=`` ``0.3``, alpha ``=`` ``0.6``)`` ``+`\
`  ``coord_equal``(``)`` ``+`\
`  ``theme_bw``(``)`

![Example Harmony UMAP of epithelial cells, coloured by the six k-means
clusters.](images/epithelial-umap.png)

Example Harmony UMAP of epithelial cells, coloured by the six k-means
clusters.

k-means labels are arbitrary. `E1` in your session is not `E1` in the
example figure. The next section chooses the cluster from a screen,
which you report alongside the label.

## Aggregate sample-by-cluster pseudobulk

`aggregateAcrossCells.se()` sums raw counts within each sample and
cluster. The column `ncells` is the number of cells that entered each
pseudobulk. Profiles with fewer than 10 cells are dropped.
`log_cell_count` is a scaled depth covariate so profiles built from
different numbers of cells are comparable.

Precomputed `get_pseudobulk()` is sample by broad cell type. It is the
right object for a cell-type-level question. This article needs sample
by cluster, so the aggregation is done here.

\
`pb_se`` ``<-`` ``scrapper``::``aggregateAcrossCells.se``(`\
`  ``sce_epi``,`\
`  factors ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    sample_id ``=`` ``sce_epi``$``sample_id``,`\
`    epi_cluster ``=`` ``sce_epi``$``epi_cluster`\
`  ``)``,`\
`  counts.name ``=`` ``"ncells"`\
`)`\
\
`SummarizedExperiment``::`[`assayNames`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``)``[`\
`  ``SummarizedExperiment``::`[`assayNames`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``)`` ``==`` ``"sums"`\
`]`` ``<-`` ``"counts"`\
\
`pb_se`` ``<-`` ``pb_se`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``ncells`` ``>=`` ``min_cells_pb``)`` ``|>`\
`  ``dplyr``::`[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``log_cell_count ``=`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(`[`scale`](https://rdrr.io/r/base/scale.html)`(`[`log1p`](https://rdrr.io/r/base/Log.html)`(``ncells``)``)``)``)`\
\
`pb_se`

Confirm that `ncells` is the size of the aggregation group. Row sums of
the count matrix are library size, which is orders of magnitude larger,
and using them as the depth covariate answers a different question.

## Choose a cluster and fit the age model

[`select_de_cluster()`](https://myushen.github.io/cellNexusWorkshop/reference/select_de_cluster.md)
keeps clusters in which both age classes are represented and the design
can separate age from study. It prefers the cluster supported by the
most studies that contribute both age classes. The default screen asks
for 150 cells and 5 samples in each class. A 4000-cell cap may need the
slightly lower thresholds below. If the function stops, read the printed
screen and lower `min_samples_per_class`, or raise
`CELLNEXUS_WORKSHOP_MAX_CELLS` and rebuild `pb_se`.

\
`de_design`` ``<-`` ``~`` ``age_binary`` ``+`` ``log_cell_count`` ``+`` ``dataset_id`\
`de_contrast`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"age_binary"``, ``"old"``, ``"young"``)``)`\
`de_log2_fold_change_threshold`` ``<-`` ``1`\
`de_signature_max_genes`` ``<-`` ``25L`\
\
`epi_cluster_screen`` ``<-`` `[`select_de_cluster`](https://myushen.github.io/cellNexusWorkshop/reference/select_de_cluster.md)`(`\
`  ``pb_se``,`\
`  ``de_design``,`\
`  min_cells_per_class ``=`` ``100L``,`\
`  min_samples_per_class ``=`` ``4L`\
`)`\
`epi_cluster_de_top`` ``<-`` ``epi_cluster_screen``$``cluster`\
`epi_cluster_screen``$``screen`

> **Exercise 8.** From `epi_cluster_screen$screen`, which cluster was
> chosen? How many studies in that cluster contribute both age classes
> (`n_datasets_spanning_age`)? Quote the screen table with the result.
> The label alone is not reproducible.

Solution

\
`epi_cluster_de_top`\
`epi_cluster_screen``$``screen`

The chosen row is the first row of the screen after eligible clusters
are sorted by `n_datasets_spanning_age`, then by the number of
pseudobulk profiles. Quote `n_datasets_spanning_age` for that row. A
value of 1 means the age term rests on a single study.

The contrast is old versus young: positive log fold changes are higher
in donors over 50. DESeq2 tests against a log2 fold change of 1, so the
p-value asks whether the change exceeds a doubling. Say that when you
report the p-values. `sfType = "poscounts"` estimates size factors on
sparse pseudobulk, where the default median-of-ratios estimator can find
no gene detected in every profile.

Genes are kept when they reach 5 counts per million in at least half of
the profiles, evaluated with respect to age, so a gene expressed in one
age class survives.

\
`pb_se`` ``<-`` ``pb_se``[``, ``pb_se``$``epi_cluster`` ``==`` ``epi_cluster_de_top``]`\
`pb_se`` ``<-`` ``pb_se``[``, ``Matrix``::`[`colSums`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html)`(`[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``, ``"counts"``)``)`` ``>`` ``0``]`\
\
[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``, ``"counts"``)`` ``<-`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``, ``"counts"``)``)`\
`pb_se``$``sizeFactor`` ``<-`` ``NULL`\
\
`pb_se``$``age_binary`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`\
`  `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``pb_se``$``age_binary``, ``"old"``, ``"young"``)``,`\
`  levels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"young"``, ``"old"``)`\
`)`\
\
`pb_se`` ``<-`` ``pb_se`` ``|>`\
`  ``keep_abundant``(`\
`    formula_design ``=`` ``~`` ``age_binary``,`\
`    minimum_proportion ``=`` ``0.5``,`\
`    minimum_count_per_million ``=`` ``5`\
`  ``)`` ``|>`\
`  ``test_differential_expression``(`\
`    ``de_design``,`\
`    contrasts ``=`` ``de_contrast``,`\
`    method ``=`` ``"DESeq2"``,`\
`    omit_contrast_in_colnames ``=`` ``TRUE``,`\
`    test_above_log2_fold_change ``=`` ``de_log2_fold_change_threshold``,`\
`    sfType ``=`` ``"poscounts"`\
`  ``)`` ``|>`\
`  ``adjust_abundance``(`\
`    .factor_unwanted ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``log_cell_count``, ``dataset_id``)``,`\
`    .factor_of_interest ``=`` ``age_binary``,`\
`    method ``=`` ``"limma_remove_batch_effect"`\
`  ``)`

`adjust_abundance()` is for the figure. The test used the unadjusted
counts. Batch-adjusted values stay out of any further test.

\
`de_results`` ``<-`` `[`rowData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``pb_se``)`` ``|>`\
`  ``tibble``::`[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``rownames ``=`` ``".feature"``)`` ``|>`\
`  ``dplyr``::`[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(`\
`    gene_symbol ``=`` ``AnnotationDbi``::``mapIds``(`\
`      ``org.Hs.eg.db``,`\
`      keys ``=`` ``stringr``::``str_replace``(``.feature``, ``"\\.[0-9]+$"``, ``""``)``,`\
`      keytype ``=`` ``"ENSEMBL"``,`\
`      column ``=`` ``"SYMBOL"``,`\
`      multiVals ``=`` ``"first"`\
`    ``)`` ``|>`\
`      `[`unname`](https://rdrr.io/r/base/unname.html)`(``)`` ``|>`\
`      ``dplyr``::`[`coalesce`](https://dplyr.tidyverse.org/reference/coalesce.html)`(``.feature``)`\
`  ``)`\
\
`signature_up`` ``<-`` ``de_results`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``pvalue`` ``<`` ``0.05``, ``log2FoldChange`` ``>`` ``0``)`` ``|>`\
`  ``dplyr``::`[`slice_min`](https://dplyr.tidyverse.org/reference/slice.html)`(``pvalue``, n ``=`` ``de_signature_max_genes``)`\
\
`signature_down`` ``<-`` ``de_results`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``pvalue`` ``<`` ``0.05``, ``log2FoldChange`` ``<`` ``0``)`` ``|>`\
`  ``dplyr``::`[`slice_min`](https://dplyr.tidyverse.org/reference/slice.html)`(``pvalue``, n ``=`` ``de_signature_max_genes``)`\
\
`signature_up`

The signature uses a nominal p-value cutoff, capped at 25 genes in each
direction. It is a gene set for scoring. A claim about an individual
gene needs the false-discovery column as well. Report both.

\
`marker_genes`` ``<-`` ``dplyr``::`[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(`\
`  ``dplyr``::`[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``signature_up``, direction ``=`` ``"up"``)``,`\
`  ``dplyr``::`[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``signature_down``, direction ``=`` ``"down"``)`\
`)`` ``|>`\
`  ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``.feature``, ``gene_symbol``, ``direction``)`\
\
`marker_plot_df`` ``<-`` ``pb_se``[``marker_genes``$``.feature``, ``]`` ``|>`\
`  `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
`  ``dplyr``::`[`left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html)`(``marker_genes``, by ``=`` ``".feature"``)`

\
`ggplot``(``marker_plot_df``, ``aes``(``age_bin``, `[`log1p`](https://rdrr.io/r/base/Log.html)`(``counts_adjusted``)``, fill ``=`` ``age_bin``)``)`` ``+`\
`  ``geom_boxplot``(``outlier.size ``=`` ``0.4``, width ``=`` ``0.65``)`` ``+`\
`  ``facet_wrap``(``~`` ``gene_symbol`` ``+`` ``direction``, scales ``=`` ``"free_y"``, ncol ``=`` ``6``)`` ``+`\
`  ``theme_bw``(``)`` ``+`\
`  ``theme``(``legend.position ``=`` ``"none"``, strip.text ``=`` ``element_text``(``size ``=`` ``7``)``)`` ``+`\
`  ``labs``(`\
`    x ``=`` ``"Age"``,`\
`    y ``=`` ``"log1p(counts_adjusted)"``,`\
`    title ``=`` `[`paste`](https://rdrr.io/r/base/paste.html)`(``"Cluster"``, ``epi_cluster_de_top``)`\
`  ``)`

![Example marker boxplots. Each point is a donor-level pseudobulk
profile in the cluster the screen selected. Abundance is log1p of the
batch-adjusted count.](images/marker-boxplot.png)

Example marker boxplots. Each point is a donor-level pseudobulk profile
in the cluster the screen selected. Abundance is log1p of the
batch-adjusted count.

> **Exercise 9.** In `de_contrast`, which age class is the numerator?
> Check one gene in `signature_up` against its boxplot and confirm the
> direction matches.

Solution

`c("age_binary", "old", "young")` tests old against young. Positive
`log2FoldChange` means higher abundance in donors over 50. The boxplot
for a gene in `signature_up` should sit higher in the `50-79` bin. If it
sits lower, the factor levels and the contrast disagree.

## Score the signature on single cells

If time remains. The score is the mean rank of the up genes minus the
mean rank of the down genes, on the cells of the cluster that was
tested. A shift between age classes is expected, because the genes were
chosen on these donors. Treat it as a consistency check: a reversed or
empty shift means a reversed contrast, the wrong cluster, or an
identifier mismatch. Report how many signature genes were found in the
matrix. Ensembl identifiers with a version suffix are a common way to
lose most of the set. Generalisation needs donors held out before gene
selection, or another cohort.

\
[`library`](https://rdrr.io/r/base/library.html)`(``singscore``)`\
\
`sce`` ``<-`` ``sce_epi`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``epi_cluster`` ``==`` ``epi_cluster_de_top``)`\
\
`rank_mat`` ``<-`` ``rankGenes``(`[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(`[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``sce``, ``"logcounts"``)``)``)`\
`sce``$``aging_signature_score`` ``<-`` ``simpleScore``(`\
`  ``rank_mat``,`\
`  upSet ``=`` ``signature_up``$``.feature``,`\
`  downSet ``=`` ``if`` ``(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``signature_down``)`` ``>`` ``0``)`` ``signature_down``$``.feature`` ``else`` ``NULL`\
`)``$``TotalScore`\
\
[`mean`](https://rdrr.io/r/base/mean.html)`(``signature_up``$``.feature`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``sce``)``)`\
[`mean`](https://rdrr.io/r/base/mean.html)`(``signature_down``$``.feature`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``sce``)``)`

\
`sce`` ``|>`\
`  ``ggplot``(``aes``(``UMAP1``, ``UMAP2``, colour ``=`` ``aging_signature_score``)``)`` ``+`\
`  ``geom_point``(``size ``=`` ``0.3``, alpha ``=`` ``0.6``)`` ``+`\
`  ``scale_colour_viridis_c``(``option ``=`` ``"magma"``)`` ``+`\
`  ``facet_wrap``(``~`` ``age_bin``, ncol ``=`` ``2``)`` ``+`\
`  ``coord_equal``(``)`` ``+`\
`  ``theme_bw``(``)`` ``+`\
`  ``labs``(``colour ``=`` ``"Ageing score"``, title ``=`` `[`paste`](https://rdrr.io/r/base/paste.html)`(``"Cluster"``, ``epi_cluster_de_top``)``)`

\
`sce`` ``|>`\
`  ``ggplot``(``aes``(``aging_signature_score``, fill ``=`` ``age_bin``, colour ``=`` ``age_bin``)``)`` ``+`\
`  ``geom_density``(``alpha ``=`` ``0.35``, linewidth ``=`` ``0.6``)`` ``+`\
`  ``theme_bw``(``)`` ``+`\
`  ``labs``(``x ``=`` ``"Ageing signature score"``, y ``=`` ``"Density"``, fill ``=`` ``"Age"``, colour ``=`` ``"Age"``)`

## What you now have

One session, from a cellNexus query to two donor-level answers:

- whether immune proportions differ between donors under and over 50, at
  fixed study (`sccomp_res`)
- whether epithelial transcription in one cluster differs between those
  same donors (`de_results`)

Save the script from
`system.file("scripts", "workshop_respiratory_ageing.R", package = "cellNexusWorkshop")`.
Changing the dataset list, the cell-type vectors, or the age cut is how
the same workflow moves to another tissue or another question.

Solutions:
[Solutions](https://myushen.github.io/cellNexusWorkshop/articles/solutions.md).
