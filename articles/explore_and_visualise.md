# 2. Explore and visualise

Continue in the session from [Query the
atlas](https://myushen.github.io/cellNexusWorkshop/articles/query_the_atlas.md).
This block downloads a capped SingleCellExperiment, checks that immune
and epithelial cells separate, and shows the same cells as a Seurat
object.

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellNexus``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`cellNexusWorkshop`](https://myushen.github.io/cellNexusWorkshop)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ggplot2`](https://ggplot2.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scater``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``scran``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`harmony`](https://github.com/immunogenomics/harmony)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``SingleCellExperiment``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``tidySingleCellExperiment``)`\
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

## Download a SingleCellExperiment

`get_single_cell_experiment()` reads the count matrices for the cells in
a metadata table and returns one `SingleCellExperiment`. The workshop
cap keeps the download inside a session. The example figures on this
site use 12,000 cells. Set `CELLNEXUS_WORKSHOP_MAX_CELLS=12000` before
this chunk to match them. Allow several minutes, and about 32 GB of RAM
for the larger cap.

Cells are subsampled within sample, and the number of samples is capped,
so a few large libraries do not dominate the highly variable genes.

\
`max_cells_qc`` ``<-`` `[`as.integer`](https://rdrr.io/r/base/integer.html)`(`[`Sys.getenv`](https://rdrr.io/r/base/Sys.getenv.html)`(``"CELLNEXUS_WORKSHOP_MAX_CELLS"``, ``"4000"``)``)`\
`cells_per_sample_qc`` ``<-`` ``80L`\
`n_samples_qc`` ``<-`` `[`max`](https://rdrr.io/r/base/Extremes.html)`(``1L``, ``max_cells_qc`` `[`%/%`](https://rdrr.io/r/base/Arithmetic.html)` ``cells_per_sample_qc``)`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`samples_qc`` ``<-`` ``dplyr``::`[`distinct`](https://dplyr.tidyverse.org/reference/distinct.html)`(``obs_lung``, ``sample_id``)`\
`samples_qc`` ``<-`` ``dplyr``::`[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)`(``samples_qc``, n ``=`` `[`min`](https://rdrr.io/r/base/Extremes.html)`(``n_samples_qc``, `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``samples_qc``)``)``)`\
\
`sampled_qc`` ``<-`` ``obs_lung`` ``|>`\
`  ``dplyr``::`[`semi_join`](https://dplyr.tidyverse.org/reference/filter-joins.html)`(``samples_qc``, by ``=`` ``"sample_id"``)`` ``|>`\
`  ``dplyr``::`[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``sample_id``)`` ``|>`\
`  ``dplyr``::`[`slice`](https://dplyr.tidyverse.org/reference/slice.html)`(`[`sample.int`](https://rdrr.io/r/base/sample.html)`(``dplyr``::`[`n`](https://dplyr.tidyverse.org/reference/context.html)`(``)``, `[`min`](https://rdrr.io/r/base/Extremes.html)`(``cells_per_sample_qc``, ``dplyr``::`[`n`](https://dplyr.tidyverse.org/reference/context.html)`(``)``)``)``)`` ``|>`\
`  ``dplyr``::`[`ungroup`](https://dplyr.tidyverse.org/reference/group_by.html)`(``)`\
\
`sce_qc`` ``<-`` ``sampled_qc`` ``|>`\
`  ``get_single_cell_experiment``(``assays ``=`` ``"counts"``)`` ``|>`\
`  ``logNormCounts``(``)`\
\
`sce_qc`

With `tidySingleCellExperiment` loaded, the object prints as a tibble
and accepts `dplyr` verbs and `ggplot2`. The underlying assay is
unchanged.

\
`sce_qc`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``cell_type_unified_ensemble``, sort ``=`` ``TRUE``)`

## Integrate studies and colour by CD45

Highly variable genes are selected per sample
(`modelGeneVar(block = sample_id)`), so genes that vary within a donor
rank above genes that merely differ between libraries (Lun et al. 2016).
Harmony then removes the study axis from the PCA (Korsunsky et al.
2019). The UMAP is drawn from that embedding.

CD45 (`PTPRC`, Ensembl `ENSG00000081237`) marks haematopoietic cells
independently of the harmonised label. Immune cells should occupy the
high-CD45 region.

\
`hvgs_qc`` ``<-`` ``sce_qc`` ``|>`\
`  ``modelGeneVar``(``block ``=`` ``sce_qc``$``sample_id``)`` ``|>`\
`  ``getTopHVGs``(``n ``=`` ``2000``)`\
\
`sce_qc`` ``<-`` ``sce_qc`` ``|>`\
`  ``runPCA``(``ncomponents ``=`` ``30``, subset_row ``=`` ``hvgs_qc``)`` ``|>`\
`  ``RunHarmony``(`\
`    group.by.vars ``=`` ``"dataset_id"``,`\
`    dims.use ``=`` ``1``:``30``,`\
`    verbose ``=`` ``FALSE``,`\
`    reduction.save ``=`` ``"HARMONY"`\
`  ``)`` ``|>`\
`  ``runUMAP``(``dimred ``=`` ``"HARMONY"``, n_neighbors ``=`` ``15``)`\
\
`cd45_id`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  `[`grep`](https://rdrr.io/r/base/grep.html)`(``"^ENSG00000081237(\\.|$)"``, `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``sce_qc``)``, value ``=`` ``TRUE``)``,`\
`  `[`intersect`](https://generics.r-lib.org/reference/setops.html)`(``"PTPRC"``, `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``sce_qc``)``)`\
`)`\
\
`if`` ``(`[`length`](https://rdrr.io/r/base/length.html)`(``cd45_id``)``)`` ``{`\
`  ``sce_qc``$``cd45`` ``<-`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``logcounts``(``sce_qc``)``[``cd45_id``[[``1``]``]``, ``]``)`\
`}`` ``else`` ``{`\
`  ``sce_cd45`` ``<-`` ``sampled_qc`` ``|>`\
`    ``get_single_cell_experiment``(``assays ``=`` ``"counts"``, features ``=`` ``"ENSG00000081237"``)`\
`  ``sce_qc``$``cd45`` ``<-`` ``NA_real_`\
`  ``i`` ``<-`` `[`match`](https://rdrr.io/r/base/match.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sce_qc``)``, `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sce_cd45``)``)`\
`  ``sce_qc``$``cd45``[``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``i``)``]`` ``<-`` `[`log1p`](https://rdrr.io/r/base/Log.html)`(`\
`    `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``SummarizedExperiment``::`[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``sce_cd45``, ``"counts"``)``[``1``, ``i``[``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``i``)``]``]``)`\
`  ``)`\
`  `[`rm`](https://rdrr.io/r/base/rm.html)`(``sce_cd45``)`\
`}`

\
`sce_qc`` ``|>`\
`  ``ggplot``(``aes``(``UMAP1``, ``UMAP2``, colour ``=`` ``cd45``)``)`` ``+`\
`  ``geom_point``(``size ``=`` ``0.3``, alpha ``=`` ``0.6``)`` ``+`\
`  ``scale_colour_viridis_c``(``option ``=`` ``"magma"``)`` ``+`\
`  ``coord_equal``(``)`` ``+`\
`  ``theme_bw``(``)`` ``+`\
`  ``labs``(``colour ``=`` ``"CD45 (PTPRC)"``)`

![Example Harmony UMAP, coloured by CD45. Immune cells occupy the
high-CD45 region.](images/cd45-umap.png)

Example Harmony UMAP, coloured by CD45. Immune cells occupy the
high-CD45 region.

The same grammar colours the embedding by the harmonised label.

\
`sce_qc`` ``|>`\
`  ``ggplot``(``aes``(``UMAP1``, ``UMAP2``, colour ``=`` ``cell_type_unified_ensemble``)``)`` ``+`\
`  ``geom_point``(``size ``=`` ``0.3``, alpha ``=`` ``0.6``)`` ``+`\
`  ``coord_equal``(``)`` ``+`\
`  ``theme_bw``(``)`` ``+`\
`  ``labs``(``colour ``=`` ``"Cell type"``)`

> **Exercise 4.** Colour the same UMAP by `dataset_id`. After Harmony,
> study identity should be mixed within a cell-type region. Islands that
> follow study boundaries are under-correction. Harmony adjusts the
> embedding used for the picture. It does not adjust the counts, so the
> models in the next two articles still include `dataset_id`.

Solution

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

> **Exercise 5.** Using `sce_qc`, count cells per
> `cell_type_unified_ensemble`. Which of
> [`immune_cell_types()`](https://myushen.github.io/cellNexusWorkshop/reference/immune_cell_types.md)
> are missing from this subsample, and what does that imply for a
> composition estimate built from these cells alone?

Solution

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

## The same cells as a Seurat object

`get_seurat()` is `get_single_cell_experiment()` followed by the Seurat
conversion. Pass a small metadata table. A handful of cells spread
across many files downloads a handful of files.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyseurat`](https://github.com/stemangiola/tidyseurat)`)`\
\
`one_study`` ``<-`` ``obs_lung`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``dataset_id`` ``==`` ``dplyr``::`[`first`](https://dplyr.tidyverse.org/reference/nth.html)`(``dataset_id``)``)`` ``|>`\
`  ``dplyr``::`[`slice_head`](https://dplyr.tidyverse.org/reference/slice.html)`(``n ``=`` ``100``)`\
\
`seurat_demo`` ``<-`` ``one_study`` ``|>`\
`  ``get_seurat``(``assays ``=`` ``"counts"``)`\
\
`seurat_demo`

`tidyseurat` gives that object the same tibble interface you just used
on the SingleCellExperiment. The rest of the workshop stays with
SingleCellExperiment, because the pseudobulk tools speak that language.

Next: [Compositional
analysis](https://myushen.github.io/cellNexusWorkshop/articles/compositional_analysis.md).
