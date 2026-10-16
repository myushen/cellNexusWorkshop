# 3. Compositional analysis

Continue in the same session. This block fits a donor-level model of
immune cell-type proportions with `sccomp` (Mangiola et al. 2023).

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`cellNexusWorkshop`](https://myushen.github.io/cellNexusWorkshop)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``sccomp``)`\
\
`if`` ``(``!`[`exists`](https://rdrr.io/r/base/exists.html)`(``"obs_lung"``)``)`` ``{`\
`  `[`library`](https://rdrr.io/r/base/library.html)`(``cellNexus``)`\
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

## Why a composition model

A proportion vector sums to one, so an increase in one cell type forces
a decrease somewhere else. Counts are overdispersed relative to a
binomial. One dissociated sample can move an entire effect. `sccomp`
models the sample composition jointly, estimates a mean-variability
association, and has an explicit outlier pass.

The denominator is the immune compartment. Epithelial recovery depends
on how hard the tissue was digested. A study that recovered many
pneumocytes would push every immune proportion down, and that shift
would enter the model as if it were age. Epithelial cells stay in
`obs_lung` for the expression article.

\
`obs_composition`` ``<-`` ``obs_lung`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``cell_type_unified_ensemble`` `[`%in%`](https://rdrr.io/r/base/match.html)` `[`immune_cell_types`](https://myushen.github.io/cellNexusWorkshop/reference/immune_cell_types.md)`(``)``)`` ``|>`\
`  ``# sccomp stores covariates in one column. A logical cannot share that`\
`  ``# column with a character study id. "FALSE" is the reference level.`\
`  ``dplyr``::`[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``age_binary ``=`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``age_binary``)``)`\
\
`obs_composition`` ``|>`\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``cell_type_unified_ensemble``, name ``=`` ``"n_cells"``)`` ``|>`\
`  ``dplyr``::`[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(``dplyr``::`[`desc`](https://dplyr.tidyverse.org/reference/desc.html)`(``n_cells``)``)`

## Fit

The formula is `~ age_binary + dataset_id`. `sample_id` is the
replicate. `cell_type_unified_ensemble` is the group. `dataset_id` is a
fixed effect, so the age contrast is within study. That is why Exercise
3 matters: a study with one age class competes with the age term.

`bimodal_mean_variability_association = TRUE` lets rare and abundant
cell types follow different variability regimes.

The session default is `pathfinder`, which is fast enough to finish in
the room. A result you would report uses Hamiltonian Monte Carlo
(`inference_method = "hmc"`). Pathfinder intervals are a working
approximation. Refit with `"hmc"` before you quote an effect, and read
the sampler diagnostics on that fit. Set
`CELLNEXUS_WORKSHOP_INFERENCE=hmc` for that fit. The first call compiles
the Stan model and can take several minutes on top of the fit.

\
`inference_method`` ``<-`` `[`Sys.getenv`](https://rdrr.io/r/base/Sys.getenv.html)`(``"CELLNEXUS_WORKSHOP_INFERENCE"``, ``"pathfinder"``)`\
\
`sccomp_res`` ``<-`` ``obs_composition`` ``|>`\
`  ``sccomp_estimate``(`\
`    formula_composition ``=`` ``~`` ``age_binary`` ``+`` ``dataset_id``,`\
`    sample ``=`` ``"sample_id"``,`\
`    cell_group ``=`` ``"cell_type_unified_ensemble"``,`\
`    bimodal_mean_variability_association ``=`` ``TRUE``,`\
`    inference_method ``=`` ``inference_method`\
`  ``)`` ``|>`\
`  ``sccomp_remove_outliers``(``)`` ``|>`\
`  ``sccomp_test``(``)`

Outlier removal refits the model with aberrant samples down-weighted.
Removal comes before the test, because the test reads the refitted
model. Look at which samples were flagged before accepting the
correction. A flag concentrated in one study is a finding of its own.

\
`sccomp_res`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``factor`` ``==`` ``"age_binary"``)`` ``|>`\
`  ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``cell_type_unified_ensemble``, ``dplyr``::`[`starts_with`](https://tidyselect.r-lib.org/reference/starts_with.html)`(``"c_"``)``, ``c_FDR``)`

Columns whose names start with `c_` carry the age effect, its interval,
`c_FDR`, and sampler diagnostics. On an HMC fit, `c_rhat` near 1 and a
large effective sample size are the checks you want. The effect is on
the model’s composition scale. Exponentiating it as a fold change in
proportion is the wrong transformation. Say “the proportion of this cell
type”.

An effect is a change in share. Because shares sum to one, one real
change produces mirror-image changes in other labels. The label set is
part of the result: this fit is a shift inside the immune compartment.

## Plot

`age_binary` has two levels, so the summary is a boxplot.
`remove_unwanted_effects = TRUE` shows proportions after the study term
has been regressed out, which is the visual counterpart of the model.
The blue box is the posterior predictive check: overlap between the
observed points and the predictive distribution is evidence that the
model describes these data.

\
`sccomp_boxplot``(`\
`  ``sccomp_res``,`\
`  factor ``=`` ``"age_binary"``,`\
`  significance_threshold ``=`` ``0.05``,`\
`  remove_unwanted_effects ``=`` ``TRUE`\
`)`

![Example sccomp boxplots from a Hamiltonian Monte Carlo fit. Each point
is one sample, after the study term has been
removed.](images/composition-boxplot.png)

Example sccomp boxplots from a Hamiltonian Monte Carlo fit. Each point
is one sample, after the study term has been removed.

> **Exercise 6.** Which cell types have `c_FDR < 0.05` for
> `factor == "age_binary"`? For each, is the estimated effect higher or
> lower in donors over 50?

Solution

\
`sccomp_res`` ``|>`\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``factor`` ``==`` ``"age_binary"``, ``c_FDR`` ``<`` ``0.05``)`` ``|>`\
`  ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``cell_type_unified_ensemble``, ``dplyr``::`[`starts_with`](https://tidyselect.r-lib.org/reference/starts_with.html)`(``"c_"``)``, ``c_FDR``)`

Read the sign of the effect column together with the boxplot. The
contrast is donors over 50 against the reference level `"FALSE"` (50 and
under), because `age_binary` was stored as character with `"FALSE"`
first.

> **Exercise 7.** Write one sentence on why `epithelial`, `pneumocyte`,
> and `secretory` are absent from this fit. What question would the
> model answer if they were added back?

Solution

Epithelial, pneumocyte, and secretory labels are left out so the
proportions are shares of the immune compartment. Dissociation
efficiency changes epithelial yield from study to study, and that
technical shift would move every immune proportion. Adding the
epithelial labels back answers a different question: change in share
among immune and epithelial cells together. Report the label set with
the effects.

If the fit is too slow for the time left, refit on a shorter label list
(for example macrophage, CD14 monocytes, NK, naive and effector-memory T
cells, naive B cells) and say so. The denominator changes, so the
effects are a different question.

Next: [Gene
expression](https://myushen.github.io/cellNexusWorkshop/articles/gene_expression.md).
