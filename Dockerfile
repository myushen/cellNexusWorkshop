FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

RUN Rscript -e "install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', 'https://cloud.r-project.org'))" \
 && Rscript -e "cmdstanr::install_cmdstan(cores = 2)"

RUN Rscript -e "BiocManager::install(c( \
      'dplyr', 'ggplot2', 'patchwork', 'stringr', 'Matrix', \
      'SummarizedExperiment', 'SingleCellExperiment', \
      'scater', 'scran', 'scuttle', 'scrapper', 'bluster', 'harmony', \
      'DESeq2', 'org.Hs.eg.db', 'AnnotationDbi', 'singscore', \
      'tidySingleCellExperiment', 'tidySummarizedExperiment', \
      'tidybulk', 'tidyseurat', 'tidyomics', 'Seurat', \
      'remotes', 'pkgdown', 'knitr', 'rmarkdown', 'dbplyr', 'duckdb' \
    ), ask = FALSE, update = FALSE)"

RUN Rscript -e "remotes::install_github(c( \
      'MangiolaLaboratory/cellNexus', \
      'MangiolaLaboratory/sccomp' \
    ), upgrade = 'never')"

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "remotes::install_local('/home/rstudio', dependencies = FALSE, upgrade = 'never', build_vignettes = TRUE)"
