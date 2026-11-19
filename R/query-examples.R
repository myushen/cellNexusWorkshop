library(cellNexus)

# Load cellnexus metadata
metadata <- get_metadata()

# QC
metadata <- metadata |> 
  filter(empty_droplet == FALSE, alive == TRUE, scDblFinder.class != "doublet")

# Filter t cd4 cells. View the number of datasets contains cd4 t cells
metadata |>
  filter(cell_type_unified_ensemble == "t cd4") |>
  distinct(dataset_id) |>
  collect()

# Filter t cd4 cells. Explore which study has the most single cell annotated as cd4 T 
metadata |>
  filter(cell_type_unified_ensemble == "t cd4") |>
  count(dataset_id) |>
  arrange(desc(n))

# Join annotations. For Mac user
census_metadata <- cellNexus:::get_census_metadata("2024-07-01")
con <- dbplyr::remote_con(metadata)
duckdb::duckdb_register_arrow(con, "census_metadata", census_metadata)
tbl(con, "census_metadata")
metadata <- metadata |>
  left_join(
    tbl(con, "census_metadata") |>
      select(observation_joinid, dataset_id,
             tissue, self_reported_ethnicity,
             cell_type, assay, disease, sex),
    by = c("dataset_id", "observation_joinid")
    )
    
# Join annotations. For windows user
url <- "https://github.com/myushen/cellNexusWorkshop/releases/download/census-data-for-windows-v1/census_metadata.parquet"
dest <- "census_metadata.parquet"
if (!file.exists(dest)) {
  download.file(
    url,
    destfile = dest,
    mode = "wb",
  )
}
census_metadata <- arrow::read_parquet("census_metadata.parquet")
    
metadata <- metadata |>
      left_join(census_metadata |>
                  select(observation_joinid, dataset_id,
                         tissue, self_reported_ethnicity,
                         cell_type, assay, disease, sex),
                by = c("dataset_id", "observation_joinid"),
                copy = T
      )

# Define Lung cohort
## Calculate the number of cells for each disease type in lung tissue
metadata |>
  filter(tissue == "lung") |>
  count(disease)

## Calculate the number of cells for healthy European people sequenced by 10x 3' v2
metadata |>
  filter(tissue == "lung", disease == "normal",
         self_reported_ethnicity == "European",
         assay == "10x 3' v2") |>
  count(disease)

## Define the cohort
lung_cohort <- metadata |>
  filter(tissue == "lung", disease == "normal",
         self_reported_ethnicity == "European",
         assay == "10x 3' v2") |>
  count(disease)

# View annotations
lung_cohort |> glimpse()

# Or load cohort to disk and view annotations - run with cautious! Do not make crash memory ;)
lung_cohort |> collect() |> View()

# Download expression data for lung cohort
sce <- lung_cohort |> get_single_cell_experiment()
sce_cpm <- lung_cohort |> get_single_cell_experiment(assay = "cpm")
pb <- lung_cohort |> get_pseudobulk(as_SummarizedExperiment = T)

# Now its your turn to explore!
my_metadata <- get_metadata() |> ...



      