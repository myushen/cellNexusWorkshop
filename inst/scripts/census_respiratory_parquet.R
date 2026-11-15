# Census obs slice for the respiratory workshop
#
# Run this on Linux or macOS, from the package root:
#   Rscript inst/scripts/census_respiratory_parquet.R
#
# cellxgene.census does not run on Windows, and the workshop virtual
# machine is Windows. This script reads the CELLxGENE Census once and
# writes the rows the query chapter joins. The workshop then reads
# inst/extdata/census_respiratory_2024-07-01.parquet instead of opening
# the Census.
#
# The slice matches the query chapter:
#   snapshot 2024-07-01, the same snapshot as cellNexus::get_metadata()
#   is_primary_data == TRUE, the same row filter as
#     cellNexus:::get_census_metadata()
#   the nine dataset_id values used in vignettes/query_the_atlas.Rmd
#   the eight columns that join adds: observation_joinid, dataset_id,
#     tissue, self_reported_ethnicity, cell_type, assay, disease, sex
#
# Tissue is not filtered here. The chapter counts tissue and disease
# before it keeps the respiratory labels.
#
# The value_filter string is R syntax. tiledbsoma parses it with
# str2lang(), so `and` is not valid and `%in%` is required.
#
# Requires cellxgene.census and arrow. Reinstall the package after this
# script runs so system.file() can find the new parquet.

library(cellxgene.census)
library(arrow)

census_version <- "2024-07-01"

dataset_ids <- c(
  "9f222629-9e39-47d0-b83f-e08d610c7479",
  "093d3bfe-6f0f-4ac0-a7a1-829f94d0a49f",
  "8c42cfd0-0b0a-46d5-910c-fc833d83c45e",
  "53d208b0-2cfd-4366-9866-c3c6114081bc",
  "3de0ad6d-4378-4f62-b37b-ec0b75a50d94",
  "d8da613f-e681-4c69-b463-e94f5e66847f",
  "e04daea4-4412-45b5-989e-76a9be070a89",
  "2adb1f8a-a6b1-4909-8ee8-484814e2d4bf",
  "2672b679-8048-4f5e-9786-f1b196ccfd08"
)

value_filter <- sprintf(
  "is_primary_data == TRUE & dataset_id %%in%% c(%s)",
  paste(sprintf("\"%s\"", dataset_ids), collapse = ", ")
)

census <- open_soma(census_version = census_version)
on.exit(census$close(), add = TRUE)

census_metadata <- census$
  get("census_data")$
  get("homo_sapiens")$
  get("obs")$
  read(
    value_filter = value_filter,
    column_names = c(
      "observation_joinid", "dataset_id", "tissue",
      "self_reported_ethnicity", "cell_type", "assay",
      "disease", "sex"
    )
  )$
  concat()

out <- "inst/extdata/census_respiratory_2024-07-01.parquet"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write_parquet(census_metadata, out, compression = "zstd", compression_level = 19)

message(
  "Wrote ", nrow(census_metadata), " rows to ", out,
  " (", round(file.info(out)$size / 1024^2, 1), " MB)"
)
