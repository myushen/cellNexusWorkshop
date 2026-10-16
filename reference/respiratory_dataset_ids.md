# Dataset identifiers for the respiratory reference cohort.

Nine curated respiratory studies from the Human Cell Atlas. Pinning the
identifiers keeps the workshop cohort stable when the atlas gains
studies. Identifiers use the hyphenated form stored by `get_metadata()`.
[`annotate_cell_nexus()`](https://myushen.github.io/cellNexusWorkshop/reference/annotate_cell_nexus.md)
rewrites hyphens to underscores so the same identifier can be placed in
a model formula.

## Usage

``` r
respiratory_dataset_ids()
```

## Value

A character vector of `dataset_id` values.

## Examples

``` r
respiratory_dataset_ids()
#> [1] "9f222629-9e39-47d0-b83f-e08d610c7479"
#> [2] "093d3bfe-6f0f-4ac0-a7a1-829f94d0a49f"
#> [3] "8c42cfd0-0b0a-46d5-910c-fc833d83c45e"
#> [4] "53d208b0-2cfd-4366-9866-c3c6114081bc"
#> [5] "3de0ad6d-4378-4f62-b37b-ec0b75a50d94"
#> [6] "d8da613f-e681-4c69-b463-e94f5e66847f"
#> [7] "e04daea4-4412-45b5-989e-76a9be070a89"
#> [8] "2adb1f8a-a6b1-4909-8ee8-484814e2d4bf"
#> [9] "2672b679-8048-4f5e-9786-f1b196ccfd08"
```
