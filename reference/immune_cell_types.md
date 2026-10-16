# Immune labels used as the composition denominator.

Composition is modelled inside the immune compartment. Epithelial yield
in dissociated lung depends on the digestion protocol, so epithelial
labels are kept out of this list and analysed in the expression article.

## Usage

``` r
immune_cell_types()
```

## Value

A character vector of `cell_type_unified_ensemble` labels.

## Examples

``` r
immune_cell_types()
#>  [1] "cd16 mono"       "nk"              "mait"            "cd14 mono"      
#>  [5] "cdc"             "b memory"        "pdc"             "nkt"            
#>  [9] "cd4 th1/th17 em" "plasma"          "cd4 th17 em"     "b naive"        
#> [13] "macrophage"      "cd4 th2 em"      "cd8 tcm"         "cd8 naive"      
#> [17] "cd8 tem"         "tgd"             "cd4 fh em"       "cd4 th1 em"     
#> [21] "cd4 tem"         "cd4 naive"       "mast"            "treg"           
#> [25] "cd4 tcm"         "granulocyte"    
```
