# Prepare a neurogam result for model-independent plotting

Internal adapter. A saved result may lack its model field, which is not
needed for plots built from saved prediction/cluster tables. P2.1
validation remains strict; any repair is applied to this function's
local copy only.

## Usage

``` r
neurogam_plot_view(x)
```
