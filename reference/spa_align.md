# Stratigraphic Plug Alignment (SPA)

Linearly interpolates plug-based measurements (e.g., TOC, porosity, XRD)
onto a high-resolution reference depth grid (e.g., XRF). The procedure
uses base R's [`approx()`](https://rdrr.io/r/stats/approxfun.html) with
`rule = 1` (default) to prevent extrapolation beyond the observed depth
range, ensuring stratigraphically consistent alignment of all datasets.

## Usage

``` r
spa_align(
  ref,
  ...,
  depth_col = "Depth_m",
  rule = 1,
  add_suffix = TRUE,
  trim = FALSE
)
```

## Arguments

- ref:

  A data.frame containing the reference depth grid and (optionally)
  high-resolution variables (e.g., XRF). Must contain the depth column
  specified in `depth_col`.

- ...:

  One or more named data.frames containing plug-based measurements to be
  interpolated (e.g., `xrd = xrd_df`, `plugs = plug_df`).

- depth_col:

  A character string giving the name of the depth column shared by all
  input datasets. Defaults to `"Depth_m"`. The depth column may use any
  unit (e.g., meters, feet, centimeters); `"Depth_m"` is only a column
  label and does not require depths to be in meters. However, all input
  datasets must use the same depth unit for interpolation to be
  meaningful.

- rule:

  Integer passed to [`approx()`](https://rdrr.io/r/stats/approxfun.html)
  (default `1`). `rule = 1` prevents extrapolation outside the observed
  depth range. `rule = 2` performs constant extrapolation (end-point
  hold).

- add_suffix:

  Logical; if `TRUE`, variable names are suffixed with the dataset name
  (e.g., `TOC_plugs`, `Quartz_xrd`).

- trim:

  Logical; if `TRUE`, drops rows where any interpolated column is `NA`.
  Default `FALSE` (preserves all reference depths).

## Value

A data.frame containing the reference depth grid and interpolated
variables aligned to the same resolution.

## Details

SPA is intended for vertically ordered core or log data, where
measurements are indexed by depth along a stratigraphic profile.

## Examples

``` r
ref <- data.frame(Depth_m = 0:10, Ca = runif(11, 100, 200))
xrd <- data.frame(Depth_m = c(2, 5, 7), Quartz = c(54, 60, 58))
plugs <- data.frame(Depth_m = c(3, 7, 9), TOC = c(3.0, 3.3, 3.5))

# Default: preserves reference depth grid
aligned <- spa_align(ref, xrd = xrd, plugs = plugs)
head(aligned)
#>   Depth_m       Ca Quartz_xrd TOC_plugs
#> 1       0 160.0761         NA        NA
#> 2       1 115.7208         NA        NA
#> 3       2 100.7399         54        NA
#> 4       3 146.6393         56     3.000
#> 5       4 149.7777         58     3.075
#> 6       5 128.9767         60     3.150

# Overlap-only alignment
aligned_overlap <- spa_align(ref, xrd = xrd, plugs = plugs, trim = TRUE)
```
