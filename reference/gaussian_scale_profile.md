# Profile-likelihood interval for a Gaussian smoothing scale

Fixes one scale at each grid value and reoptimizes all remaining
conductance and measurement parameters on the fitted graph.

## Usage

``` r
gaussian_scale_profile(fit, layer, n = 25, level = 0.95)
```

## Arguments

- fit:

  An unslimmed Gaussian scale fit with its original model closures.

- layer:

  Name of the raster layer whose scale is profiled.

- n:

  Number of grid points across the fitted scale bounds, at least five.
  The fitted estimate is also included if it is not already on the grid.

- level:

  Confidence level between zero and one.

## Value

A list with layer, estimate, level, interval, grid, and fit. The grid
contains scale values and profile log likelihoods. The returned fit
stores the interval, so `confint(result$fit)` uses it at the matching
level.

## Details

The interval is the connected likelihood-ratio support region containing
the fitted estimate, using the one-degree-of-freedom chi-square cutoff.
Grid crossings are refined by root finding. Increase `n` to check for
additional support regions. Limits stop at the original bounds. At a
boundary the usual chi-square calibration is an approximation. Failed
optimizations stop the calculation instead of silently dropping profile
points. Both conductance and final measurement-parameter convergence are
required. This can be substantially slower than Wald intervals.

## Worked example

See
[`vignette("v07-gaussian-scale-optimization", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v07-gaussian-scale-optimization.md)
for scale profiles, forest-effect plots, and interpretation in map
units.

## See also

[`gaussian_scale_summary`](https://wpeterman.github.io/terradish/reference/gaussian_scale_summary.md),
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)

Other spatial scale:
[`gaussian_scale_summary()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_summary.md),
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`terradish_scale_optim()`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md)

## Examples

``` r
# \donttest{
set.seed(24)
r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
                 ymin = 0, ymax = 10, crs = "EPSG:3857")
terra::values(r) <- rnorm(100)
names(r) <- "x"
graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 98, 5)),
                             directions = 4)
factory <- gaussian_smoothed_loglinear_conductance(graph,
  sigma_lower = 0.4, sigma_upper = 1.5)
response <- simulate_covariance_response(c(x = 0.7, sigma.x = 0.8),
  ~x, graph, conductance_model = factory, tau = 1, sigma = 0.1,
  nu = 10000, seed = 14)$covariance
fit <- terradish(response ~ x, graph, conductance_model = factory,
  measurement_model = wishart_covariance, nu = 10000, optimizer = "newton")
profile <- gaussian_scale_profile(fit, "x", n = 5)
profile$interval
#>     lower     upper 
#> 0.7807026 0.8396359 
confint(profile$fit)
#>             2.5 %    97.5 %
#> x       0.6780223 0.7157837
#> sigma.x 0.7807026 0.8396359
#> attr(,"note")
#> [1] "Stored Gaussian profile intervals are used where available; remaining scale intervals are truncated Wald intervals."
# The stored profile interval replaces the matching scale Wald interval.
# }
```
