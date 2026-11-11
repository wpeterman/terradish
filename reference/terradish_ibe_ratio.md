# Express environmental effects in resistance-distance units

Reports the resistance-distance equivalent of one unit of environmental
difference, using the environmental coefficient divided by the
resistance coefficient in the additive measurement model.

## Usage

``` r
terradish_ibe_ratio(fit)
```

## Arguments

- fit:

  A fitted model from
  [`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
  or
  [`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
  or a fit with no environmental terms.

## Value

A data frame with covariate, ratio, SE, lower, upper, and family.
Intervals use a 95 percent delta-method approximation based on the joint
covariance of both coefficients, including fitted-conductance
uncertainty. Fits without environmental terms return a zero-row data
frame.

## Details

The ratios are \\\gamma_k/\beta\\ for MLPE and \\\lambda_k/\tau\\ for
Wishart models. They require the same environmental transform and
scaling for comparison. Their interpretation agrees exactly for the same
covariance-derived distance and approximately for distances that scale
linearly with it. Separate fits estimate their own conductance surfaces,
so their resistance distances can differ. Absolute measurement
coefficients are not comparable across families.

A zero resistance coefficient gives undefined ratios and uncertainty,
returned as NA. A zero Wishart kernel coefficient receives a one-sided
upper limit and no symmetric lower limit. Boundary inference remains
approximate. The distance-equivalent ratio interpretation also appears
in BEDASSLE, although that model uses geographic rather than resistance
distance.

## References

Bradburd GS, Ralph PL, Coop GM (2013). Disentangling the effects of
geographic and ecological isolation on genetic differentiation.
Evolution, 67, 3258–3273.
[doi:10.1111/evo.12193](https://doi.org/10.1111/evo.12193) .

## See also

[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

## Examples

``` r
# A small model-based example, with fixed information supplied explicitly.
r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
names(r) <- "x"
sites <- terra::xyFromCell(r, c(1, 8, 20, 37, 57, 64))
graph <- conductance_surface(r, sites)
simulated <- simulate_covariance_response(c(x = 0.6), ~x, graph,
                                          nu = 1000, seed = 93)
response <- simulated$covariance
fit <- terradish(response ~ x, graph, measurement_model = wishart_covariance,
                 nu = 1000)
environment <- pairwise_covariates(
  pairwise_endpoint_covariates(data.frame(climate = c(2, 4, 1, 7, 5, 9))))
measurement <- wishart_covariates(environment, model = "wishart_covariance")
kernel <- attr(measurement, "kernel_covariates")[, , 1]
response <- response + 0.1 * kernel
joint <- terradish(response ~ x, graph, measurement_model = measurement,
                   nu = 1000)
terradish_ibe_ratio(joint)
#>                covariate      ratio        SE      lower     upper
#> 1 lambda_absdiff_climate 0.09621109 0.0126249 0.07146673 0.1209554
#>               family
#> 1 wishart_covariance
# Each ratio is the resistance-distance equivalent of one unit of difference.
```
