# Rescale the effective information in a Wishart fit

Changes `nu` without refitting the coefficients. Wishart objectives and
their derivatives are proportional to this user-supplied quantity.

## Usage

``` r
terradish_rescale_nu(fit, nu)
```

## Arguments

- fit:

  A full or slim fitted Wishart model.

- nu:

  New finite, positive value of `nu`.

## Value

A fitted object with rescaled likelihood, curvature, AIC, and comparison
metadata. Estimates are unchanged; `nu_original` records the original
value across repeated rescalings. Standard errors computed by
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
and [`vcov()`](https://rdrr.io/r/stats/vcov.html) use the rescaled
curvature.

## Details

This is a sensitivity calculation, not an estimate of effective
information. It does not remedy misspecification or dependent histories.
If `nu` is divided by 20, standard errors grow by `sqrt(20)`; the
likelihood and likelihood-ratio statistics are divided by 20. AIC
changes because its likelihood term changes while its parameter penalty
stays fixed. Report conclusions over a scientifically plausible range.

## See also

[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)

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
sensitivity <- terradish_rescale_nu(fit, nu = 50)
rbind(original = coef(fit), sensitivity = coef(sensitivity))
#>                     x
#> original    0.6761321
#> sensitivity 0.6761321
rbind(original = sqrt(diag(vcov(fit))),
      sensitivity = sqrt(diag(vcov(sensitivity))))
#>                      x
#> original    0.08944225
#> sensitivity 0.39999790
```
