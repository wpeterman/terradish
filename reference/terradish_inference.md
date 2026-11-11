# Standard errors and confidence intervals for conductance parameters

[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) report the
uncertainty of the fitted conductance parameters from the curvature
(observed information) of the log-likelihood at the estimate. They work
on a full or slim fit.

## Usage

``` r
# S3 method for class 'terradish'
vcov(object, ...)

# S3 method for class 'terradish'
nobs(object, ...)

# S3 method for class 'terradish'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fitted `terradish` object.

- ...:

  Unused arguments for generic compatibility.

- parm:

  Names or indices of parameters to include.

- level:

  Confidence level between zero and one. This is the quantity that
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  calls `conf.level` and that
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  and [`plot()`](https://rspatial.github.io/terra/reference/plot.html)
  call `quantile`.

## Value

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the named
covariance matrix of the conductance parameters, the inverse of the
observed information, on the scale of
[`coef()`](https://rdrr.io/r/stats/coef.html).
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per parameter and two columns holding the lower and upper
Wald limits, labeled by their percentiles.
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
focal sites.

## Details

**Reading the output.** Both functions are on the scale of
[`coef()`](https://rdrr.io/r/stats/coef.html): log conductance per unit
of each covariate (per standard deviation for covariates scaled with
[`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)),
and map units for Gaussian `sigma.<layer>` smoothing scales.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns Wald
intervals, the estimate plus or minus a normal quantile times its
standard error. `sqrt(diag(vcov(fit)))` gives the standard errors that
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
prints. Exponentiate a log-conductance estimate or interval endpoint to
read it as a conductance multiplier: `exp(confint(fit))` is the interval
for the ratio of conductance per unit increase in the covariate, and an
interval that excludes 1 on that scale excludes zero on the log scale.
Do not exponentiate rows for Gaussian smoothing scales.

The intervals are conditional on the fitted model. For Wishart fits they
also depend on the supplied `nu`: standard errors scale as
`1 / sqrt(nu)` (see
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)).
Uncertainty in the nuisance parameters is reported separately by
[`summary()`](https://rspatial.github.io/terra/reference/summary.html).

Gaussian scale intervals are Wald intervals truncated at the fitted
bounds. Stored profile intervals take precedence when their confidence
level matches. The observation count is the number of focal sites
(individuals or populations), the sample size that AICc and BIC use in
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md).
Pairwise observations are not counted, because pairs that share a site
are not independent.

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`slim_terradish`](https://wpeterman.github.io/terradish/reference/slim_terradish.md),
[`gaussian_scale_profile`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md),
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
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
coef(fit)       # Relative log conductance per unit of x.
#>         x 
#> 0.6761321 
vcov(fit)       # Joint uncertainty in conductance parameters.
#>             x
#> x 0.007999916
confint(fit)    # A model-based interval conditional on the supplied nu.
#>       2.5 %    97.5 %
#> x 0.5008285 0.8514357
exp(confint(fit))  # The same interval as a conductance multiplier.
#>      2.5 %   97.5 %
#> x 1.650088 2.343008
nobs(fit)       # Number of focal sites, the sample size for AICc and BIC.
#> [1] 6
```
