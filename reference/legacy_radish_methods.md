# Legacy radish S3 methods

Deprecated S3 methods retained for legacy objects that still dispatch on
the older `"radish"` class names.

## Usage

``` r
# S3 method for class 'radish'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'radish'
summary(object, conf.level = 0.95, ...)

# S3 method for class 'summary.radish'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'radish'
coef(object, ...)

# S3 method for class 'radish'
fitted(object, type = c("response", "distance", "covariance"), ...)

# S3 method for class 'radish'
simulate(object, nsim = 1, seed = NULL, method = "permutation", ...)

# S3 method for class 'radish'
anova(object, ..., alternative = NULL)

# S3 method for class 'radish'
logLik(object, ...)

# S3 method for class 'radish'
AIC(object, ..., k = 2)

# S3 method for class 'radish'
residuals(object, ...)
```

## Arguments

- x:

  A fitted legacy `radish` object.

- digits:

  Number of digits to print.

- ...:

  Additional arguments passed through to generic methods.

- object, alternative:

  Fitted legacy `radish` objects.

- conf.level:

  Confidence level for nuisance-parameter confidence intervals reported
  by
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html).
  These intervals are conditional on the fitted conductance surface.

- signif.stars:

  Should significance stars be printed in coefficient tables?

- type:

  Type of fitted values to return.

- nsim:

  Number of simulated response matrices to generate.

- seed:

  Optional random seed used by
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- method:

  Simulation method.

- k:

  Penalty multiplier supplied to
  [`AIC()`](https://rdrr.io/r/stats/AIC.html).

## Value

- [`print()`](https://rdrr.io/r/base/print.html) returns its input
  invisibly after printing it.

- [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  returns a list of class `summary.terradish` and `summary.radish`,
  containing conductance and nuisance-parameter estimates, uncertainty
  tables, convergence diagnostics, log-likelihood, and model-comparison
  statistics. Conductance coefficients are on the log-conductance scale;
  nuisance-parameter uncertainty is conditional on the fitted
  conductance surface.

- [`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
  vector of fitted conductance coefficients, or an empty vector when
  they are unidentified.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns a
  numeric focal-site matrix of expected responses, resistance distances,
  or graph covariances, selected by `type`.

- [`simulate()`](https://rdrr.io/r/stats/simulate.html) returns a
  numeric response matrix for `nsim = 1`, or a three-dimensional array
  with simulations along the third dimension otherwise.

- [`anova()`](https://rdrr.io/r/stats/anova.html) returns a matrix of
  class `anova`, comparing two nested fits by log-likelihood, parameter
  count, likelihood-ratio statistic, and its reference-distribution
  p-value.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html) returns a numeric
  scalar of class `logLik`, with parameter-count and focal-site-count
  attributes.

- [`AIC()`](https://rdrr.io/r/stats/AIC.html) returns a numeric scalar
  equal to \\-2\ell + kK\\, where \\K\\ is the number of estimated
  parameters. Smaller values favor a fit when response and likelihood
  are comparable.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html) returns a
  numeric matrix of observed minus fitted responses, centered for
  covariance-response models.

See
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
for interpretation and comparison requirements.

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
