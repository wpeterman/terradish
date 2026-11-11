# Log-linear conductance model with spline terms

Returns a function of class `"terradish_conductance_model"` that
represents a log-linear conductance surface after expanding selected
numeric covariates into spline basis columns. The result is an
unpenalized regression spline with a fixed number of coefficients,
fitted through the selected terradish measurement likelihood.

## Usage

``` r
smooth_loglinear_conductance(
  formula,
  x,
  df = 4L,
  basis = c("ns", "bs"),
  degree = 3L,
  intercept = FALSE
)
```

## Arguments

- formula:

  Model formula describing which spatial covariates drive conductance.
  Smooth terms are written as `s(variable, df = 4)` or
  `s(variable, k = 4)`. Ordinary linear terms can be mixed with smooth
  terms, e.g. `~ forestcover + s(altitude, df = 4)`.

- x:

  Data frame of spatial covariates extracted from a
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
  object (typically `surface$x`).

- df:

  Default degrees of freedom for smooth terms that do not supply `df` or
  `k`.

- basis:

  Default spline basis. `"ns"` uses
  [`ns`](https://rdrr.io/r/splines/ns.html) natural splines; `"bs"` uses
  [`bs`](https://rdrr.io/r/splines/bs.html) B-splines.

- degree:

  B-spline polynomial degree used when `basis = "bs"`.

- intercept:

  Logical; include intercept columns within each spline basis? The
  default is `FALSE`, matching terradish's usual convention of omitting
  a global conductance intercept.

## Value

A function of class `"terradish_conductance_model"` that accepts a
numeric vector of conductance parameters `theta` and returns conductance
values plus derivative functions used internally by the optimizer.
Derivatives are with respect to the expanded basis columns.

## Details

The model first expands all `s()` terms into basis columns and then
applies the same positive log-link used by
[`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md):

\$\$C_i = \exp(B_i \theta)\$\$

Each spline column is centered over the active graph cells. The fitted
knots, boundary knots, degree, and column centers are retained and
reused for plotting, focal-support prediction, and coarse warm starts.
Centering changes the arbitrary conductance reference level but
preserves the likelihood and fitted conductance coefficients.
Uncertainty bands are relative to the landscape-mean spline contribution
to log conductance.

where \\B_i\\ is the row of the expanded spline/parametric design matrix
for grid cell \\i\\. The same conductance factory works with MLPE and
Wishart measurement models.

Smoothing parameters are not estimated in this first implementation. The
effective smoothness is controlled by fixed basis dimension through `df`
or `k`; use small values for stable first-pass model comparison. AIC and
AICc at nominal Wishart information can over-select flexible curves.
Prefer spatial CV for selecting conductance terms, and examine
sensitivity to `nu` before interpreting uncertainty. Tails beyond the
sampled covariate range are poorly identified. A selected curve can
describe how a process maps onto the graph without identifying an
ecological response mechanism.

`summary(fit)$spline_monotonicity` describes each fitted smooth over the
range at focal sites. It reports whether the curve is monotone, its
direction, and the number of derivative sign changes. This is a shape
diagnostic, not an uncertainty test.

## Worked example

See
[`vignette("v06-spline-conductance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v06-spline-conductance.md)
for fitted conductance curves, genetic-response plots, and model
comparison.

## See also

[`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md),
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other conductance models:
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md),
[`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md)

## Examples

``` r
x <- data.frame(altitude = seq(-1, 1, length.out = 8),
                forestcover = c(0.2, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2, 0.1))
model <- smooth_loglinear_conductance(~ forestcover + s(altitude, df = 3), x)
theta <- attr(model, "default")
fit <- model(theta)
fit$conductance
#> [1] 1 1 1 1 1 1 1 1
```
