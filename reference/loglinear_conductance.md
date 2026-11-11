# Log-linear conductance model

Returns a function of class `"terradish_conductance_model"` that
represents a log-linear mapping from spatial covariates to conductance.
This is the recommended conductance model for most applications.

## Usage

``` r
loglinear_conductance(formula, x)
```

## Arguments

- formula:

  Model formula describing which spatial covariates drive conductance.
  The left-hand side is ignored; only the right-hand side terms are used
  (e.g. `~ altitude + forestcover`).

- x:

  Data frame of spatial covariates extracted from a
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
  object (typically `surface$x`). Every variable named on the right-hand
  side of `formula` must be a column of `x`.

## Value

A function of class `"terradish_conductance_model"` that accepts a
numeric vector of conductance parameters `theta` and returns a list with
elements `conductance` (a vector of per-cell conductance values),
`confint` (a function for confidence intervals), and derivative
functions used internally by the optimizer.

## Details

The conductance at grid cell `i` is:

\$\$C_i = \exp(\theta_1 x\_{i1} + \theta_2 x\_{i2} + \ldots)\$\$

where \\x\_{ij}\\ is the value of covariate \\j\\ at cell \\i\\ and
\\\theta_j\\ is the corresponding conductance parameter.

The intercept is intentionally omitted: multiplying all conductances by
a common factor rescales resistance inversely, and the measurement
model's resistance coefficient absorbs that factor. The absolute
conductance level is therefore not identifiable.

**Interpreting \\\theta\\:**

- \\\theta_j \> 0\\: higher values of covariate \\j\\ increase fitted
  relative conductance and tend to lower graph resistance distance,
  conditional on the other terms and graph domain.

- \\\theta_j \< 0\\: higher values lower fitted relative conductance.

- \\\theta_j = 0\\: the fitted conductance is locally unchanged by that
  term, conditional on the model.

- A one-standard-deviation increase in covariate \\j\\ multiplies fitted
  conductance by \\\exp(\theta_j)\\ when the covariate was standardized
  and the term is linear.

These are model-based associations, not direct estimates of movement
rate, habitat suitability, or a causal landscape effect.

The exponential link guarantees strictly positive conductances for any
real \\\theta\\.

Categorical covariates must be stored as `factor` columns in `x` (see
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
for how to encode them). They are dummy-coded using the default R
contrasts via
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html), with one
level dropped as a reference. In-formula transformations such as
`I(x^2)` and interaction terms `x * z` are supported.

## See also

[`smooth_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md),
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other conductance models:
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`smooth_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md),
[`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md)

## Examples

``` r
x <- data.frame(altitude = c(-1, 0, 1), forestcover = c(0.2, 0.6, 0.5))
model <- loglinear_conductance(~ altitude + forestcover, x)

# Evaluate at specific parameter values
fit <- model(c(altitude = 0.3, forestcover = -0.2))
fit$conductance  # per-cell conductance values
#> [1] 0.7117703 0.8869204 1.2214028

# Interaction and polynomial terms work too
x2 <- data.frame(altitude = c(-1, -0.25, 0.5, 1),
                 fc = c(0.2, 0.6, 0.5, 0.8))
model2 <- loglinear_conductance(~ altitude + I(altitude^2) + fc, x2)
```
