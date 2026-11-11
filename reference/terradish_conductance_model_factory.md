# Conductance models: the factory interface

A conductance model maps the spatial covariates at each raster cell to a
positive conductance. Conductance models are created by factories, which
are functions of class `"terradish_conductance_model_factory"`. Pass a
factory to the `conductance_model` argument of
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Details

**Choosing a model.**

- [`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md):

  Log-linear conductance in the formula terms, including factors,
  interactions, and polynomial terms. The default.

- [`smooth_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md):

  Log-linear conductance with fixed-degree regression splines for `s()`
  terms.

- [`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md):

  A constructor that returns a factory. Selected rasters are Gaussian
  smoothed and each smoothing scale is estimated with the coefficients.

**Interface.** A factory has the signature `factory(formula, x)`, where
`formula` is the conductance formula (its left-hand side is ignored) and
`x` is the data frame of covariates stored in a
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
object (`surface$x`).
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
calls the factory for you. It returns a function of class
`"terradish_conductance_model"` that takes the conductance parameter
vector `theta` and returns a list with `conductance` (one positive value
per graph vertex), `confint` (a function giving pointwise intervals from
a parameter covariance matrix), and the derivative functions
`df__dtheta`, `df__dx`, `d2f__dtheta_dtheta`, and `d2f__dtheta_dx` that
the optimizer uses. The returned function carries a `"default"`
attribute with named starting values, and may carry `"lower"` and
`"upper"` bounds. The legacy class names
`"radish_conductance_model_factory"` and `"radish_conductance_model"`
are also carried for compatibility.

All built-in models use a log link, so every coefficient is a change in
log conductance and [`exp()`](https://rdrr.io/r/base/Log.html) of a
coefficient is a conductance multiplier. No intercept is estimated,
because a common conductance multiplier is absorbed by the measurement
model.

## Worked examples

See
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md)
for a log-linear model,
[`vignette("v05-categorical-conductance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v05-categorical-conductance.md)
for class contrasts,
[`vignette("v06-spline-conductance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v06-spline-conductance.md)
for nonlinear curves and their plots, and
[`vignette("v07-gaussian-scale-optimization", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v07-gaussian-scale-optimization.md)
for estimated smoothing scales.

## See also

[`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md),
[`smooth_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md),
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other conductance models:
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md),
[`smooth_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md)
