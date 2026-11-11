# Print, summarize, compare, and extract from a fitted model

S3 methods for working with fitted objects returned by
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
and summaries returned by
[`summary()`](https://rdrr.io/r/base/summary.html).

## Usage

``` r
# S3 method for class 'terradish'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'terradish'
summary(object, conf.level = 0.95, ...)

# S3 method for class 'summary.terradish'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'terradish'
coef(object, ...)

# S3 method for class 'terradish'
fitted(object, type = c("response", "distance", "covariance"), ...)

# S3 method for class 'terradish'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  method = c("permutation", "parametric"),
  ...
)

# S3 method for class 'terradish'
anova(object, ..., alternative = NULL)

# S3 method for class 'terradish'
logLik(object, ...)

# S3 method for class 'terradish'
AIC(object, ..., k = 2)

# S3 method for class 'terradish'
residuals(object, ...)
```

## Arguments

- x:

  A fitted `terradish` object.

- digits:

  Number of digits to print.

- ...:

  Additional arguments passed through to generic methods.

- object, alternative:

  Fitted `terradish` objects.

- conf.level:

  Confidence level for nuisance-parameter confidence intervals reported
  by
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html).
  These intervals are conditional on the fitted conductance surface.
  [`confint()`](https://rdrr.io/r/stats/confint.html) calls the same
  quantity `level`;
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  and [`plot()`](https://rspatial.github.io/terra/reference/plot.html)
  call it `quantile`.

- signif.stars:

  Should significance stars be printed in coefficient tables?

- type:

  Type of fitted values to return: `"response"` (default) gives the
  expected response matrix on the scale of the data (genetic distance or
  covariance, according to the measurement model), `"distance"` gives
  the fitted resistance distances among focal sites, and `"covariance"`
  gives the fitted graph covariance (the focal-site block of the
  generalized inverse of the graph Laplacian).

- nsim:

  Number of simulated response matrices to generate.

- seed:

  Optional random seed used by
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- method:

  Simulation method. Only `"permutation"` (default) is implemented: it
  adds the residual matrix, with sites permuted jointly in rows and
  columns, to the fitted response. `"parametric"` is reserved and
  currently stops with an error.

- k:

  Penalty multiplier supplied to
  [`AIC()`](https://rdrr.io/r/stats/AIC.html).

## Value

- [`print()`](https://rdrr.io/r/base/print.html) returns its input
  invisibly after printing it.

- [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  returns a list of class `summary.terradish` (also `summary.radish`).
  Its main elements are `ztable` (conductance estimates on the
  log-conductance scale, standard errors, z values, and p values),
  `vcov` and `vcor` (their covariance and correlation matrices), `phi`
  and `phi_table` (nuisance-parameter estimates, standard errors, and
  intervals on the measurement model's scale, conditional on the fitted
  conductance surface), `phi_vcov`, `phi_at_bound` (flags for kernel
  weights at zero), `convergence`, `boundary` and
  `no_structure_boundary`, `loglik`, `df`, `aic`, `gradnorm`, `fcalls`,
  `iters`, `diagnostics`, `call`, and `dim`. `ibe_ratio` holds the
  output of
  [`terradish_ibe_ratio`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md)
  (zero rows without environmental terms). Gaussian scale fits add
  `sigma_table`, with one row per smoothing scale: the estimate and
  standard error in map units, the fitting bounds, and a `near_bound`
  flag for estimates within one percent of a bound. Spline fits add
  `spline_monotonicity`.

- [`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
  vector of fitted conductance coefficients on the log-conductance
  scale, or an empty vector when they are unidentified.

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
  scalar of class `logLik`, with parameter-count (`df`) and
  focal-site-count (`nobs`) attributes.

- [`AIC()`](https://rdrr.io/r/stats/AIC.html) returns \\-2\ell + kK\\,
  where \\K\\ is the number of estimated parameters; the default `k = 2`
  gives Akaike's Information Criterion. Smaller values favor a fit when
  response and likelihood are comparable.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html) returns a
  numeric matrix of observed minus fitted responses, centered on both
  site dimensions for covariance-response models.

## Details

**Reading
[`summary()`](https://rspatial.github.io/terra/reference/summary.html).**
The coefficient table reports each conductance coefficient on the
log-conductance scale with its Wald standard error, z value, and p
value; [`exp()`](https://rdrr.io/r/base/Log.html) of an estimate is the
multiplicative change in conductance per unit of the covariate. Use
[`confint()`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
for intervals on the conductance coefficients. The nuisance table
reports the measurement-model parameters with conditional Wald
intervals.

Parameter names repeat across measurement models with different
meanings. In
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) and
[`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
`tau` is the log residual precision. In
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
and
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
`tau` is the nonnegative multiplier on the graph covariance and `sigma`
is the log nugget variance. See the parameter glossary in
[`terradish-package`](https://wpeterman.github.io/terradish/reference/terradish-package.md).

Nuisance-parameter estimates and intervals from
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
remain on the parameterization used by the measurement model. For
example, MLPE `tau` is log precision and its stored `rho` is
unconstrained; the actual shared-site correlation is
\\\text{plogis}(\rho)/2\\. Consult the selected measurement-model help
page before interpreting `phi`. For a verified numerical zero
residual-variance limit in a Wishart fit,
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
retains the finite log-scale representation of `sigma` but suppresses
its Wald standard error and interval. A zero residual variance does not
imply a zero graph contribution. Covariance-response residuals are
centered on both site dimensions: \\H(S-\Sigma)H\\, where \\H\\ is the
site-centering matrix. Fitted covariance matrices retain their full
\\n\times n\\ form.

[`anova()`](https://rdrr.io/r/stats/anova.html) is a likelihood-ratio
test for nested conductance formulas fitted to the same response, sites,
graph domain, conductance-model factory, and measurement model; Wishart
fits must also use the same `nu`. Both fits must record successful outer
and final nuisance convergence; missing diagnostics require refitting.
The function checks formula-term nesting, but the analyst must still
verify that the full measurement and conductance parameterization is
substantively nested. Its chi-squared reference distribution assumes the
tested parameters are regular interior parameters. It can be unreliable
when a null value is on a parameter boundary or other regularity
conditions fail; use a justified alternative or simulation calibration
in those cases. One boundary case is handled automatically: when the
larger model adds only Wishart environmental kernel weights (`lambda_*`,
which are bounded at zero),
[`anova()`](https://rdrr.io/r/stats/anova.html) warns and reports a
chi-bar-square reference with binomial weights instead of the ordinary
chi-square. The table heading names the reference distribution that was
used.

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
for fitting,
[`conductance`](https://wpeterman.github.io/terradish/reference/conductance.md)
for the fitted surface as a raster,
[`plot.terradish`](https://wpeterman.github.io/terradish/reference/plot.terradish.md)
for the five diagnostic plot types,
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md)
for ranking a model set, and
[`mlpe_response_change`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md)
for putting an MLPE covariate association on the response scale. See
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md)
for a walkthrough of reading these outputs and
[`vignette("v04-model-comparison", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v04-model-comparison.md)
for [`anova()`](https://rdrr.io/r/stats/anova.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) in context.

Other model fitting:
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md),
[`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md),
[`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md),
[`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md)

## Examples

``` r
# \donttest{
data(melip)
covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)
names(covariates) <- c("altitude", "forestcover")
# Coarsened so the example runs quickly; name layers, then scale.
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords[keep])

fit <- terradish(Fst ~ altitude + forestcover, data = surface,
                 conductance_model = loglinear_conductance,
                 measurement_model = mlpe)

fit                    # coefficients, log-likelihood, and AIC
#> Conductance surface estimated by maximum likelihood
#> Call:   terradish(formula = Fst ~ altitude + forestcover, data = surface, 
#>     conductance_model = loglinear_conductance, measurement_model = mlpe)
#> 
#> Coefficients:
#>    altitude  forestcover  
#>    -0.79207     -0.09613  
#> 
#> Loglikelihood: 1495.046 (6 degrees freedom)   AIC: -2978.092 
#> Convergence: converged | projected gradient: 1.93e-07 
summary(fit)           # adds standard errors, Wald tests, and the phi table
#> Conductance surface with 2531 vertices (31 focal) estimated by maximum likelihood
#> Call:   terradish(formula = Fst ~ altitude + forestcover, data = surface, 
#>     conductance_model = loglinear_conductance, measurement_model = mlpe)
#> 
#> Loglikelihood: 1495.046 (6 degrees freedom)
#> AIC: -2978.092 
#> 
#> Number of function calls: 19 
#> Number of optimization steps: 6 
#> Convergence: converged | projected gradient: 1.93e-07 
#> Norm of gradient at MLE: 2.57979e-07 
#> 
#> Nuisance parameters (conditional on fitted conductance surface; measurement-model scale):
#>        Estimate  Std. Error  Lower 95%  Upper 95%
#> alpha  -0.04900   0.01571    -0.07979   -0.01821 
#> beta    0.41808   0.01458     0.38951    0.44665 
#> tau     5.56654   0.23580     5.10439    6.02870 
#> rho     2.06189   0.27586     1.52121    2.60256 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> altitude    -0.79207    0.15007  -5.278 1.31e-07 ***
#> forestcover -0.09613    0.19259  -0.499    0.618    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Correlation of Coefficients:
#>              altitude
#> forestcover 0.7879281

# Conductance coefficients are on the LOG scale, so exponentiate for a
# multiplicative reading: exp(theta) is the change in conductance per
# one-unit (one standard deviation, here) increase in the covariate.
coef(fit)
#>    altitude forestcover 
#> -0.79206795 -0.09612631 
exp(coef(fit))
#>    altitude forestcover 
#>   0.4529072   0.9083493 

logLik(fit)            # carries a df attribute, so AIC()/BIC() work
#> 'log Lik.' 1495.046 (df=6)
AIC(fit)               # compare only fits using the same response and likelihood family
#>       tau 
#> -2978.092 
head(fitted(fit)[, 1:4])     # fitted genetic distances
#>              [,1]         [,2]        [,3]         [,4]
#> [1,] -0.048999521 -0.005447103  0.14508089  0.040612727
#> [2,] -0.005447103 -0.048999521  0.15684956  0.052361591
#> [3,]  0.145080892  0.156849563 -0.04899952  0.091377571
#> [4,]  0.040612727  0.052361591  0.09137757 -0.048999521
#> [5,]  0.062776441  0.074524882  0.10982244  0.005662658
#> [6,] -0.008355229  0.007553196  0.13834571  0.033749003
head(residuals(fit)[, 1:4])  # observed minus fitted
#>              [,1]         [,2]        [,3]         [,4]
#> [1,]  0.048999521  0.014693842 -0.04113520  0.007279262
#> [2,]  0.014693842  0.048999521 -0.04828265 -0.003839283
#> [3,] -0.041135203 -0.048282645  0.04899952 -0.038938063
#> [4,]  0.007279262 -0.003839283 -0.03893806  0.048999521
#> [5,]  0.012324963  0.003186550 -0.04396978  0.032500892
#> [6,]  0.036041844  0.024135824 -0.01236317  0.031867978

# A nested comparison: is the interaction worth its parameter?
fit_int <- terradish(Fst ~ altitude * forestcover, data = surface,
                     conductance_model = loglinear_conductance,
                     measurement_model = mlpe)
anova(fit, fit_int)
#> Likelihood ratio test: chi-square (1 df)
#> Null: ~ altitude + forestcover
#> Alt: ~ altitude + forestcover + altitude:forestcover
#>      logLik Df  ChiSq Df(ChiSq)  Pr(>Chi)    
#> Null 1495.0  6                               
#> Alt  1510.6  7 31.124         1 2.421e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
