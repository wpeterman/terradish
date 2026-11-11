# Wishart measurement model with environmental covariance kernels

Creates a Wishart measurement model whose fitted covariance includes the
resistance-implied covariance plus one or more positive semidefinite
kernels built from site-level environmental covariates. Use this as the
Wishart analogue of
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
when your measurement model is
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
(squared-distance responses) or
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
(covariance-matrix responses).

## Usage

``` r
wishart_covariates(
  x,
  coords = NULL,
  transform = c("absdiff", "sqdiff", "euclidean", "manhattan"),
  scale = FALSE,
  model = c("wishart_covariance", "generalized_wishart")
)
```

## Arguments

- x:

  Site-level covariates. Supported inputs are the same as
  [`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md):
  a numeric vector, matrix, data frame, or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).
  A previously constructed `terradish_pairwise_covariates` object is
  accepted unchanged.

- coords:

  Required when `x` is a raster. Focal-point coordinates in the same
  projection as `x`; accepts the same inputs as
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- transform:

  Pairwise environmental difference: `"absdiff"` (default), `"sqdiff"`,
  `"euclidean"`, or `"manhattan"`. See
  [`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md).
  Ignored for an existing pairwise-covariate object.

- scale:

  Logical. Standardize site-level covariates to zero mean and unit
  variance before constructing kernels? Recommended when covariates use
  different units or scales.

- model:

  Which Wishart likelihood to use. `"generalized_wishart"` is
  appropriate only when `S` is an admissible squared-distance matrix
  coherently related to a centered positive-semidefinite covariance
  matrix. An arbitrary distance matrix, including an unchecked
  \\F\_{ST}\\ matrix, is not sufficient. `"wishart_covariance"` is
  appropriate when `S` is a **covariance** matrix (e.g. from
  [`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html)).
  In both cases the same \\\Sigma\\ parameterization is used; they
  differ in how \\\Sigma\\ is compared to the observed data.

## Value

A function of class `"terradish_measurement_model"` suitable for the
`measurement_model` argument of
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
and
[`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md).
The nuisance parameter vector \\\phi\\ estimated by
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
contains `tau`, one `lambda_<covariate>` per environmental kernel, and
`sigma` (in that order). The function stores the kernel array in
attribute `"kernel_covariates"` and supports site-subsetting for model
fitting through
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Details

A zero kernel coefficient sets the boundary flag while retaining
conductance inference when the resistance coefficient is positive.
Summaries mark zero kernel coefficients and give one-sided upper limits
instead of symmetric intervals. Tests adding only kernel coefficients
use a chi-bar-square reference with binomial weights. This reference is
exact for information-orthogonal kernels and approximate otherwise. MLPE
environmental coefficients are unconstrained and retain ordinary
chi-square tests.

`wishart_covariates()` is the Wishart analogue of
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
but environmental covariates enter the model through positive
semidefinite kernels. For each pairwise difference matrix \\D_k\\, the
kernel is \\K_k=-\frac12 H D_k H\\, where \\H\\ centers sites. A kernel
with a negative eigenvalue beyond numerical tolerance is rejected. The
`"sqdiff"` transform gives the centered outer product of the original
covariate. The fitted covariance is:

\$\$\Sigma = \tau E(\theta) + \sum_k \lambda_k K_k + \exp(\sigma) I\$\$

where \\E(\theta)\\ is the resistance-implied covariance (generalized
inverse of the graph Laplacian at conductance parameters \\\theta\\),
\\\tau\\ is the weight on conductance-implied covariance, \\\lambda_k
\geq 0\\ are nonnegative kernel weights, and \\\exp(\sigma)\\ is
diagonal residual variance not represented by the graph or environmental
kernels.

The nuisance parameters estimated alongside \\\theta\\ are:

- `tau`:

  Nonnegative scale on the resistance-implied covariance \\E\\. A value
  near zero indicates no detectable contribution from
  conductance-implied covariance under this model.

- `lambda_<covariate>`:

  Nonnegative weight per unit of the corresponding pairwise
  environmental difference. A negative association is represented at the
  zero bound; MLPE environmental slopes are unconstrained. The weight is
  not a causal IBE effect.

- `sigma`:

  Log-scale nugget: identity component added as \\\exp(\sigma) I\\,
  representing diagonal variance not captured by the other covariance
  components. It is not a smoothing scale or direct drift or
  population-size estimate.

This design naturally preserves positive definiteness of \\\Sigma\\ when
the nugget is positive and all kernel weights are nonnegative. It is
evaluated on site contrasts, as in
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).
The implied squared-distance mean is \\\tau R\_{ij}+\sum_k\lambda_k
D\_{k,ij} +2\exp(\sigma)\\, using the same additive environmental terms
as MLPE. The likelihoods and coefficient constraints still differ. Use
`wishart_covariates()` when the response satisfies the selected Wishart
model's requirements and you can justify an effective degrees-of-freedom
value `nu`; use
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
when `nu` is unknown or when the distance-regression framework better
matches the response. Neither formulation cleanly separates causal IBE
from IBR when predictors and kernels are spatially correlated. Scaling
is applied to site values before the pairwise transform; kernels are not
normalized afterward. Site subsets retain these original pairwise values
and rebuild centered kernels without restandardizing the subset. In an
unpublished development null simulation, AIC selected an environmental
kernel in 19 of 40 runs at the SNP-count value of `nu`, and in none at a
replicate-calibrated value. These are results of that design, not
general error-rate guarantees. Use
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
for sensitivity and `terradish_cv_folds(..., nuisance = "fixed")` to
assess predictive support for a measurement-model extension.

## Numerical nugget limits

Fitting uses the same small-nugget safeguards as
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
retaining supplied bounds and fixed parameters. A verified numerical
zero residual-variance limit does not remove the graph or environmental
kernel coefficients.
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
omits Wald uncertainty for the log-scale nugget at that limit.

## Worked example

See
[`vignette("v08-ibe-ibr-workflow", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v08-ibe-ibr-workflow.md)
for environmental covariance kernels and their distinction from MLPE
slopes.

## Interpreting pairwise effects

Construct one
[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
object to use the same environmental transforms and geographic distances
in both likelihood families. Their fitted conductance surfaces are
conditional on these terms, and conductance signs can change when the
mean structure changes. Report fits with and without pairwise effects. A
uniform `~ 1` conductance model with pairwise terms represents IBD plus
IBE, not IBE alone.

A positive environmental coefficient means differentiation increases
with environmental difference conditional on resistance. It does not
distinguish assortative mating, dispersal filtering, local adaptation,
or site variance that rises at environmental extremes. MLPE coefficients
are unconstrained; Wishart kernel weights are nonnegative. Compare their
effects in resistance-distance units with
[`terradish_ibe_ratio`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
which divides by the resistance coefficient and propagates joint
uncertainty. Ratios are undefined when that denominator is zero.

## Effective Wishart information

For Wishart likelihoods, `nu` scales the log likelihood. Point estimates
do not depend on it, but standard errors, likelihood-ratio tests, and
information criteria do. It is a dispersion parameter that you set,
rather than the number of markers. Linkage and the population history
shared by loci can make effective information much smaller than the SNP
count. In unpublished forward-time simulations run during development,
standard errors at the SNP count were three to five times smaller than
the spread among replicate histories. Those model-based errors
understated uncertainty.

Use
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
to compare conductance formulas, and
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
to report inference across plausible values. For the same successful
folds and nuisance settings, changing `nu` rescales Wishart CV scores
without changing their ranking. Neither procedure identifies the true
effective information or repairs an unsuitable response matrix.

## References

McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica
Sinica 19:631-649.

## See also

[`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md),
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)

## Examples

``` r
library(terra)

data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)
names(covariates) <- c("altitude", "forestcover")
# Coarsened so the example runs quickly; name layers, then scale.
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
coords <- coords[keep]
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords)

# Build an altitude kernel for a covariance-Wishart model.
g_wc <- wishart_covariates(terra::unwrap(melip.altitude), coords = coords,
                           model = "wishart_covariance", scale = TRUE)
inherits(g_wc, "terradish_measurement_model")  # TRUE
#> [1] TRUE

# Simulate an admissible covariance response, then fit the extended model.
# The simulator takes the nugget variance on its natural scale (0.2); the
# fitted `sigma` is its logarithm.
# \donttest{
simulated <- simulate_covariance_response(
  theta = c(altitude = 0.2, forestcover = -0.1),
  formula = ~ altitude + forestcover,
  data = surface, tau = 0.8, sigma = 0.2, nu = 40, seed = 1
)
fit_joint_w <- terradish(
  simulated$covariance ~ altitude + forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_wc,
  nu                = 40
)
# phi table shows tau, lambda_altitude, and sigma (log nugget).
summary(fit_joint_w)
#> Conductance surface with 2531 vertices (31 focal) estimated by maximum likelihood
#> Call:   terradish(formula = simulated$covariance ~ altitude + forestcover, 
#>     data = surface, conductance_model = loglinear_conductance, 
#>     measurement_model = g_wc, nu = 40)
#> 
#> Loglikelihood: 13.50245 (5 degrees freedom)
#> AIC: -17.00491 
#> 
#> Number of function calls: 12 
#> Number of optimization steps: 4 
#> Convergence: converged | projected gradient: 7.51e-13 
#> Norm of gradient at MLE: 9.585887e-13 
#> 
#> Nuisance parameters (conditional on fitted conductance surface; measurement-model scale):
#>                                         Estimate   Std. Error  Lower 95%
#> tau                                      0.990540   0.166938    0.663347
#> lambda_absdiff_spatial_template_lowres   0.017644   0.014093   -0.009977
#> sigma                                   -1.625156   0.118411   -1.857237
#>                                         Upper 95%
#> tau                                      1.317733
#> lambda_absdiff_spatial_template_lowres   0.045265
#> sigma                                   -1.393075
#> 
#> Resistance-distance equivalent of one unit of environmental difference:
#>                               covariate     ratio         SE       lower
#>  lambda_absdiff_spatial_template_lowres 0.0178124 0.01543092 -0.01243165
#>       upper             family
#>  0.04805646 wishart_covariance
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> altitude    -0.04533    0.12550  -0.361    0.718
#> forestcover -0.11802    0.19563  -0.603    0.546
#> 
#> Correlation of Coefficients:
#>              altitude
#> forestcover 0.4338603
# }
```
