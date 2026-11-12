# Fit a conductance model by maximum likelihood

Uses maximum likelihood to fit a parameterized conductance surface to
genetic data, given a function relating spatial data to conductance (a
"conductance model") and a function linking the graph kernel to an
observed genetic distance or covariance matrix (a "measurement model").

## Usage

``` r
terradish(
  formula,
  data,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  nu = NULL,
  theta = NULL,
  leverage = FALSE,
  nonnegative = TRUE,
  conductance = TRUE,
  optimizer = c("newton", "bfgs", "auto"),
  verbose = FALSE,
  control = NewtonRaphsonControl(ctol = 1e-06, ftol = 1e-06),
  validate = FALSE,
  cores = 1L,
  curvature = c("exact", "gauss_newton"),
  solver = c("direct", "auto", "amg"),
  solver_control = NULL,
  approximation = c("none", "landmark", "coarse_raster"),
  approximation_control = NULL,
  measurement_control = NULL,
  slim = FALSE
)
```

## Arguments

- formula:

  A formula with an observed genetic distance or covariance matrix on
  the left-hand side and conductance covariates from `data` on the
  right-hand side.

- data:

  An object of class `terradish_graph` (see
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md))

- conductance_model:

  A function of class `terradish_conductance_model_factory` (see
  [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md))

- measurement_model:

  A function of class `terradish_measurement_model` (see
  [`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md))

- nu:

  Effective Wishart degrees of freedom passed to measurement models that
  require it, such as
  [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
  and
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).
  It must be supplied and is not estimated. It controls effective
  information, not a marker count. Report sensitivity across plausible
  values with
  [`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md).
  Ignored by non-Wishart measurement models.

- theta:

  Optional starting values for the conductance parameters, on the scale
  reported by [`coef()`](https://rdrr.io/r/stats/coef.html) (log
  conductance per unit of each covariate; map units for Gaussian
  `sigma.<layer>` terms) and in the order of the model terms. The
  default `NULL` uses the conductance model's default start, which is
  zero for every log-linear coefficient.

- leverage:

  Logical. Compute leverage (influence) arrays? When `TRUE`, the fit
  stores `leverage$S`, a focal-site by focal-site by parameter array
  giving the first-order change in each conductance estimate per unit
  change in each entry of the response matrix, and `leverage$X`, the
  corresponding array for the spatial covariates (the change in each
  estimate per unit change in a covariate value at a graph vertex).
  Large absolute values flag observations or cells that drive an
  estimate. Unavailable for Gaussian scale-aware models.

- nonnegative:

  Force regression-like `measurement_model` to have nonnegative slope?

- conductance:

  Retained for backward compatibility. Only `conductance = TRUE` is
  currently implemented.

- optimizer:

  The optimization algorithm to use: `newton` uses the exact Hessian
  during optimization, with computational cost that grows linearly with
  the number of parameters; `bfgs` uses an approximate Hessian during
  optimization with cheaper iterations but often more steps; and `auto`
  selects `bfgs` when there are more than three conductance parameters
  and `newton` otherwise, unless the conductance model declares a
  preferred optimizer (the Gaussian scale-aware model always prefers
  `bfgs`). Regardless of the optimizer, `terradish()` evaluates the
  fitted model once more at the final parameter values using the
  requested `curvature`, so standard errors and summaries use the
  selected final derivative machinery.

- verbose:

  Logical. Should the optimizer report its progress? The default `FALSE`
  fits silently. Set `TRUE` to see the iteration count, objective value,
  gradient norm, and Hessian determinant at each step, which is what you
  want when a fit fails to converge or you are tuning `control`. The
  trace is emitted with
  [`message()`](https://rdrr.io/r/base/message.html), so you can also
  capture or silence it with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html). Setting
  `verbose` explicitly overrides the `verbose` field of `control`;
  leaving it alone lets a hand-built `control` object speak for itself.
  The much noisier line-search trace stays separate: turn it on with
  `control = NewtonRaphsonControl(ls.control = HagerZhangControl(verbose = TRUE))`.

- control:

  A list containing options for the optimization routine (see
  [`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  for list)

- validate:

  Deprecated and ignored. Numerical validation of leverage is not
  implemented in `terradish()`. Any value other than `FALSE` gives a
  deprecation warning and the fit proceeds unchanged. The argument will
  be removed in a future release.

- cores:

  Number of worker processes. `1` evaluates serially. Parallel execution
  is experimental and uses PSOCK workers on every operating system. Only
  Hessian and partial solves run in parallel; AMG and cached CHOLMOD
  ignore this setting. Multiple workers can be slower on small graphs. A
  one-time message explains these limits when applicable. On Windows,
  one PSOCK cluster is retained for the duration of the fit instead of
  being rebuilt for every derivative evaluation.

- curvature:

  Curvature used for optimization steps and for the returned covariance
  matrix. `"exact"` (default) uses the exact Hessian. `"gauss_newton"`
  uses an information-based approximation that drops residual-weighted
  second-derivative terms and requires only first derivatives of the
  conductance model. It equals the exact Hessian only when those omitted
  terms are zero or negligible; convergence alone does not guarantee
  agreement. For
  [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) the
  two can differ materially: in the package's tests on the `melip` data,
  Gauss-Newton standard errors differed from the exact ones by 3 to 14
  percent in three configurations and by about a factor of two in a
  fourth, so report `"exact"` standard errors for MLPE fits. Standard
  errors from
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  are then the asymptotic information-based errors. With
  `leverage = TRUE` the leverage diagnostics inherit the same
  approximation.

- solver:

  Linear-system solver. `"direct"` uses sparse Cholesky, `"amg"` uses
  algebraic multigrid, and `"auto"` chooses AMG above the large-graph
  threshold regardless of the number of right-hand sides. AMG iteration
  counts can rise with conductance contrast.

- solver_control:

  Optional named list of solver settings passed to
  [`terradish_algorithm`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md).
  For `solver = "direct"`, supported entries include `factorization`,
  `supernodal_min_vertices`, `supernodal_max_rhs`, `perm`, and
  `solve_backend`; the latter can be `"matrix"` or the experimental
  `"cholmod_cpp"` and `"cholmod_cpp_cached"` backends.
  `reference = "auto"` grounds a maximum-conductance vertex on
  high-contrast graphs; `"last"` retains the original reference. This
  changes the numerical solve while preserving the same graph model.
  Direct solves also reject relative system residuals exceeding
  `residual_tolerance = 1e-8`; see
  [`terradish_algorithm`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md).
  For `solver = "auto"`, supported selection entries include
  `auto_direct_max_vertices`, `auto_amg_min_vertices`. The legacy
  `auto_direct_max_rhs` setting is accepted but no longer affects
  selection. For `solver = "amg"` or `"auto"`, `terradish()` also
  understands an adaptive schedule with entries such as `adaptive`,
  `tol_early`, `tol_mid`, `tol_final`, `maxit_early`, `maxit_mid`,
  `maxit_final`, and `warmup_evals`. AMG controls also support
  `reuse_preconditioner` to reuse multigrid hierarchy information across
  nearby optimization steps when possible, and
  `reuse_preconditioner_max_age` to periodically rebuild a reused
  hierarchy. Direct supernodal factorizations can benefit from a
  threaded BLAS, but those thread counts are controlled by the external
  R/BLAS build rather than by `terradish()` itself. If a fit stops with
  a CHOLMOD or non-positive-definite reduced Laplacian message, the
  current conductance values made the graph numerically singular. Common
  fixes are to scale raster covariates, start conductance parameters
  closer to zero, check retained raster cells for missing or infinite
  values, and retry with `solver = "auto"` or `solver = "amg"` while
  diagnosing the problem.

- approximation:

  Exploratory approximation used during optimization. `"none"` uses the
  full focal set throughout. `"landmark"` optimizes first on a
  space-filling subset of focal populations and then refines on the full
  likelihood when supported by the measurement model. `"coarse_raster"`
  optimizes first on an aggregated raster and then refines on the
  full-resolution graph. If multiple coarse factors are supplied, they
  are evaluated from coarsest to finest before the final full-resolution
  stage. This is an opt-in warm-start strategy, not a replacement for
  the exact full-resolution likelihood. Exact refinement is required.

- approximation_control:

  Optional named list controlling the landmark or coarse-raster
  approximation. Landmark entries include `n_landmarks`, `fraction`,
  `min_landmarks`, `max_landmarks`, `method` (`"spacefill"`, `"random"`,
  or `"sequential"`), `seed`, and `exact_refine`. For
  `approximation = "coarse_raster"`, supported entries include `factor`,
  `aggregate_fun`, `directions`, `exact_refine`, and `refine_control`;
  this path requires `conductance_surface(..., saveStack = TRUE)`.
  `factor` may be a single positive integer or a vector such as
  `c(8, 4, 2)` for a multilevel coarse-to-fine warm start. The final
  reported fit is still evaluated on the full data. `refine_control` may
  be a
  [`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  list used only for that final exact full-resolution refinement; by
  default the main `control` is reused. For
  `measurement_model = leastsquares`, landmark exact-refinement is
  guarded and the approximation stage is skipped.

- measurement_control:

  Optional
  [`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  object for the inner nuisance-parameter profile. This is separate from
  `control`, which governs optimization of the conductance parameters.
  The fitted object records the inner optimizer's convergence code and
  iteration count for the final evaluation. A failed nuisance warm start
  is retried from the model's default start, retaining any fixed values
  and custom bounds. Warnings from the failed attempt remain visible;
  the retry is accepted only if it converges. Standard MLPE nuisance
  fits scale the resistance slope internally; their convergence
  criterion is therefore independent of resistance units.

- slim:

  Logical. If `TRUE`, apply
  [`slim_terradish`](https://wpeterman.github.io/terradish/reference/slim_terradish.md)
  to the returned fit, removing model closures and leverage arrays. Use
  this for durable storage only when later conductance prediction and
  model reevaluation are unnecessary.

## Value

An object of class `c("terradish", "radish")`, a list. Most users need
only the methods:
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
for estimates, standard errors, and nuisance-parameter intervals,
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) for conductance
coefficients,
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
and [`plot()`](https://rspatial.github.io/terra/reference/plot.html) for
maps and curves, and [`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`anova()`](https://rdrr.io/r/stats/anova.html), or
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md)
for comparison (see
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)).
The list elements are:

- call, formula:

  The matched call and the model formula.

- dim:

  Named vector with the number of graph vertices (raster cells), focal
  sites, and graph edges.

- convergence:

  List describing the conductance optimization: `code` (0 converged, 1
  iteration limit reached, 2 stalled or line search failed), `message`,
  `iterations`, `max_abs_projected_gradient`, `criterion`, `boundary`
  (`TRUE` when an estimate sits on a bound or no resistance structure
  was detected), and `restarted`. Check `object$convergence$code == 0`
  before interpreting estimates.

- fit:

  The final graph and likelihood evaluation. `object$fit$phi` holds the
  nuisance (measurement-model) estimates on the measurement model's
  internal scale, for example log precision;
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  reports them with intervals. `object$fit$fitted` is the expected
  response matrix among focal sites. `object$fit$subproblem$convergence`
  is the nuisance optimizer's code (0 means converged), which should
  also be 0 before estimates are interpreted.

- loglik, df, aic:

  Maximized log-likelihood, number of estimated parameters (conductance
  coefficients plus nuisance parameters), and \\AIC = -2\ell + 2K\\.

- mle:

  Conductance estimates (`theta`, on the log-conductance scale per unit
  of each covariate as supplied to the graph), and the gradient and
  Hessian of the log-likelihood at the estimate. Elements ending in
  `_internal` are on the optimizer's scale.

- submodels:

  The conductance-model and measurement-model functions used for the
  fit;
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  and [`plot()`](https://rspatial.github.io/terra/reference/plot.html)
  reuse them.

- comparison:

  Measurement family, `nu`, graph fingerprint, and response columns that
  [`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md)
  and [`anova()`](https://rdrr.io/r/stats/anova.html) use to check that
  fits are comparable.

- cost, diagnostics, approximation:

  Newton steps and function calls, solver timing and evaluation
  diagnostics, and a record of any landmark or coarse-raster warm start.

- gaussian_scale_info, spline_monotonicity:

  Metadata for Gaussian smoothing and spline conductance terms, or
  `NULL`.

- leverage:

  Leverage arrays for the response (`S`) and spatial covariates (`X`)
  when `leverage = TRUE`, otherwise `NULL` entries.

- storage:

  Present only when `slim = TRUE`; records the model closures and
  leverage arrays that were removed.

## Details

![terradish logo](figures/terradish-sticker.png) A parameterized
conductance surface maps raster covariates to positive vertex
conductance. Vertex values define weighted edges and a graph Laplacian
\\L(\theta)\\. Its generalized inverse \\E(\theta) = L(\theta)^+\\ is
the shared graph-process kernel. Distance likelihoods use resistance
distances derived from \\E(\theta)\\, whereas covariance likelihoods use
the kernel directly. The selected probability model is called the
measurement model. These alternatives describe how related genetic
summaries are observed; they are not competing biological process
theories.

For example, a log-linear choice of conductance model is:

`vertex_conductance[i] = exp(covariates[i,] %*% theta)`

where `theta` are unknown parameters, `covariates` is a design matrix
with each row representing a vertex in the graph (cell in the raster).
The conductance of a given edge (an off-diagonal entry in the graph
Laplacian) is:

`edge_conductance[i,j] = vertex_conductance[i] + vertex_conductance[j]`

if argument `conductance = TRUE`, and

`edge_conductance[i,j] = 1/(1/vertex_conductance[i] + 1/vertex_conductance[j])`

otherwise. This alternative parameterization is not currently
implemented, so `conductance` must be `TRUE`. This differs by a factor
of two from CIRCUITSCAPE, where the edge conductance/resistance is the
average of the vertex conductance/resistance.

`terradish` estimates `theta` by maximizing the selected measurement
likelihood for the observed response. The optimization is done via
Newton's method (default; requires computation of Hessian during
optimization), via the BFGS algorithm (requires gradient only during
optimization) if `optimizer = "bfgs"`, or via a simple parameter-count
heuristic if `optimizer = "auto"`. Fitting one parameterized surface to
the genetic data, instead of regressing on separately computed
resistance distances, follows Peterman and Pope (2021). In all cases,
the final returned fit is evaluated at the optimized parameters using
the requested `curvature` for inference.

Conductance coefficients are conditional model associations. Their sign
and magnitude depend on covariate coding, other terms, graph
construction, response, measurement likelihood, and candidate set. They
are not direct estimates of causal effects, habitat suitability,
migration, dispersal, or organismal movement rates.

For an explanation of how categorical spatial covariates are handled,
see `details` of
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
and the examples below. The dummy coding of categorical covariates is
done by the function passed to `conductance_model` (e.g.
[`terradish::loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md)).

**Large rasters.** For large continuous rasters, two opt-in helpers can
reduce runtime. First, build the graph with
`conductance_surface(..., crop_buffer = )` when focal sites occupy only
part of the raster. This removes vertices outside the buffered sampling
extent before fitting. Cropping changes resistance: an unpublished
development check on one landscape found increases of 15–23% with a
two-cell buffer and 1–2% with a ten-cell buffer. These values are
case-specific. Refit with a larger buffer to assess sensitivity on your
landscape. Second, use `approximation = "coarse_raster"` with
`approximation_control = list(factor = c(4, 2), exact_refine = TRUE)` to
optimize on one or more aggregated rasters before refining on the
original graph. With `exact_refine = TRUE`, the returned coefficients
and likelihood are from the full-resolution graph; the coarse fits are
only starting values. To keep the full-resolution cleanup deliberately
short, pass `refine_control = NewtonRaphsonControl(maxit = 2, ...)`
inside `approximation_control`. Exact refinement is required.

For spline terms,
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
reports monotonicity over the covariate range at the focal sites,
including the number of derivative sign changes. A monotone curve
increases or decreases throughout that range; a nonmonotone curve
changes direction. This describes the fitted curve, not a confidence
statement about its shape outside or within that range.

**Gaussian scale-aware conductance.**
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md)
can also use `approximation = "coarse_raster"`. In that case the
Gaussian smoother is rebuilt on each aggregated raster, while the
full-resolution sigma bounds and unit conversion are preserved so a
warm-started `sigma` continues to mean the same map-unit distance at
every stage.

Currently, all of the built-in conductance models in `terradish` use the
default contrast coding (e.g. for a categorical covariate with `K`
factors, the estimated parameters are the `K-1` mean differences against
a reference category). The intercept is excluded if it is not
identifiable.

If the resistance contribution is zero, or the formula specifies uniform
conductance, conductance coefficients are not identified and their
influence, gradient, and Hessian are unavailable. A zero environmental
kernel weight does not by itself remove conductance inference.

**Convergence and starting values.** Inspect `fit$convergence` before
interpreting estimates. Code 0 requires a largest absolute projected
gradient below `ctol`, or an objective change below `ftol` together with
a projected gradient below `sqrt(ctol)`. Code 1 means the iteration
limit was reached; code 2 means a stall or failed line search. The
projected gradient ignores components pointing outside active parameter
bounds. The record includes the criterion, iterations, gradient,
boundary status, and whether a restart was attempted. A small objective
change alone is not convergence. Also inspect
`fit$fit$subproblem$convergence` for nuisance optimization. Numerical
convergence does not establish model adequacy.

If a user-supplied starting point reaches a no-structure boundary,
fitting retries the default starting point and keeps the better
likelihood. Report starting-value sensitivity when distinct solutions
remain plausible.

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

Peterman WE, Pope NS (2021). The use and misuse of regression models in
landscape genetic analyses. Molecular Ecology, 30, 37–47.
[doi:10.1111/mec.15716](https://doi.org/10.1111/mec.15716) .

Clarke RT, Rothery P, Raybould AF (2002). Confidence limits for
regression relationships between distance matrices: estimating gene flow
with distance. Journal of Agricultural, Biological, and Environmental
Statistics, 7, 361–372.
[doi:10.1198/108571102320](https://doi.org/10.1198/108571102320) .

Hanks EM, Hooten MB (2013). Circuit theory and model-based inference for
landscape connectivity. Journal of the American Statistical Association,
108, 22–33.
[doi:10.1080/01621459.2012.724647](https://doi.org/10.1080/01621459.2012.724647)
.

## See also

[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
to build `data`,
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
and
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
to read a fit,
[`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md)
and
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md)
for the model choices.

Other model fitting:
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md),
[`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md),
[`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md),
[`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)

## Examples

``` r
# \donttest{
library(terra)

data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)

# Name the layers first, then scale them. Scaling avoids numeric overflow
# and makes each coefficient a per-standard-deviation effect.
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

# A few melip sites share a raster cell; keep one site per cell.
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
coords <- coords[keep]
Fst <- melip.Fst[keep, keep]

# Build the graph (eight neighbors by default).
surface <- conductance_surface(covariates, coords)

# MLPE (the default) accounts for pairs that share a sampling site.
fit_mlpe <- terradish(Fst ~ altitude + forestcover, data = surface,
                      measurement_model = mlpe)
fit_mlpe$convergence$code   # 0 means converged; check before reading on
#> [1] 0
summary(fit_mlpe)
#> Conductance surface with 22443 vertices (34 focal) estimated by maximum likelihood
#> Call:   terradish(formula = Fst ~ altitude + forestcover, data = surface, 
#>     measurement_model = mlpe)
#> 
#> Loglikelihood: 1825.936 (6 degrees freedom)
#> AIC: -3639.873 
#> 
#> Number of function calls: 22 
#> Number of optimization steps: 6 
#> Convergence: converged | projected gradient: 4.82e-08 
#> Norm of gradient at MLE: 6.341795e-08 
#> 
#> Nuisance parameters (conditional on fitted conductance surface; measurement-model scale):
#>        Estimate  Std. Error  Lower 95%  Upper 95%
#> alpha  -0.04170   0.00853    -0.05842   -0.02498 
#> beta    0.33838   0.01138     0.31608    0.36068 
#> tau     6.64870   0.16712     6.32116    6.97625 
#> rho     0.65243   0.26065     0.14157    1.16330 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> altitude    -0.40667    0.06327  -6.428  1.3e-10 ***
#> forestcover  0.20608    0.10103   2.040   0.0414 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Correlation of Coefficients:
#>             altitude
#> forestcover 0.167607
exp(coef(fit_mlpe))         # conductance multiplier per one SD
#>    altitude forestcover 
#>   0.6658661   1.2288485 

# Least squares ignores shared-site dependence; useful for exploration.
fit_ls <- terradish(Fst ~ altitude + forestcover, data = surface,
                    measurement_model = leastsquares)

# Conductance surface with a 95% pointwise interval
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

# Test for an interaction using a likelihood-ratio test
fit_mlpe_interaction <- terradish(Fst ~ forestcover * altitude,
                                  data = surface, measurement_model = mlpe)
anova(fit_mlpe, fit_mlpe_interaction)
#> Likelihood ratio test: chi-square (1 df)
#> Null: ~ altitude + forestcover
#> Alt: ~ forestcover + altitude + forestcover:altitude
#>      logLik Df  ChiSq Df(ChiSq)  Pr(>Chi)    
#> Null 1825.9  6                               
#> Alt  1839.8  7 27.627         1 1.471e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Test against the uniform-conductance (isolation by distance) model
fit_mlpe_ibd <- terradish(Fst ~ 1, data = surface, measurement_model = mlpe)
anova(fit_mlpe, fit_mlpe_ibd)
#> Likelihood ratio test: chi-square (2 df)
#> Null: ~ 1
#> Alt: ~ altitude + forestcover
#>      logLik Df  ChiSq Df(ChiSq)  Pr(>Chi)    
#> Null 1783.4  4                               
#> Alt  1825.9  6 85.133         2 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Categorical covariates must be factor-valued layers; see the Details of
# ?conductance_surface.
forestcover <- terra::unwrap(melip.forestcover)
forestcover_class <- cut(terra::values(forestcover)[, 1],
                         breaks = c(0, 1/6, 1/3, 1))
forestcover_cat <- terra::setValues(forestcover,
                                    as.numeric(forestcover_class))
forestcover_cat <- terra::as.factor(forestcover_cat)
RAT <- levels(forestcover_cat)[[1]]
RAT$VALUE <- levels(forestcover_class) # explicitly define level names
levels(forestcover_cat) <- RAT

covariates_cat <- c(forestcover_cat, covariates[["altitude"]])
names(covariates_cat) <- c("forestcover", "altitude")
surface_cat <- conductance_surface(covariates_cat, coords)

# Coefficients for forestcover are contrasts against the first class.
fit_mlpe_cat <- terradish(Fst ~ forestcover + altitude, surface_cat,
                          measurement_model = mlpe)
summary(fit_mlpe_cat)
#> Conductance surface with 22443 vertices (34 focal) estimated by maximum likelihood
#> Call:   terradish(formula = Fst ~ forestcover + altitude, data = surface_cat, 
#>     measurement_model = mlpe)
#> 
#> Loglikelihood: 1824.911 (6 degrees freedom)
#> AIC: -3637.823 
#> 
#> Number of function calls: 22 
#> Number of optimization steps: 6 
#> Convergence: converged | projected gradient: 9.35e-08 
#> Norm of gradient at MLE: 1.010877e-07 
#> 
#> Nuisance parameters (conditional on fitted conductance surface; measurement-model scale):
#>        Estimate   Std. Error  Lower 95%  Upper 95%
#> alpha  -0.042416   0.008777   -0.059619  -0.025212
#> beta    0.474134   0.015981    0.442813   0.505456
#> tau     6.598406   0.171054    6.263146   6.933665
#> rho     0.729229   0.260313    0.219026   1.239433
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> forestcover   0.2105     0.1447   1.455    0.146    
#> altitude     -0.4231     0.0634  -6.674 2.49e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Correlation of Coefficients:
#>          forestcover
#> altitude   0.1410089
# }
```
