# Simulation-based power and parameter recovery for sampling designs

Runs a fast simulation study for covariance-response `terradish` fits.
Candidate sampling designs are created by subsetting the focal sites in
an existing
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
object, responses are generated with
[`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md),
and one or more candidate conductance models are refit to each simulated
covariance matrix. Gaussian smoothing widths are simulated and reported
in map units.

## Usage

``` r
covariance_response_power(
  theta,
  formula,
  data,
  sample_sizes,
  strategies = c("spacefill", "random"),
  conductance_model = loglinear_conductance,
  fit_models = NULL,
  tau = 1,
  sigma = 0,
  nu,
  nsim = 20L,
  n_designs = 1L,
  seed = NULL,
  theta_scale = c("external", "internal"),
  effect_parameters = NULL,
  alpha = 0.05,
  conductance_cor_threshold = 0.9,
  optimizer = c("auto", "newton", "bfgs"),
  control = NewtonRaphsonControl(maxit = 8, verbose = FALSE),
  solver = c("auto", "direct", "amg"),
  solver_control = NULL,
  approximation = c("none", "landmark", "coarse_raster"),
  approximation_control = NULL,
  nonnegative = TRUE,
  leverage = FALSE,
  cores = 1L,
  verbose = FALSE,
  nu_fit = nu
)
```

## Arguments

- theta:

  True conductance parameters. By default these are interpreted on the
  same external scale returned by
  [`coef()`](https://rdrr.io/r/stats/coef.html) for fitted `terradish`
  models. For most conductance models the internal and external scales
  are identical. For
  [`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
  the external scale reports `sigma.*` parameters in map units.

- formula:

  Formula for the true conductance model. The left hand side, if
  supplied, is ignored.

- data:

  A
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
  object containing all candidate focal sites.

- sample_sizes:

  Integer vector of focal sample sizes to assess.

- strategies:

  Sampling strategies. Supported values are `"spacefill"`, `"random"`,
  and `"sequential"`.

- conductance_model:

  Conductance-model factory used to simulate the true response.

- fit_models:

  Optional named list of candidate models to refit. Each element may
  contain `formula`, `conductance_model`, `theta`, `optimizer`,
  `control`, `solver`, `solver_control`, `approximation`,
  `approximation_control`, `nonnegative`, and `leverage`. When omitted,
  the true model is refit.

- tau:

  Nonnegative scaling applied to the conductance-implied covariance
  matrix during simulation.

- sigma:

  Nonnegative nugget variance added to the covariance diagonal during
  simulation, on its natural (variance) scale. The fitted Wishart
  `sigma` is the logarithm of this quantity.

- nu:

  Effective Wishart information used to simulate the response. Larger
  values reduce simulated sampling noise. The fitting value defaults to
  it through `nu_fit = nu`; neither is an estimate of information from
  real data. Report a sensitivity analysis rather than substituting the
  marker count. Simulation requires `nu` to be at least the largest
  assessed covariance-matrix dimension.

- nsim:

  Number of covariance-response simulations per sampling design.

- n_designs:

  Number of independent site designs per sample size for the `"random"`
  strategy. Deterministic strategies are evaluated once.

- seed:

  Optional random seed for reproducible site and response simulation.

- theta_scale:

  Is `theta` supplied on the `"external"` or `"internal"`
  conductance-model scale?

- effect_parameters:

  Optional character vector of parameters used in the all-effects
  recovery summary. Defaults to the nonzero true parameters.

- alpha:

  Significance level for Wald-style parameter recovery summaries.

- conductance_cor_threshold:

  Threshold used for the conductance-surface recovery summary.

- optimizer, control, solver, solver_control, approximation,
  approximation_control, nonnegative, leverage, cores:

  Defaults passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
  for candidate models unless a model-specific value is supplied in
  `fit_models`.

- verbose:

  Logical; print progress messages?

- nu_fit:

  Effective Wishart degrees of freedom used to fit each simulation.
  Defaults to the generating `nu`. Increasing this value narrows fitted
  intervals without adding independent information. Compare coverage as
  well as power when exploring this sensitivity.

## Value

An object of class `"terradish_covariance_power"` containing:

- `results`:

  One row per fitted model and simulated response.

- `summary`:

  Scenario-level fit, model-selection, and conductance-recovery
  summaries.

- `parameter_results`:

  One row per assessed parameter and simulation replicate.

- `parameter_summary`:

  Parameter-level power, bias, RMSE, and coverage summaries.

- `settings`:

  Simulation settings and true parameters, including `nu` for simulation
  and `nu_fit` for inference.

Power and coverage are internal to the assumed response model. Matching
simulation and fitting values cannot diagnose incorrect information in a
real data set. Nonconverged fits count as failed fits and
non-detections; both outer and final nuisance optimization must record
code 0. Missing convergence diagnostics also count as failed fits.
Inspect `fit_rate` before interpreting power or coverage.

## Details

This helper avoids forward-time genetic simulation. It assumes a
proposed conductance surface and draws covariance responses from the
same Wishart likelihood fitted by
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).
That makes it useful for approximate power and sample-size screening
before investing in slower individual-based simulations.

Sampling designs are efficient because the raster graph is not rebuilt
for each sample size. The helper subsets `data$demes` and the
corresponding right-hand-side columns, so all candidate designs share
the same graph, covariate matrix, and conductance-model basis expansion.

Fits that land on the no-signal boundary are retained as non-detections
in parameter-level power summaries. Their coefficient estimates and
standard errors are reported as `NA`, matching explicit fitting
failures.

For spline conductance models, individual basis coefficients are not
usually the most interpretable recovery target. Inspect
`parameter_summary` for basis-level behavior, but prefer
`summary$conductance_power`, `summary$mean_conductance_cor`, and
model-selection rates when judging whether the non-linear conductance
surface was recovered.

## See also

[`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md),
[`print.terradish_covariance_power`](https://wpeterman.github.io/terradish/reference/print.terradish_covariance_power.md)

Other simulation and study design:
[`print.terradish_covariance_power()`](https://wpeterman.github.io/terradish/reference/print.terradish_covariance_power.md),
[`simulate_covariance_response()`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)

## Examples

``` r
# \donttest{
# Coarsened: a power sweep refits many times, so cell count is the
# dominant cost. Use the full resolution for a real design study.
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

spline_model <- smooth_loglinear_conductance(
  ~ forestcover + s(altitude, df = 4), surface$x)
theta <- attr(spline_model, "default")
theta["forestcover"] <- -0.15
theta[grepl("^s\\(altitude\\)", names(theta))] <-
  c(-0.65, 1.1, -0.9, 0.45)

gaussian_model <- gaussian_smoothed_loglinear_conductance(
  surface, scale_vars = "altitude")

power <- covariance_response_power(
  theta = theta,
  formula = ~ forestcover + s(altitude, df = 4),
  data = surface,
  sample_sizes = c(8, 12),
  strategies = "spacefill",
  conductance_model = smooth_loglinear_conductance,
  fit_models = list(
    spline = list(formula = ~ forestcover + s(altitude, df = 4),
                  conductance_model = smooth_loglinear_conductance),
    gaussian = list(formula = ~ altitude + forestcover,
                    conductance_model = gaussian_model)
  ),
  tau = 0.8,
  sigma = 0.12,
  nu = 150,
  nsim = 3,          # raise to 200+ for a real power estimate
  n_designs = 1,
  seed = 1,
  control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
)
power$summary
#>    strategy sample_size    model n_replicates  fit_rate conductance_power
#> 1 spacefill           8 gaussian            3 0.6666667                 0
#> 3 spacefill           8   spline            3 0.0000000                 0
#> 2 spacefill          12 gaussian            3 0.3333333                 0
#> 4 spacefill          12   spline            3 0.0000000                 0
#>   mean_conductance_cor mean_conductance_spearman selected_AIC_rate
#> 1            0.4212063                 0.4399261         0.6666667
#> 3                   NA                        NA         0.0000000
#> 2            0.8064656                 0.7996254         0.3333333
#> 4                   NA                        NA         0.0000000
#>   selected_AICc_rate mean_elapsed conductance_cor_threshold all_parameter_power
#> 1          0.6666667    0.6566667                       0.9                   0
#> 3          0.0000000    0.3966667                       0.9                   0
#> 2          0.3333333    0.7066667                       0.9                   0
#> 4          0.0000000    0.5000000                       0.9                   0
power$parameter_summary
#>     strategy sample_size    model     parameter  true n_replicates  fit_rate
#> 1  spacefill           8 gaussian   forestcover -0.15            3 0.6666667
#> 3  spacefill           8   spline   forestcover -0.15            3 0.0000000
#> 5  spacefill           8   spline s(altitude).1 -0.65            3 0.0000000
#> 7  spacefill           8   spline s(altitude).2  1.10            3 0.0000000
#> 9  spacefill           8   spline s(altitude).3 -0.90            3 0.0000000
#> 11 spacefill           8   spline s(altitude).4  0.45            3 0.0000000
#> 2  spacefill          12 gaussian   forestcover -0.15            3 0.3333333
#> 4  spacefill          12   spline   forestcover -0.15            3 0.0000000
#> 6  spacefill          12   spline s(altitude).1 -0.65            3 0.0000000
#> 8  spacefill          12   spline s(altitude).2  1.10            3 0.0000000
#> 10 spacefill          12   spline s(altitude).3 -0.90            3 0.0000000
#> 12 spacefill          12   spline s(altitude).4  0.45            3 0.0000000
#>    power correct_sign_rate significant_rate  coverage mean_estimate       bias
#> 1      0         0.3333333                0 0.6666667    0.23425530 0.38425530
#> 3      0         0.0000000                0 0.0000000            NA         NA
#> 5      0         0.0000000                0 0.0000000            NA         NA
#> 7      0         0.0000000                0 0.0000000            NA         NA
#> 9      0         0.0000000                0 0.0000000            NA         NA
#> 11     0         0.0000000                0 0.0000000            NA         NA
#> 2      0         0.3333333                0 0.3333333   -0.06667512 0.08332488
#> 4      0         0.0000000                0 0.0000000            NA         NA
#> 6      0         0.0000000                0 0.0000000            NA         NA
#> 8      0         0.0000000                0 0.0000000            NA         NA
#> 10     0         0.0000000                0 0.0000000            NA         NA
#> 12     0         0.0000000                0 0.0000000            NA         NA
#>          rmse
#> 1  0.45536774
#> 3          NA
#> 5          NA
#> 7          NA
#> 9          NA
#> 11         NA
#> 2  0.08332488
#> 4          NA
#> 6          NA
#> 8          NA
#> 10         NA
#> 12         NA
# }
```
