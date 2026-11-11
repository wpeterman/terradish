# Spatial cross-validation of conductance formulas on a fixed graph

Fits one or more conductance formulas to each training fold while
retaining the complete landscape graph, then evaluates the trained
coefficients on the held-out focal sites. Measurement parameters can be
reprofiled on each test fold or held at their training estimates for
predictive scoring.

## Usage

``` r
terradish_cv_folds(
  data,
  formulas,
  folds,
  model = "mlpe",
  nu = NULL,
  conductance_model = loglinear_conductance,
  baseline = TRUE,
  keep_fits = c("none", "slim", "full"),
  checkpoint = NULL,
  resume = TRUE,
  cores = 1L,
  measurement_control = NULL,
  nuisance = c("reprofile", "fixed"),
  ...
)

# S3 method for class 'terradish_cv_folds'
print(x, ...)
```

## Arguments

- data:

  A prebuilt `terradish_graph`, usually returned by
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- formulas:

  A formula or named list of formulas. All formulas must use the same
  response matrix and likelihood family.

- folds:

  Fold membership for each focal site. Use
  [`terradish_folds`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
  or supply study-specific labels. A list of fold vectors requests
  repeated cross-validation.

- model:

  Measurement model function or one of `"mlpe"`
  ([`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md)),
  `"ls"`
  ([`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md)),
  or `"wishart"`. The string `"wishart"` selects
  [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
  which expects a squared-distance response. For a covariance-matrix
  response pass the function itself, `model = wishart_covariance`. With
  fixed nuisance parameters, a named list matching the formula names can
  compare models within one family.

- nu:

  Effective Wishart degrees of freedom when required by the selected
  measurement model.

- conductance_model:

  One conductance-model factory or a named list matching the formula
  names, allowing log-linear, Gaussian, and spline models in the same
  comparison.

- baseline:

  Should each test fold also be scored at zero conductance coefficients,
  the uniform-conductance baseline on the same graph?

- keep_fits:

  Whether to retain no training fits, slim fits, or full fits.

- checkpoint:

  Optional file path. Completed fold results are written atomically
  after every successful or failed fit so a long run can be inspected or
  resumed.

- resume:

  Logical. If `TRUE` and `checkpoint` exists, continue a compatible run
  without reevaluating completed formula-fold combinations.

- cores:

  Number of workers used inside each model fit. Fold evaluation is
  serial to avoid nested parallel clusters.

- measurement_control:

  Optional
  [`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  object for profiling nuisance parameters.

- nuisance:

  `"reprofile"` (default) estimates measurement parameters separately on
  each held-out block. `"fixed"` holds every measurement parameter,
  including environmental coefficients, at its training estimate.

- ...:

  Additional arguments passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).
  An explicit `optimizer` is used for each training fit. Otherwise the
  training call tries Newton and retries with BFGS if Newton raises an
  error. For [`print()`](https://rdrr.io/r/base/print.html), further
  arguments passed to the data-frame print method.

- x:

  An object of class `"terradish_cv_folds"`, for
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `"terradish_cv_folds"`. Its `folds` element records
the site assignments, `results` contains one row per formula and fold,
and `summary` reports summed and mean held-out scores by formula,
common-fold totals, paired differences from the best model, the standard
error of their mean, and a one-standard-error flag. Repeated runs also
report the SD of repeat totals over common folds. Optional retained fits
are stored in `fits`. [`print()`](https://rdrr.io/r/base/print.html)
returns its input invisibly.

Columns of `results`: `model`, `fold`, `repetition`, `n_train` and
`n_test` (sites), `loglik` (held-out log-likelihood of the trained
model; larger is better), `baseline_loglik` (the same for uniform
conductance), `loglik_gain` (their difference), `nuisance_iterations`
and `nuisance_convergence`, `elapsed_seconds`, `error` and
`baseline_error` (messages for failed fits or scores), `common_fold`
(`TRUE` when every model succeeded on that fold), and
`paired_difference` (the model's `loglik` minus the best model's on the
same fold).

Columns of `summary`, which
[`print()`](https://rdrr.io/r/base/print.html) displays with the best
model first: `folds`, `failed`, and `n_folds_ok` (fold counts; read
these first), `total_loglik` and `mean_loglik` over all of a model's
folds (`NA` total when any fold failed), `common_folds_total` and
`n_common_folds` (the total over folds on which every model succeeded;
models are ranked by this column), `between_repeat_sd`,
`total_loglik_gain` and `mean_loglik_gain` (improvement over uniform
conductance), `mean_difference` and `se_difference` (mean paired
difference from the best model and its standard error across folds), and
`within_one_se` (`TRUE` when that mean difference is no larger than its
standard error, so the model is not clearly worse than the best one).

## Details

The graph domain and raster covariates remain fixed. Only the focal-site
entries of the graph (its `demes` index and the matching solver
right-hand sides), the response matrix, and any measurement-model site
attributes are subset. This evaluates transfer to held-out sites within
the observed landscape rather than rebuilding a new landscape for every
fold.

Reprofiled scores compare conductance formulas under one measurement
model; they must not be compared across calls with different measurement
models or effective degrees of freedom. Fixed-parameter scores are
predictive marginal densities and can compare measurement models within
one likelihood family. Neither mode scores pairs spanning training and
test sites. Spatial folds score only within-cluster pairs. Wishart score
gains scale with `nu`, although ideal fitted rankings are invariant to
its scale.

Checkpoints include a signature of the response, graph, focal cells,
covariates, model definitions, formulas, `nu`, and fold assignments. For
factor covariates the signature retains class placement, levels, and
contrasts. Uniform `~1` formulas are scored with unit conductance. A
mismatched signature stops resumption. Nonconverged fits and failed
scores count as failed folds. Training fits must record code 0 for both
outer and final nuisance optimization; missing diagnostics also fail the
fold. Checkpoints predating this convergence requirement are
incompatible and must be recomputed (use a new checkpoint or
`resume = FALSE`). With fixed nuisance parameters, no-structure boundary
fits remain scoreable because predictions do not depend on conductance.
A failed uniform-baseline score is recorded separately: model scores
remain usable, but baseline gains are unavailable for that fold. Models
with failed folds receive NA totals and an explicit success count.
Rankings use only folds successful for every model. Paired differences
and their standard errors are descriptive; repeated folds share
observations and are not independent replicates.

`loglik_gain` is the trained model's held-out log-likelihood minus the
uniform-conductance baseline log-likelihood. Larger values indicate
improved prediction over that baseline. Raw log-likelihood and gain are
comparable only among models using the same response representation,
measurement likelihood, folds, graph, and effective degrees of freedom.
They are not a common scoring rule for comparing MLPE,
generalized-Wishart distance, and covariance-Wishart models.

## Worked examples

See
[`vignette("v04-model-comparison", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v04-model-comparison.md)
for pond comparisons and
[`vignette("v09-individual-salamanders", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v09-individual-salamanders.md)
for individual sampling, matched markers, and spatial prediction plots.

## See also

[`terradish_folds`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

## Examples

``` r
# \donttest{
r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8,
                 ymin = 0, ymax = 8, crs = "EPSG:32617")
terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
names(r) <- "x"
sites <- terra::xyFromCell(r, c(1, 5, 8, 12, 20, 27, 33, 37, 45, 53, 57, 64))
surface <- conductance_surface(r, sites)
response <- simulate_covariance_response(c(x = 0.8), ~x, surface,
  nu = 1000, seed = 93)$covariance
folds <- terradish_folds(sites, k = 3, method = "random", seed = 1)
cv <- terradish_cv_folds(
  surface,
  list(uniform = response ~ 1, gradient = response ~ x),
  folds = folds,
  model = wishart_covariance, nu = 1000, nuisance = "fixed",
  control = NewtonRaphsonControl(maxit = 100))
# Check fold success counts (failed, n_folds_ok) before comparing the
# common-fold predictive scores. The best model is listed first.
cv
#> Fixed-domain terradish cross-validation
#> Measurement model: wishart_covariance 
#> Repeats: 1 | nuisance: fixed 
#> Ranks use common successful folds; differences are paired against the best model.
#> 
#>     model folds failed n_folds_ok total_loglik common_folds_total
#>  gradient     3      0          3     6462.632           6462.632
#>   uniform     3      0          3     6324.266           6324.266
#>  n_common_folds between_repeat_sd mean_loglik total_loglik_gain
#>               3                NA    2154.211          138.3657
#>               3                NA    2108.089            0.0000
#>  mean_loglik_gain se_difference mean_difference within_one_se
#>          46.12189       0.00000         0.00000          TRUE
#>           0.00000      17.04508       -46.12189         FALSE
cv$summary[, c("model", "n_folds_ok", "common_folds_total",
               "mean_difference", "se_difference", "within_one_se")]
#>      model n_folds_ok common_folds_total mean_difference se_difference
#> 1 gradient          3           6462.632         0.00000       0.00000
#> 2  uniform          3           6324.266       -46.12189      17.04508
#>   within_one_se
#> 1          TRUE
#> 2         FALSE
# Random folds illustrate the call; choose ecological spatial blocks in practice.
# }
```
