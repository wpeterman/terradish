# Print and read a power assessment

Compact display for the object returned by
[`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md):
the settings the sweep ran under, followed by the scenario-level summary
table.

## Usage

``` r
# S3 method for class 'terradish_covariance_power'
print(x, ...)
```

## Arguments

- x:

  An object of class `terradish_covariance_power`, from
  [`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md).

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the printed output.

## Details

**How to read the output.** The header records what the sweep did: how
many simulated data sets per design (`nsim`), which sample sizes, and
which site-selection strategies. Judge the numbers below it against that
`nsim`: with a small `nsim` the Monte Carlo error on a rate is
substantial, so treat a demonstration run as indicative only.

The table beneath is `x$summary`, one row per scenario (strategy by
sample size by fitted model). `fit_rate` is the proportion of
simulations that produced a usable fit; a value well below 1 can reflect
weak information, numerical difficulty, boundary fits, or a candidate
model that is poorly matched to the scenario. `conductance_power` is the
proportion exceeding the specified conductance-correlation threshold,
`mean_conductance_cor` summarizes correlation with the generating
conductance surface, and the model-selection columns report how often
each candidate had the lowest AIC or AICc among usable fits.

Per-coefficient results (power, bias, RMSE, and interval coverage) are
not printed, because there is one row per parameter per scenario. Read
them from `x$parameter_summary`. Power there is the proportion of
simulations in which the coefficient was both correctly signed and
significant, so it is the quantity to compare against a target such as
0.8.

## See also

[`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md)
for the sweep itself and what each returned component contains,
[`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)
for a single simulated data set, and
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
for fitting one. See
[`vignette("v10-simulation-design", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v10-simulation-design.md)
for a worked study-design walkthrough.

Other simulation and study design:
[`covariance_response_power()`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md),
[`simulate_covariance_response()`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)

## Examples

``` r
# \donttest{
data(melip)
covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")
# coarsened, and only a handful of simulations, so the example is quick
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
coords <- terra::unwrap(melip.coords)
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
surface <- conductance_surface(covariates, coords[keep])

power <- covariance_response_power(
  theta        = c(forestcover = 0.5, altitude = -0.4),
  formula      = ~ forestcover + altitude,
  data         = surface,
  sample_sizes = c(8, 12),
  strategies   = "spacefill",
  tau = 1, sigma = 0.2, nu = 100,
  nsim = 3,          # raise to 200+ for a real power estimate
  seed = 1,
  control = NewtonRaphsonControl(maxit = 8)
)

power                     # settings header plus the scenario summary
#> Covariance-response power assessment
#> Simulations per design: 3 
#> Sample sizes: 8, 12 
#> Strategies: spacefill 
#> 
#>   strategy sample_size model n_replicates  fit_rate conductance_power
#>  spacefill           8 truth            3 0.6666667         0.3333333
#>  spacefill          12 truth            3 1.0000000         1.0000000
#>  mean_conductance_cor mean_conductance_spearman selected_AIC_rate
#>             0.9331250                 0.9641851         0.6666667
#>             0.9942349                 0.9961422         1.0000000
#>  selected_AICc_rate mean_elapsed conductance_cor_threshold all_parameter_power
#>           0.6666667    0.5200000                       0.9           0.6666667
#>           1.0000000    0.5766667                       0.9           1.0000000
power$parameter_summary   # per-coefficient power, bias, RMSE, coverage
#>    strategy sample_size model   parameter true n_replicates  fit_rate     power
#> 1 spacefill           8 truth    altitude -0.4            3 0.6666667 0.6666667
#> 3 spacefill           8 truth forestcover  0.5            3 0.6666667 0.6666667
#> 2 spacefill          12 truth    altitude -0.4            3 1.0000000 1.0000000
#> 4 spacefill          12 truth forestcover  0.5            3 1.0000000 1.0000000
#>   correct_sign_rate significant_rate  coverage mean_estimate        bias
#> 1         0.6666667        0.6666667 0.6666667    -0.3627751  0.03722491
#> 3         0.6666667        0.6666667 0.6666667     0.8952228  0.39522278
#> 2         1.0000000        1.0000000 1.0000000    -0.4488384 -0.04883836
#> 4         1.0000000        1.0000000 1.0000000     0.5140838  0.01408379
#>         rmse
#> 1 0.05984194
#> 3 0.44084634
#> 2 0.05480261
#> 4 0.08743964
# }
```
