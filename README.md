# terradish <img src="man/figures/terradish-sticker.png" align="right" height="200"/>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21225712.svg)](https://doi.org/10.5281/zenodo.21225712)

terradish estimates landscape conductance surfaces from genetic data by maximum likelihood on a raster graph. It fits log-linear conductance models and can optionally estimate each covariate's Gaussian smoothing scale and spline shape. Models fit pairwise genetic distances through the MLPE likelihood, or allele-frequency covariance and squared distances derived from it through a Wishart likelihood on site contrasts. Both families accept the same pairwise environmental covariates alongside resistance. Exact sparse Cholesky and algebraic multigrid solvers handle large rasters. Fixed-domain spatial cross-validation compares conductance formulas.

Estimated coefficients describe relative conductance and are conditional on the measurement-model terms. They do not directly estimate movement, migration, or causal landscape effects. For Wishart likelihoods, the user-supplied information parameter `nu` controls inferential precision. Using the SNP count overstated information by roughly an order of magnitude in the package's forward-time validation experiments.

## Installation

The supported core requires landgraph version 0.0.3 or later. Install a compatible companion version before terradish:

```r
remotes::install_github("wpeterman/landgraph")
remotes::install_github("wpeterman/terradish")
```

For local development with both repositories checked out side by side, install landgraph from its checkout first. The terradish CRAN release depends on that companion version being available on CRAN.

## Quick start

```r
library(terradish)
library(terra)
data(melip)

# Unwrap the portable example data and standardize each raster layer.
altitude <- terra::unwrap(melip.altitude)
forestcover <- terra::unwrap(melip.forestcover)
sites <- terra::unwrap(melip.coords)
covariates <- c(altitude, forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)
graph <- conductance_surface(covariates, sites)  # eight neighbors by default

# MLPE accounts for shared sampling sites among pairwise genetic distances.
fit <- terradish(melip.Fst ~ forestcover + altitude, graph,
                 measurement_model = mlpe)
fit$convergence
summary(fit)
coef(fit)
confint(fit)
vcov(fit)
plot(fit, type = "surface", data = graph)
plot(fit, type = "fit")
```

Check convergence before interpreting the surface. Code 0 meets the stopping rule; code 1 reaches the iteration limit; code 2 reports a stall or failed line search. Also inspect `fit$fit$subproblem` for nuisance optimization. Numerical convergence does not establish model adequacy.

## Interpreting estimates

For a standardized linear covariate, `exp(coef(fit))` is the conductance multiplier for a one-standard-deviation increase, holding other terms fixed. A common conductance scale is absorbed by the measurement model's resistance coefficient, so absolute conductance is unidentified.

Gaussian `sigma` describes raster smoothing in map units. The smoothed layer is re-standardized at every candidate sigma, so its coefficient is per standard deviation of that layer. The half-cell lower bound still smooths, and the upper bound is limited by the retained kernel window. Dispersal averaging can favor positive scales without a separate ecological scale of effect. Use `gaussian_scale_summary()`, `summary(fit)$sigma_table`, and `gaussian_scale_profile()` to examine units, bounds, and uncertainty.

Spline models are unpenalized regression splines with fixed degrees of freedom. Their centered bases and knots are retained for prediction. `summary(fit)$spline_monotonicity` describes shape over the focal-site covariate range. Tails outside that range are poorly identified.

| Response | Measurement model | Main requirement |
|---|---|---|
| Pairwise genetic distances, including ratio-estimator FST | `mlpe` or `leastsquares` | A defensible mean structure; MLPE models shared-site correlation |
| Allele-frequency covariance | `wishart_covariance` | Coherent covariance and explicit effective information |
| Squared distances derived from covariance | `generalized_wishart` | Admissible geometry and explicit effective information |

The two Wishart forms use the same site-contrast likelihood for matching covariance and distance representations. Passing `check_distance_response()` does not give an FST ratio estimator a Wishart sampling model. Grouped covariance should use the coherent `gower` diagonal from landgraph; `within` is on a different scale. Rare-variant weighting and unequal group sizes require sensitivity checks.

### What `nu` changes

Wishart point estimates are invariant to rescaling `nu` at the optimum. Standard errors scale as `1 / sqrt(nu)`, while likelihood-ratio statistics and the likelihood contribution to information criteria scale with `nu`. Linkage and shared population history can make effective information much smaller than the marker count.

Use `terradish_rescale_nu(fit, nu = ...)` to report conclusions over plausible values. Dividing `nu` by 20 increases standard errors by `sqrt(20)` without refitting estimates. AIC at nominal marker information can over-select environmental terms and flexible curves. Use spatial CV for predictive support and rescaling for inferential sensitivity.

## Combine IBE and IBR

```r
# Raster covariates come first; coordinates come second.
environment <- pairwise_endpoint_covariates(covariates, sites)
lonlat <- terra::crds(terra::project(sites, "EPSG:4326"))
local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
                     mean(lonlat[, 2]), mean(lonlat[, 1]))
projected_sites <- terra::project(sites, local_crs)
pairwise <- pairwise_covariates(environment,
  geographic_100km = dist(terra::crds(projected_sites)) / 100000)

joint <- terradish(melip.Fst ~ forestcover + altitude, graph,
                   measurement_model = mlpe_covariates(pairwise))
terradish_ibe_ratio(joint)

# The same object can enter a model for a coherent covariance response.
wishart_measurement <- wishart_covariates(pairwise, model = "wishart_covariance")
```

Geographic differences in this example are in units of 100 km. The IBE:IBR ratio expresses an environmental difference in resistance-distance units with a joint uncertainty interval. Raw coefficients are not directly comparable across families. Ratios require matching transforms and scaling and are undefined when the resistance coefficient is zero.

Environmental terms can reflect assortative mating, dispersal filtering, local adaptation, or environmental patterns in site variance. Report conductance fits with and without pairwise terms because signs can change. A uniform `~ 1` surface with environmental terms represents **IBD + IBE**, not IBE alone.

## Model comparison

`terradish_folds()` constructs folds; `terradish_cv_folds()` keeps the landscape graph fixed while withholding sites. It accepts per-formula conductance factories and repeated folds. Reprofiled nuisance parameters compare conductance formulas under one measurement model. Set `nuisance = "fixed"` to retain every measurement parameter from training and score predictive densities when comparing measurement extensions.

Read success counts first. Failed folds make full totals unavailable, and rankings use only folds successful for every candidate. Paired differences and their standard errors are descriptive because folds share data. Checkpoints reject incompatible input or model changes.

MLPE and Wishart scores do not share a common predictive scale. Within one Wishart comparison, changing `nu` rescales scores without changing ranks when the same folds succeed. Likelihood-ratio tests require nested models, matching responses and graphs, and fixed `nu`. Environmental Wishart weights need a boundary reference, not an ordinary chi-square test.

## Prediction and large landscapes

New-landscape prediction reuses fitted log-linear terms, spline knots and centers, or Gaussian post-smoothing scaling. Apply stored input scaling with `scale_covariates(new_rasters, reference = original_scaled_rasters)` before building a new graph. Slim fits retain covariance and interval methods but cannot predict after their model closures are removed.

`solver = "auto"` selects AMG above the large-graph threshold regardless of the number of right-hand sides. AMG iteration counts can increase with conductance contrast. Experimental parallel derivatives use PSOCK workers for Hessian and partial solves; AMG and cached CHOLMOD ignore `cores`.

Cropping changes resistance. One audit found increases of 15–23% with a two-cell buffer and 1–2% with ten cells. These are case-specific results: compare larger-buffer refits. Landmark and coarse-raster starts require exact refinement; `terradish_grid()` supports exact evaluation only. `terradish_assess_settings()` is a speed probe, not a convergence or statistical check.

## Experimental features

Directed conductance, hierarchical fields, site-specific drift terms, pair subsets, Kron reduction, block-CG, and legacy CV remain on the `experimental` branch. Their unresolved limitations are documented in that branch's `EXPERIMENTAL.md`.

```r
remotes::install_github("wpeterman/terradish@experimental")
```

## Changes in 1.0.0

- The supported core concentrates on log-linear, Gaussian, and spline conductance with MLPE and contrast Wishart likelihoods.
- Graphs use eight directions by default. Gaussian alignment, kernel bounds, spline prediction, and environmental kernels are corrected.
- Convergence records, covariance and interval methods, Gaussian profiles, IBE ratios, and `nu` sensitivity are available.
- Fixed-domain CV supports per-formula factories, fixed nuisance prediction, repeated folds, failure accounting, and compatible checkpoints.
- Stored scaling, new-landscape prediction, combined pairwise covariates, and `nu_fit` sensitivity are supported.
- Research prototypes and unsupported solvers have moved off the supported core. See `NEWS.md` for migration details.

## Vignettes and attribution

Start with `vignette("getting-started", package = "terradish")`. Other guides cover comparison, IBE and IBR, Wishart covariance, splines, Gaussian scales, large landscapes, and simulation design. Open the collection with `browseVignettes("terradish")`.

Nate Pope developed the sparse graph optimization, reverse-mode derivatives, MLPE, and generalized Wishart foundations in [radish](https://github.com/nspope/radish). terradish extends that framework with terra-native workflows and the supported features above. Selected deprecated `radish*` wrappers remain; full compatibility with every historical entry point is not promised.

Contact Bill Peterman (Peterman.73@osu.edu) or submit an issue on GitHub for bug reports and feature requests.
