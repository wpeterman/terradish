# terradish <img src="man/figures/terradish-sticker.png" align="right" height="139" alt="terradish logo" />

| CI | CRAN | License | Archive |
|:---|:---|:---|:---|
| [![R-CMD-check](https://github.com/wpeterman/terradish/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/wpeterman/terradish/actions/workflows/R-CMD-check.yaml) | [![CRAN status](https://www.r-pkg.org/badges/version/terradish)](https://CRAN.R-project.org/package=terradish) | [![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21225711.svg)](https://doi.org/10.5281/zenodo.21225711) |
| [![codecov](https://codecov.io/gh/wpeterman/terradish/graph/badge.svg)](https://app.codecov.io/gh/wpeterman/terradish) [![pkgdown](https://github.com/wpeterman/terradish/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/wpeterman/terradish/actions/workflows/pkgdown.yaml) | [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/terradish)](https://CRAN.R-project.org/package=terradish) | | |

## Overview

terradish estimates landscape conductance surfaces from genetic data by maximum likelihood on a raster graph. It fits log-linear conductance models and can optionally estimate each covariate's Gaussian smoothing scale and spline shape.

Models fit pairwise genetic distances through the MLPE likelihood, or allele-frequency covariance and squared distances derived from it through a Wishart likelihood on site contrasts. Both families accept the same pairwise environmental covariates alongside resistance. Exact sparse Cholesky and algebraic multigrid solvers handle large rasters. Fixed-domain spatial cross-validation compares conductance formulas.

Estimated coefficients describe relative conductance and are conditional on the measurement-model terms. They do not directly estimate movement, migration, or causal landscape effects. For Wishart likelihoods, the user-supplied information parameter `nu` controls inferential precision. Using the SNP count overstated information by roughly an order of magnitude in the package's forward-time validation experiments.

## Citation

To cite terradish or acknowledge its use, run `citation("terradish")` for the current reference and a BibTeX entry. The archived releases share the DOI [10.5281/zenodo.21225711](https://doi.org/10.5281/zenodo.21225711).

Fitting one parameterized conductance surface to genetic data, in place of regressing on separate resistance distances, follows Peterman and Pope (2021), [doi:10.1111/mec.15716](https://doi.org/10.1111/mec.15716).

## Installation

There are two ways to install terradish:

```r
# Get the released version from CRAN
install.packages("terradish")

# Alternatively, install the development version from GitHub
# install.packages("remotes")
remotes::install_github("wpeterman/terradish")
```

terradish requires the companion package landgraph, version 0.0.3 or later. Both routes install it from CRAN automatically.

## Using terradish

The full documentation, with every guide and help page, is at <https://wpeterman.github.io/terradish/>.

A terradish analysis has three steps: build a raster graph from covariates and sampling locations with `conductance_surface()`, fit a conductance model to a genetic response with `terradish()`, and then check, compare, and map the fit.

### Quick start

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

# A few sites fall in the same raster cell. Keep one site per cell so that
# each graph vertex holds one sampled population.
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(sites)))
sites <- sites[keep]
Fst <- melip.Fst[keep, keep]
graph <- conductance_surface(covariates, sites)  # eight neighbors by default

# MLPE accounts for shared sampling sites among pairwise genetic distances.
fit <- terradish(Fst ~ forestcover + altitude, graph,
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

joint <- terradish(Fst ~ forestcover + altitude, graph,
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

Cropping changes resistance. One unpublished development check found increases of 15 to 23% with a two-cell buffer and 1 to 2% with ten cells. These are case-specific results: compare larger-buffer refits. Landmark and coarse-raster starts require exact refinement; `terradish_grid()` supports exact evaluation only. `terradish_assess_settings()` is a speed probe, not a convergence or statistical check.

## Experimental features

Directed conductance, hierarchical fields, site-specific drift terms, pair subsets, Kron reduction, block-CG, and legacy CV remain on the `experimental` branch. Their unresolved limitations are documented in that branch's `EXPERIMENTAL.md`.

```r
remotes::install_github("wpeterman/terradish@experimental")
```

## Guides

The guides are numbered in reading order. Start with `vignette("v01-getting-started", package = "terradish")`, or with `vignette("v03-mlpe-pairwise-distance", package = "terradish")` if your data are a pairwise genetic distance matrix such as FST. Open the collection with `browseVignettes("terradish")`.

| Guide | Question it answers |
|:---|:---|
| `v01-getting-started` | How do I run one complete analysis, from genotypes to a conductance map? |
| `v02-wishart-covariance` | Which genetic covariance response should I model, and is it admissible? |
| `v03-mlpe-pairwise-distance` | How do I fit and read a model of pairwise genetic distances? |
| `v04-model-comparison` | Which candidate model do information criteria and spatial cross-validation support? |
| `v05-categorical-conductance` | How do I fit and interpret land-cover classes? |
| `v06-spline-conductance` | Does a curved association predict better than a straight one? |
| `v07-gaussian-scale-optimization` | Over what neighborhood does a covariate act? |
| `v08-ibe-ibr-workflow` | Do environmental differences matter beyond resistance? |
| `v09-individual-salamanders` | How do individual genotypes reveal a nonlinear resistance response? |
| `v10-simulation-design` | How many sites would a study need? |
| `v11-large-landscapes` | Which solver, storage, and domain settings suit a large raster? |
| `v12-data-generation` | How were the bundled salamander data simulated? |

Changes in each release are listed in `NEWS.md` and in the website changelog.

## Attribution

Nate Pope developed the sparse graph optimization, reverse-mode derivatives, MLPE, and generalized Wishart foundations in [radish](https://github.com/nspope/radish). terradish extends that framework with terra-native workflows and the supported features above. Selected deprecated `radish*` wrappers remain; full compatibility with every historical entry point is not promised.

Report bugs and request features through the [issue tracker](https://github.com/wpeterman/terradish/issues).
