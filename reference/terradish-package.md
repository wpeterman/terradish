# terradish: Fast Gradient-Based Optimization of Resistance Surfaces

Estimates landscape conductance surfaces from genetic data by maximum
likelihood on a raster graph. It fits log-linear conductance models and
can optionally estimate each covariate's Gaussian smoothing scale and
spline shape. Models fit pairwise genetic distances through the MLPE
likelihood, or allele-frequency covariance and squared distances derived
from it through a Wishart likelihood on site contrasts. Both families
accept the same pairwise environmental covariates alongside resistance.
Exact sparse Cholesky and algebraic multigrid solvers handle large
rasters; fixed-domain spatial cross-validation compares conductance
formulas.

## Details

![terradish logo](figures/terradish-sticker.png)

**Where to start.**
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
builds the graph,
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
fits the model, and the methods documented in
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
([`summary()`](https://rspatial.github.io/terra/reference/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`plot()`](https://rspatial.github.io/terra/reference/plot.html),
[`anova()`](https://rdrr.io/r/stats/anova.html)) read the result.
Conductance parameters are on the log scale, so `exp(coef(fit))` is the
multiplicative change in conductance per unit of a covariate. These
coefficients describe conditional associations within the fitted graph
and candidate set. They do not, by themselves, identify causal
mechanisms, absolute migration or dispersal rates, or habitat
suitability. For Wishart models, the user-supplied effective degrees of
freedom, `nu`, sets inferential precision. Using the SNP count
substantially overstated information in the package's forward-time
validation experiments. Report sensitivity across plausible values with
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md).

**Main components.**

- Conductance models:

  [`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md)
  (the default),
  [`smooth_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md)
  for spline responses, and
  [`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md)
  for estimating a covariate's Gaussian raster-smoothing scale inside
  the model.

- Measurement models:

  [`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md)
  and [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md)
  for a genetic distance matrix;
  [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
  and
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
  for a Wishart likelihood, which need the effective degrees of freedom
  `nu`;
  [`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md)
  to verify a generalized-Wishart squared-distance response;
  [`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
  and
  [`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)
  to add isolation-by-environment terms.

- Pairwise inputs and inference:

  [`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
  combines geographic distances and environmental differences;
  [`terradish_ibe_ratio`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md)
  expresses their fitted effects in resistance-distance units.
  [`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
  provides covariance matrices and confidence intervals, and
  [`gaussian_scale_profile`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md)
  profiles smoothing-scale uncertainty.

- Comparison and diagnostics:

  [`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md),
  [`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
  [`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
  and
  [`terradish_assess_settings`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md),
  a speed probe for solver and optimizer settings. It does not establish
  statistical validity.

- Preparation and design:

  [`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md),
  [`crop_to_focal_buffer`](https://wpeterman.github.io/terradish/reference/crop_to_focal_buffer.md),
  [`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md),
  and
  [`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md)
  support reproducible inputs and model-based study design.
  [`slim_terradish`](https://wpeterman.github.io/terradish/reference/slim_terradish.md)
  reduces stored fit size when later prediction is unnecessary.

**Parameter names and their scales.** Several names recur with different
meanings. Read each one on the scale of the model that reports it.

- `theta`:

  Conductance parameters, returned by
  [`coef()`](https://rdrr.io/r/stats/coef.html). Log-linear and spline
  coefficients are on the log-conductance scale, so `exp(theta)` is a
  conductance multiplier per unit of the covariate. Gaussian
  `sigma.<layer>` entries of `theta` are smoothing scales in map units.

- `phi`:

  Nuisance parameters of the measurement model, reported by
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html).
  Their names and scales depend on the model: `alpha`, `beta`, `tau`,
  `rho`, `sigma`, and environmental coefficients.

- `tau`:

  In [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md)
  and
  [`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
  the log residual precision: the residual variance is `exp(-tau)`. In
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
  [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
  and
  [`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md),
  the nonnegative multiplier on the graph covariance, on its natural
  scale.

- `sigma`:

  In fitted Wishart models, the log nugget variance: the nugget is
  `exp(sigma)`. In
  [`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)
  and
  [`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md),
  the nugget variance itself, on its natural scale. In
  [`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
  a Gaussian smoothing scale in map units (`sigma.<layer>`). In
  [`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md),
  a line-search constant unrelated to any model.

- `rho`:

  In [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md),
  the shared-site correlation on a logit scale: the correlation is
  `plogis(rho) / 2`, between 0 and 0.5. `rho_start` in
  [`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
  is on the natural scale. In
  [`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md),
  an unrelated line-search constant.

- `nu`:

  Effective Wishart degrees of freedom, supplied by the user and never
  estimated. It scales the log-likelihood, so standard errors scale as
  `1 / sqrt(nu)` while point estimates do not change. It is not the
  marker count. Unused by the distance models.

- `alpha`, `beta`:

  In the distance models, the intercept and the slope of genetic
  distance on resistance distance, in the units of the response.

- `lambda_<covariate>`:

  In
  [`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
  a nonnegative weight on an environmental covariance kernel.

**Moving from radish.** The original `radish*` function names still work
but are deprecated; see
[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md).

**Experimental features.** Directed and hierarchical models,
site-specific drift terms, selected-pair likelihoods, and other research
prototypes remain on the `experimental` branch. Their identifiability
and validation limits are described in that branch's `EXPERIMENTAL.md`.
Install it explicitly with
`remotes::install_github("wpeterman/terradish@experimental")`.

**Vignettes.** Start with
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md),
or with
[`vignette("v03-mlpe-pairwise-distance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v03-mlpe-pairwise-distance.md)
if your data are a pairwise genetic distance matrix. The others cover
model comparison, joint isolation by environment and resistance, the
covariance-based Wishart workflow, spline conductance, scale
optimization, large landscapes, and simulation-based study design. List
them all with `browseVignettes("terradish")`.

## Worked examples

[`vignette("v12-data-generation", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v12-data-generation.md),
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md),
[`vignette("v09-individual-salamanders", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v09-individual-salamanders.md),
[`vignette("v03-mlpe-pairwise-distance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v03-mlpe-pairwise-distance.md).

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

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md),
[`melip`](https://wpeterman.github.io/terradish/reference/melip.md),
[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)

## Author

**Maintainer**: Bill Peterman <peterman.73@osu.edu>
([ORCID](https://orcid.org/0000-0001-5229-9268))

Authors:

- Bill Peterman <peterman.73@osu.edu>
  ([ORCID](https://orcid.org/0000-0001-5229-9268))

- Nathaniel S. Pope <nspope@utexas.edu>
  ([ORCID](https://orcid.org/0000-0001-8409-7812))

Other contributors:

- Denis Demidov (Bundled AMGCL library in inst/include/amgcl)
  \[copyright holder\]
