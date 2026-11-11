# Package index

## Build the graph

Prepare raster covariates and connect sampling locations to a raster
graph.

- [`scale_covariates()`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)
  : Scale raster covariates and retain original-scale metadata
- [`conductance_surface()`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
  : Build the raster graph from covariates and sampling sites
- [`crop_to_focal_buffer()`](https://wpeterman.github.io/terradish/reference/crop_to_focal_buffer.md)
  : Crop rasters to a focal-site buffer

## Genetic responses

Construct and check the genetic response. The covariance and distance
constructors are re-exported from landgraph.

- [`reexports`](https://wpeterman.github.io/terradish/reference/reexports.md)
  [`cov_from_biallelic`](https://wpeterman.github.io/terradish/reference/reexports.md)
  [`cov_from_genetic_data`](https://wpeterman.github.io/terradish/reference/reexports.md)
  [`fst_from_biallelic`](https://wpeterman.github.io/terradish/reference/reexports.md)
  [`dist_from_cov`](https://wpeterman.github.io/terradish/reference/reexports.md)
  [`dist_from_biallelic`](https://wpeterman.github.io/terradish/reference/reexports.md)
  : Objects exported from other packages
- [`check_distance_response()`](https://wpeterman.github.io/terradish/reference/check_distance_response.md)
  : Check that a squared-distance response is admissible
- [`pca_dist()`](https://wpeterman.github.io/terradish/reference/pca_dist.md)
  : Pairwise genetic distances from PCA scores

## Fit and evaluate models

Fit a conductance model by maximum likelihood and work with the fitted
object.

- [`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md)
  : Fit a conductance model by maximum likelihood
- [`print(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`summary(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`print(`*`<summary.terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`coef(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`fitted(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`simulate(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`anova(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`logLik(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`AIC(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  [`residuals(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
  : Print, summarize, compare, and extract from a fitted model
- [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  : Predict the fitted conductance surface with intervals
- [`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md)
  : Resistance distances on a grid of conductance parameters
- [`slim_terradish()`](https://wpeterman.github.io/terradish/reference/slim_terradish.md)
  : Shrink a fitted model for storage

## Conductance models

How raster covariates map to cell conductance.

- [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md)
  : Conductance models: the factory interface
- [`loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md)
  : Log-linear conductance model
- [`smooth_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md)
  : Log-linear conductance model with spline terms
- [`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md)
  : Log-linear conductance model with estimated Gaussian smoothing
  scales

## Spatial scale

Estimate, profile, and summarize Gaussian smoothing scales.

- [`terradish_scale_optim()`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md)
  : Search raster smoothing scales by refitting the model
- [`gaussian_scale_profile()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md)
  : Profile-likelihood interval for a Gaussian smoothing scale
- [`gaussian_scale_summary()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_summary.md)
  : Summarize fitted Gaussian smoothing scales

## Measurement models

How resistance distance or covariance links to the observed genetic
response, including pairwise and site-level environmental covariates.

- [`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md)
  : Measurement models: the likelihood interface
- [`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md) :
  Maximum likelihood population effects (MLPE) measurement model
- [`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
  : MLPE measurement model with pairwise environmental covariates
- [`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md)
  : Change in genetic distance across an MLPE pairwise covariate
  contrast
- [`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
  : Wishart covariance measurement model
- [`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)
  : Wishart measurement model with environmental covariance kernels
- [`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
  : Generalized Wishart measurement model for squared distances
- [`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md)
  : Least-squares measurement model
- [`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
  : Combine environmental differences with other pairwise covariates
- [`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md)
  : Pairwise environmental differences from site-level covariates
- [`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md)
  : Express environmental effects in resistance-distance units

## Inference and model comparison

Standard errors and intervals, information criteria, and spatial
cross-validation.

- [`vcov(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
  [`nobs(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
  [`confint(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)
  : Standard errors and confidence intervals for conductance parameters
- [`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
  : Rescale the effective information in a Wishart fit
- [`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md)
  : Rank fitted models by AIC, AICc, or BIC
- [`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md)
  : Assign sampling sites to cross-validation folds
- [`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
  [`print(`*`<terradish_cv_folds>`*`)`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
  : Spatial cross-validation of conductance formulas on a fixed graph

## Plotting

Diagnostic and interpretive plots of a fitted model.

- [`plot(`*`<terradish>`*`)`](https://wpeterman.github.io/terradish/reference/plot.terradish.md)
  : Plot a fitted model: fit, surface, marginal curves, and smoothing
  kernels

## Simulation and study design

Simulate covariance responses and assess parameter recovery.

- [`simulate_covariance_response()`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)
  : Simulate genetic covariance matrices from a known conductance
  surface
- [`covariance_response_power()`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md)
  : Simulation-based power and parameter recovery for sampling designs
- [`print(`*`<terradish_covariance_power>`*`)`](https://wpeterman.github.io/terradish/reference/print.terradish_covariance_power.md)
  : Print and read a power assessment

## Large landscapes and solver settings

Solver benchmarks and speed settings for large rasters.

- [`terradish_assess_settings()`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md)
  : Recommend optimizer and solver settings for a data set
- [`print(`*`<terradish_setting_assessment>`*`)`](https://wpeterman.github.io/terradish/reference/print.terradish_setting_assessment.md)
  : Print and read a settings assessment
- [`terradish_solver_benchmark()`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)
  : Benchmark sparse Cholesky solver settings on a graph

## Optimizer controls

Tuning for the Newton-type optimizers and their line searches.

- [`NewtonRaphsonControl()`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  : Control settings for the Newton and BFGS optimizers
- [`ArmijoControl()`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md)
  : Control settings for Armijo backtracking line search
- [`HagerZhangControl()`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md)
  : Control settings for Hager-Zhang line search

## Data sets

Bundled example data used in the help pages and guides.

- [`ohio_ponds`](https://wpeterman.github.io/terradish/reference/ohio_ponds.md)
  : Simulated Ohio forest-pond SNP and microsatellite data
- [`individual_salamanders`](https://wpeterman.github.io/terradish/reference/individual_salamanders.md)
  : Simulated individual salamander genotypes and resistance analyses
- [`melip`](https://wpeterman.github.io/terradish/reference/melip.md) :
  Melipona stingless bee landscape and FST example data

## Compatibility

Names kept from the radish package and their terradish replacements.

- [`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
  : Deprecated radish names

## Low-level tools

Direct likelihood evaluation and small helpers. Most analyses do not
need these.

- [`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md)
  : Profile likelihood on a grid of conductance parameters
- [`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md)
  : Low-level likelihood and derivatives at one parameter value
- [`lower()`](https://wpeterman.github.io/terradish/reference/lower.md)
  : Extract the strict lower triangle of a matrix
- [`scale_to_0_1()`](https://wpeterman.github.io/terradish/reference/scale_to_0_1.md)
  : Rescale values to the unit interval
