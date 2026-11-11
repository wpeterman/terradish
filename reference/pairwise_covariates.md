# Combine environmental differences with other pairwise covariates

Combines differences calculated from site attributes with user-supplied
pairwise matrices, such as geographic distance, in one measurement-model
input. Values retain their supplied units and scaling.

## Usage

``` r
pairwise_covariates(...)
```

## Arguments

- ...:

  Objects from
  [`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
  or named symmetric numeric matrices or `dist` objects. Matrices must
  have a zero diagonal. Every input must describe the same sites in the
  same order.

## Value

A classed numeric matrix with one row per unordered site pair and one
column per covariate. Rows follow the lower-triangle order used by R's
[`dist()`](https://rdrr.io/r/stats/dist.html). Both
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
and
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)
accept this object.

## Details

No centering, scaling, or geometric correction is applied. For a Wishart
model, each supplied distance must additionally produce a positive
semidefinite centered kernel; the measurement-model constructor checks
it. Site labels, when supplied on square matrices, must agree across
inputs.

## Interpreting pairwise effects

Construct one `pairwise_covariates` object to use the same environmental
transforms and geographic distances in both likelihood families. Their
fitted conductance surfaces are conditional on these terms, and
conductance signs can change when the mean structure changes. Report
fits with and without pairwise effects. A uniform `~ 1` conductance
model with pairwise terms represents IBD plus IBE, not IBE alone.

A positive environmental coefficient means differentiation increases
with environmental difference conditional on resistance. It does not
distinguish assortative mating, dispersal filtering, local adaptation,
or site variance that rises at environmental extremes. MLPE coefficients
are unconstrained; Wishart kernel weights are nonnegative. Compare their
effects in resistance-distance units with
[`terradish_ibe_ratio`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
which divides by the resistance coefficient and propagates joint
uncertainty. Ratios are undefined when that denominator is zero.

## See also

[`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

## Examples

``` r
climate <- pairwise_endpoint_covariates(data.frame(temp = c(4, 7, 9)))
geography <- dist(cbind(c(0, 1, 3), c(0, 2, 1)))
z <- pairwise_covariates(climate, geographic = geography)
colnames(z)
#> [1] "absdiff_temp" "geographic"  
```
