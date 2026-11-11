# Wishart covariance measurement model

A function of class `"terradish_measurement_model"` that evaluates a
Wishart likelihood for an observed genetic **covariance** matrix. This
is the natural model when raw genotype data are summarized with
[`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html)
or
[`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html)
rather than as pairwise distances. For an admissible squared-distance
representation coherently related to a centered positive semidefinite
covariance matrix, use
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md).
For other pairwise distances, use
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) or
[`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md).

## Usage

``` r
wishart_covariance(
  E,
  S,
  phi,
  nu,
  gradient = TRUE,
  hessian = TRUE,
  partial = TRUE,
  nonnegative = TRUE,
  validate = FALSE
)
```

## Arguments

- E:

  Conductance-implied covariance matrix: the generalized inverse of the
  graph Laplacian at the current conductance parameters. Passed
  automatically by the optimizer; users normally do not call this
  function directly.

- S:

  Square, symmetric matrix of observed population-level genetic
  covariance. Typically the output of
  [`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html)
  or
  [`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html).
  Must have the same dimensions as `E` and should be positive
  (semi-)definite.

- phi:

  Named numeric vector of nuisance parameters `(tau, sigma)`. Omit to
  obtain least-squares starting values.

- nu:

  One finite positive number describing effective Wishart information.
  You must supply it; the fit does not estimate it. It is not the marker
  count. See the effective-information section and report sensitivity
  across plausible values with
  [`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md).

- gradient:

  Logical. Compute gradient of the negative log-likelihood with respect
  to `phi`?

- hessian:

  Logical. Compute Hessian with respect to `phi`?

- partial:

  Logical. Compute second partial derivatives with respect to `phi`,
  `E`, and `S`?

- nonnegative:

  Unused; present for interface consistency.

- validate:

  Logical. Numerically validate gradients and Hessians via numDeriv?
  Very slow; for debugging small examples only.

## Value

When `phi` is missing, a list with elements `phi` (least-squares
starting values for `(tau, sigma)`), `lower` (`c(0, -Inf)`), and
`upper`. Otherwise a list containing:

- objective:

  Negative log-likelihood.

- fitted:

  Fitted covariance matrix \\\Sigma\\ (same dimensions as `S`).

- boundary:

  Logical; `TRUE` if `tau = 0`, indicating no detectable IBR signal.

- gradient:

  Gradient with respect to `phi` (if `gradient = TRUE`).

- hessian:

  Hessian matrix with respect to `phi` (if `hessian = TRUE`).

- gradient_E:

  Gradient with respect to `E` (if `partial = TRUE`).

- partial_E:

  Jacobian of `gradient_E` with respect to `phi` (if `partial = TRUE`).

- partial_S:

  Jacobian of `gradient` with respect to the lower triangle of `S` (if
  `partial = TRUE`).

- jacobian_E:

  Function for reverse-mode AD through `E` (if `partial = TRUE`).

- jacobian_S:

  Function for reverse-mode AD through `S` (if `partial = TRUE`).

## Details

The nuisance parameters are:

- `tau`:

  Nonnegative scale applied to the conductance-implied covariance `E`. A
  value near zero signals no detectable graph-kernel contribution under
  this response and model.

- `sigma`:

  Log-scale nugget: the identity component added to the model covariance
  is \\\exp(\sigma)\\. It captures diagonal variation not represented by
  the graph covariance, including sampling variation and model mismatch.
  It is not a smoothing scale or a direct estimate of genetic drift or
  population size.

The fitted covariance is \\\Sigma = \tau E + \exp(\sigma) I\\. The
likelihood is evaluated on orthonormal site contrasts. Let \\L\\ have
\\n-1\\ orthonormal columns perpendicular to the constant vector, and
let \\\tilde\Sigma = L^T\Sigma L\\ and \\\tilde S = L^T S L\\. The
negative log-likelihood, omitting terms constant in the fitted
parameters for a fixed response and supplied \\\nu\\, is:

\$\$-\ell = \frac{\nu}{2} \left\[ \log\|\tilde\Sigma\| +
\mathrm{tr}(\tilde\Sigma^{-1} \tilde S) \right\]\$\$

The response is a sample covariance on site contrasts with effective
degrees of freedom \\\nu\\. Adding a constant to every entry of `S` or
centering it across sites does not change this objective. This is the
same objective as `generalized_wishart(dist_from_cov(S))`; the nugget
contributes \\2\exp(\sigma)\\ to the expected squared distances. Working
on contrasts discards only the overall site-mean direction, which
carries no information in a centered genetic covariance.

**Preparing the response.** Use a coherent covariance construction, such
as `diagonal = "gower"` for grouped data in landgraph. A `"within"`
diagonal is on a different scale and is unsuitable for Wishart fitting,
even when the matrix happens to be positive definite. Per-locus
standardization up-weights rare variants; filter by minor-allele
frequency and assess sensitivity. Unequal group sizes create unequal
sampling variance that a common diagonal nugget does not automatically
fix.

During fitting,
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
checks very small nugget estimates for false convergence on the log
scale. Improving finite-variance starts retain the supplied bounds. An
unrestricted nugget can instead be represented at a verified numerical
zero limit: its finite log value gives a positive variance below
floating-point covariance precision, and the one-sided variance gradient
must satisfy the stopping tolerance. This does not change the likelihood
or override fixed parameters. The fitted object records
`fit$residual_variance_zero_limit`;
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
omits Wald uncertainty for `sigma` at that limit.

## Worked examples

[`vignette("v02-wishart-covariance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v02-wishart-covariance.md),
[`vignette("v09-individual-salamanders", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v09-individual-salamanders.md).

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

Hanks EM, Hooten MB (2013). Circuit theory and model-based inference for
landscape connectivity. Journal of the American Statistical Association,
108, 22–33.
[doi:10.1080/01621459.2012.724647](https://doi.org/10.1080/01621459.2012.724647)
.

## See also

[`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html),
[`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html),
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
[`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

## Examples

``` r
E <- matrix(c(1.2, 0.3, 0.2,
              0.3, 1.5, 0.4,
              0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
Sigma <- 0.8 * E + 0.2 * diag(3)
set.seed(1)
S <- rWishart(1, df = 20, Sigma = Sigma)[, , 1] / 20
# tau multiplies the graph covariance; sigma is the LOG nugget variance.
value <- wishart_covariance(E, S, phi = c(tau = 0.8, sigma = log(0.2)),
                            nu = 20)
value$objective  # negative log-likelihood at these values
#> [1] 24.51609
value$fitted     # fitted covariance, tau * E + exp(sigma) * I
#>      [,1] [,2] [,3]
#> [1,] 1.16 0.24 0.16
#> [2,] 0.24 1.40 0.32
#> [3,] 0.16 0.32 1.08
```
