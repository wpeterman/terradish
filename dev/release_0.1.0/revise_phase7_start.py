from pathlib import Path

root = Path('dev/check/phase7-package/vignettes')
p = root / 'large-landscapes.Rmd'
s = p.read_text(encoding='utf-8').replace('assess_terradish_settings()', 'terradish_assess_settings()').replace('resistance_distances()', 'terradish_distance()')
p.write_text(s, encoding='utf-8')
p = root / 'getting-started.Rmd'
s = p.read_text(encoding='utf-8')

def between(start, end, replacement):
    global s
    a, b = s.index(start), s.index(end, s.index(start))
    s = s[:a] + replacement + '\n\n' + s[b:]

s = s.replace('`terradish` automates these steps and provides gradients, standard', '''Version 0.1.0 focuses on symmetric resistance models with log-linear,
fixed-degree spline, or jointly estimated Gaussian-scale conductance.
Environmental differences between sites can enter either MLPE or Wishart
measurement models. Directional, selected-pair, hierarchical, and other
experimental extensions are maintained on the `experimental` branch.

`terradish` automates these steps and provides gradients, standard''')
between('The **measurement model** describes', '------------------------------------------------------------------------\n\n## Installation', '''The **measurement model** links the graph to the observed genetic summary.
For pairwise genetic distances, `mlpe` accounts for a shared-site component of
dependence; `leastsquares` assumes independent errors. For a suitably constructed
covariance matrix or its corresponding squared-distance matrix, the Wishart
models use a likelihood on site contrasts. The response construction and
inferential assumptions determine which family is appropriate. A Wishart fit
also requires an analyst-supplied effective degrees of freedom, `nu`.''')
s = s.replace('`directions = 8` is typical for landscape genetics;', '`directions = 8` is the default;')
between('`mlpe` corrects for the non-independence', '### Fitting with `generalized_wishart`', '''`mlpe` models a shared-site component of dependence because pairs (A-B) and
(A-C) both involve site A. It uses maximum likelihood, not restricted maximum
likelihood. This correction does not guarantee calibrated uncertainty for every
population process. Treat `leastsquares` as a simple comparator rather than a
general solution to dependence.

`generalized_wishart` requires a squared-distance response coherently related to
a positive-semidefinite covariance matrix. Use `check_distance_response()` on
the exact response. Passing that check establishes compatible geometry, not a
valid sampling model. In particular, F~ST~ is a ratio estimator; it does not
acquire a Wishart sampling distribution merely by passing the geometry check.
For a generic genetic dissimilarity, consider MLPE and assess its assumptions.

`wishart_covariance` accepts the covariance directly, for example from
`cov_from_biallelic()`. Both Wishart interfaces evaluate the same contrast
likelihood on the site-contrast subspace when their inputs are related by
`dist_from_cov()`. Centering the covariance does not discard information used
by that likelihood. The previous full-rank covariance formulation is removed.

### What nu means

`nu` is a positive effective degrees-of-freedom value supplied by the analyst.
The fit does not estimate it. Raw SNP counts, locus counts, and allele contrasts
are not automatically effective independent observations: linkage and shared
population history can reduce the information they contain. A block jackknife
can diagnose marker dependence, but it does not capture every source of shared
history and should be treated as an upper-bound guide to precision.

For a fixed response and model, increasing `nu` leaves the theoretical optimum
unchanged while reducing interior standard errors in proportion to
`1 / sqrt(nu)`. It also magnifies likelihood differences and can favor flexible
models. State the chosen value and its rationale, compare plausible values
with `terradish_rescale_nu()`, and use held-out prediction to assess conductance
complexity. See `vignette("wishart-covariance", package = "terradish")` and
`vignette("model-comparison", package = "terradish")`.''')
s = s.replace('on the graph covariance kernel) and `sigma` (a log-scale identity', 'on the graph covariance kernel) and `sigma` (a log-scale identity')
anchor = '## Visualize the fit'
addition = '''### Uncertainty and convergence

```{r coefficient-uncertainty}
vcov(fit_mlpe)
confint(fit_mlpe)
fit_mlpe$convergence
fit_mlpe$fit$subproblem$convergence
```

**How to read the output.** `vcov()` gives the covariance of the conductance
coefficient estimates; `confint()` reports intervals on their fitted scale.
For a log-linear coefficient, exponentiate both interval endpoints to obtain
an interval for the relative conductance ratio. These are conditional model
intervals, not direct intervals for migration or movement rates. At a
no-structure boundary, conductance coefficients are unidentified and their
inference is unavailable.

### Convergence and starting values

Check convergence before interpreting estimates. Code 0 requires either a
projected gradient below `ctol`, or an objective change below `ftol` together
with a projected gradient below `sqrt(ctol)`. The projected gradient omits
components pointing outside active parameter bounds. Code 1 means the
iteration limit was reached; code 2 means a stall or failed line search. A
small likelihood change alone is insufficient. Also inspect the nested
nuisance fit's convergence code.

If a user-supplied starting point reaches a no-structure boundary, the fitter
tries the default start and retains the better likelihood. The convergence
record reports whether that restart occurred. For difficult fits, compare
scientifically plausible starting values and inspect the final likelihoods,
gradients, coefficients, and active bounds. Increasing the iteration limit can
help a progressing fit; weakening tolerances does not resolve a stalled one.
Numerical convergence alone does not establish model adequacy.

------------------------------------------------------------------------

'''
s = s.replace(anchor, addition + anchor, 1)
between('Use `support = "none"` (default)', '------------------------------------------------------------------------\n\n## Quick-reference', '''Focal support means the covariate values in the raster cells containing
sampled sites. The default `support_probs = c(0, 1)` uses their minimum and
maximum; the example uses their 1st and 99th percentiles. Values outside those
limits are set to the nearest limit before conductance is evaluated. For
splines, evaluation retains the fitted knots and basis definition.

Use `support = "none"` (default) for full-domain evaluation. Clamping with
`support = "focal"` changes the displayed or predicted surface; it does not
refit the likelihood, change the fitted graph, or remove uncertainty about
poorly sampled environments. Report the limits and compare both maps.
Predictions on new rasters must retain training transformations; use
`scale_covariates(new_covariates, reference = covariates)` before building the
new surface rather than scaling each domain independently.''')
s = s.replace('summary(fit)                                   ', 'summary(fit)                                   ')
p.write_text(s, encoding='utf-8')
