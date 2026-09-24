from pathlib import Path

p = Path('dev/check/phase7-package/vignettes/wishart-covariance.Rmd')
s = p.read_text(encoding='utf-8')

def between(start, end, replacement):
    global s
    a, b = s.index(start), s.index(end, s.index(start))
    s = s[:a] + replacement + '\n\n' + s[b:]

between('2.  Fit the resistance model', '| Model | Input', '''2. Fit the resistance model with `wishart_covariance`, using a declared
   effective degrees of freedom, `nu`.

Both Wishart interfaces use the **contrast likelihood**. For an orthonormal
site-contrast matrix `Q`, the covariance likelihood uses `Q' S Q` and the
corresponding fitted contrast covariance. Components in the constant-site
direction are excluded. Centering `S` therefore leaves this likelihood
unchanged. The previous full-rank formulation was removed because that
constant direction depends on arbitrary centering conventions and does not
provide the intended resistance information.

`dist_from_cov(S)` retains the covariance information used by this likelihood.
Covariance and matching squared-distance inputs are two representations of
the same contrast model. Neither is inherently more informative. Begin with
`vignette("getting-started", package = "terradish")` for the graph workflow.''')
s = s.replace('row/column per population) that can be passed directly to `terradish` as\nthe response `S`.', '''row/column per population). Check the construction, sampling design, site
order, and covariance geometry before using the matrix as response `S`.
The helpers are provided by the companion `landgraph` package.''')
s = s.replace('### Numeric feature matrix:', '''Biallelic standardization can give rare variants large influence. Apply and
report a defensible minor-allele-frequency filter, then examine sensitivity
to that filter. Unequal chromosome counts also produce unequal sampling
variance across populations; a common nugget does not remove that concern.
Covariance metadata describe the construction, not proof of Wishart sampling.

### Numeric feature matrix:''')
between('Wishart models treat `nu`', '### Key arguments', '''`nu` is a positive effective degrees-of-freedom value supplied by the analyst.
It is not estimated by `terradish`, and it is not fixed automatically by the
number of markers used to build `S`. Linkage creates dependence among nearby
markers; shared population history can create dependence even among unlinked
markers. Counting SNPs, loci, or within-locus allele contrasts does not by
itself establish effective replication.

A block jackknife can assess how marker blocks contribute to variability.
Because it need not capture all dependence due to shared history, regard its
implied information as an upper-bound guide rather than a complete calibration.
Describe the blocks and the uncertainty sources the procedure represents.

Declare a primary value with its rationale and a range of plausible values.
Use `terradish_rescale_nu()` for sensitivity of intervals and likelihood-based
comparisons. Use fixed-domain cross-validation to assess conductance terms
and shape. A stable result over a declared range is evidence of sensitivity
within that range, not proof that the primary value is correct.

If only a published F~ST~ matrix is available, a reported locus count does not
make that response Wishart. F~ST~ is a ratio estimator; consider MLPE for such
dissimilarities and evaluate its shared-site assumptions. See
`vignette("model-comparison", package = "terradish")` for predictive comparisons.''')
s = s.replace('`"within"` replaces the Gower diagonal with within-population variance (Dyer-style); `"gower"` leaves it as-is', '`"auto"` selects `"gower"`; `"within"` replaces the diagonal with within-population variance')
between('The `diagonal = "within"` setting', '------------------------------------------------------------------------\n\n## Set up the landscape', '''The default grouped-data construction uses `diagonal = "gower"`. Replacing
that diagonal with `diagonal = "within"` changes the covariance estimand and
can change fitted signs; the supplied release simulations recovered the
expected signs with Gower diagonals but not with the within-population
replacement. Retain `within` only for a separately justified construction,
and compare its effect explicitly. It may also produce an indefinite matrix.

```{r check-eigenvalues}
ev <- eigen(S_ms, symmetric = TRUE, only.values = TRUE)$values
range(ev)
```

**How to read the output.** A centered covariance has a zero eigenvalue in
the constant direction, so strict positive definiteness of the full matrix
is not the target. Tiny negative values may reflect rounding. Substantive
negative contrast eigenvalues require investigation; changing the diagonal
is not a neutral repair.''')
s = s.replace('# nu = 1000 is an illustrative effective marker count for this simulation.', '# nu = 1000 is an illustrative simulation parameter, not a marker-count rule.')
s = s.replace('# This is what would come from cov_from_genetic_data() applied to real genotypes.', '# This draw follows the fitted statistical model; real genotype summaries may not.')
s = s.replace('The fitted covariance model is: \\> **Sigma = tau · E(theta) + exp(sigma)\n· I**', '''The covariance kernel is `Sigma = tau * E(theta) + exp(sigma) * I`.
The likelihood evaluates its projection into the site-contrast subspace.''')
a = s.index('fit_nu <- function(nu) {')
b = s.index('# Point estimates agree', a)
s = s[:a] + '''# Reuse the fitted optimum; rescale its likelihood and uncertainty.
fit_nu <- function(nu) {
  f <- terradish_rescale_nu(fit_wc, nu = nu)
  list(theta = coef(f), se = sqrt(diag(vcov(f))))
}
f100   <- fit_nu(100)
f1000  <- fit_nu(1000)
f10000 <- fit_nu(10000)

''' + s[b:]
s = s.replace('# Point estimates agree across three orders of magnitude', '# Point estimates agree across this 100-fold range')
s = s.replace('standard error scales with 1 / sqrt(`nu`). Under the model, more', 'interior standard error scales with 1 / sqrt(`nu`). Boundary inference\nrequires separate care. Under the model, more')
between('The practical upshot:', '------------------------------------------------------------------------\n\n## Recover the true parameters', '''Report how intervals and model comparisons depend on `nu`. Rescaling
does not re-estimate conductance, so this sensitivity is inexpensive and
separates the fitted surface from the assumed precision. Use the same `nu`
for every candidate in a comparison. Neither a high AIC weight nor a small
nominal LRT p-value establishes that `nu` represents independent information.''')
s = s.replace('entire Wishart negative log-likelihood is multiplied by `nu / 2`', 'parameter-dependent Wishart objective is multiplied by `nu / 2`')
s = s.replace('```{r plot-fit, fig.cap = "***Observed vs. fitted genetic covariance.***"}', '```{r plot-fit, fig.cap = "***Observed and fitted covariance in the site-centered representation.***"}')
s = s.replace('in `S_cov`. Points near the diagonal indicate a good fit.', 'in the centered representation used by the contrast fit. Points near the\ndiagonal indicate agreement for that summary, not calibrated sampling uncertainty.')
s = s.replace('<!-- Phase 7: rewrite -->\n\n', '')
s = s.replace('Check the realized matrix or use a Gaussian distance likelihood such as `mlpe`.', 'F~ST~ is a ratio estimator without a Wishart sampling justification; passing a geometry check is insufficient. Consider `mlpe` and assess its assumptions.')
between('| Consideration | `wishart_covariance`', '------------------------------------------------------------------------\n\n## Quick-reference', '''| Consideration | Guidance |
|----|----|
| Raw genotype data available | Construct and document an appropriate covariance summary; assess filtering and unequal sampling. |
| Covariance and its matching squared distances available | Both Wishart interfaces use the same contrast likelihood. Choose the representation that makes preprocessing clearest. |
| Only a precomputed dissimilarity is available | Geometry is necessary but insufficient for Wishart sampling; consider MLPE for generic genetic distances. |
| Site-specific sampling variance matters | A common nugget does not fully address unequal sampling; examine design and sensitivity. |
| Effective information is uncertain | Declare a range of `nu`, rescale inference, and assess conductance complexity with held-out prediction. |

### Environmental differences between sites

`pairwise_covariates()` combines geographic distances and endpoint environmental
differences in one design. Both `mlpe_covariates()` and `wishart_covariates()`
accept that design, while linking it to different response models. For Wishart,
the nonnegative environmental kernel weights describe conditional contributions
to covariance structure. They do not identify a particular biological
mechanism. See `vignette("ibe-ibr-workflow", package = "terradish")` for a
worked comparison, conditional conductance estimates, and IBE:IBR ratios.''')
s = s.replace('  adding a diagnostic residual field and evaluating sensitivity to\n  spatial confounding.\n', '')
p.write_text(s, encoding='utf-8')
