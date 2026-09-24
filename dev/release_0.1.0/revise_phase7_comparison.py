from pathlib import Path

p = Path('dev/check/phase7-package/vignettes/model-comparison.Rmd')
s = p.read_text(encoding='utf-8')
s = s.replace('equal conductance. Resistance distance then equals graph distance\n(weighted only by geography, not landscape features). This is the IBD\nbaseline:', '''equal conductance. Resistance still combines all routes through the graph;
it is neither shortest-path distance nor necessarily proportional to Euclidean
distance. This uniform-surface model is the graph-based IBD baseline:''')
s = s.replace('> chi-squared reference distribution.', '''> chi-squared reference distribution. The response, site order, graph,
> pairwise design, basis definition, and supplied Wishart `nu` must agree.
> A linear term is not automatically nested in an arbitrary spline basis, and
> a fitted Gaussian scale changes the model dimension. For Wishart fits,
> the LRT magnitude depends on `nu`; interpret it descriptively at the declared
> value and report sensitivity rather than treating nominal marker counts as
> calibrated replication.''')
s = s.replace('criterion instead.', 'criterion or held-out prediction instead.', 1)
s = s.replace('### When to use AICc', '''### Effective information and over-selection

Wishart AIC can select unnecessary conductance terms when nominal `nu`
overstates effective information. AICc does not repair an incorrectly scaled
likelihood. Use the same predeclared candidates in cross-validation and
compare inference over plausible effective degrees of freedom:

```{r nu-sensitivity, eval = FALSE}
# fit_wishart is a converged Wishart fit at the declared primary nu.
sensitivity <- lapply(c(10, 30, 100), function(value) {
  terradish_rescale_nu(fit_wishart, nu = value)
})
lapply(sensitivity, confint)
```

**How to read the output.** Rescaling changes uncertainty and likelihood-based
comparisons without changing the fitted coefficients. If conclusions depend
on a high `nu`, report that dependence. Rescale every candidate to the same
value before comparing their information criteria.

### When to use AICc''')
s = s.replace('simpler models. The implementation here uses the number of selected pair\nrows for a pair-subset fit and otherwise uses the number of pairwise\nobservations', 'simpler models. The implementation uses the number of pairwise\nobservations')
s = s.replace('focal-support clamping to keep model comparison focused on supported\ncovariate ranges:', 'focal-support clamping to inspect predictions within sampled covariate\nranges. This display choice does not alter the fitted likelihood or model\ncomparison:')
s = s.replace('<!-- Phase 7: rewrite -->\n\n', '')
s = s.replace('  model = mlpe,\n  keep_fits', '  model = mlpe,\n  nuisance = "fixed",\n  keep_fits')
a = s.index('For each test fold, the function reprofiles')
b = s.index('The comparison boundary remains strict.', a)
s = s[:a] + '''With `nuisance = "fixed"`, all measurement-model coefficients are learned
from training sites and used unchanged to score held-out sites. This tests
prediction from the trained model. A separately fitted uniform-conductance
baseline provides `loglik_gain`; positive gain favors the trained landscape
surface over that baseline on the same held-out response.

The default `nuisance = "reprofile"` instead optimizes nuisance coefficients
using each test response, separately for the trained surface and baseline.
That asks whether the learned surface retains useful structure after local
recalibration. It is a different prediction target. Declare the mode and do
not compare scores across modes.

### Different models and repeated folds

Each formula can have its own conductance or measurement-model factory. For
example, compare a linear forest-cover term with a fixed-degree spline using
the same response and folds:

```{r cv-repeats, eval = FALSE}
folds_repeated <- lapply(c(42, 91, 137), function(seed) {
  terradish_folds(coords_projected, k = 5, seed = seed)
})
cv_shape <- terradish_cv_folds(
  data = surface,
  formulas = list(linear = melip.Fst ~ forestcover,
                  spline = melip.Fst ~ s(forestcover, df = 3)),
  folds = folds_repeated,
  model = mlpe,
  conductance_model = list(linear = loglinear_conductance,
                           spline = smooth_loglinear_conductance),
  nuisance = "fixed"
)
cv_shape$summary
cv_shape$results
```

**How to read the output.** Inspect convergence and score failures first.
Rankings use folds on which every candidate succeeded. `mean_difference`
and `se_difference` describe paired score differences against the best
candidate on those common folds; `within_one_se` flags a difference no larger
than its estimated standard error. Prefer a simpler model when its predictive
performance is indistinguishable under the declared design. Overlapping
training sets and repeated partitions make these standard errors descriptive,
not independent-replicate hypothesis tests. Repeats expose split sensitivity;
they do not create new independent sites.

Checkpoint files retain completed fits and scores. Resume only the identical
response, graph, formulas, model settings, and folds; changed inputs are
rejected. Keep failed baseline scores separate from failed trained-model scores:
a valid held-out score can exist even when its baseline gain is unavailable.

''' + s[b:]
a = s.index('Concordant LRT, information-criterion, and held-out results')
b = s.index('------------------------------------------------------------------------', a)
s = s[:a] + '''Use cross-validation to assess conductance terms and shape, then use `nu`
sensitivity to describe the uncertainty of Wishart inference. Treat an LRT at
a fixed `nu` as conditional, descriptive evidence. Information criteria add
a view of in-sample fit and complexity, but agreement among these summaries
does not validate the underlying effective sample size. None establishes
causality, observed movement, absolute permeability, or unique attribution to
a landscape mechanism.

''' + s[b:]
s = s.replace('  model = mlpe\n)', '  model = mlpe, nuisance = "fixed"\n)')
s = s.replace('cv$summary                        # compare within this likelihood and fold set', 'cv$summary                        # paired differences on common successful folds')
p.write_text(s, encoding='utf-8')
