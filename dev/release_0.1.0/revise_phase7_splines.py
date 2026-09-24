from pathlib import Path

p = Path('dev/check/phase7-package/vignettes/spline-conductance.Rmd')
s = p.read_text(encoding='utf-8')
s = s.replace('which is computed once and then treated as fixed.', '''which is computed once and then treated as fixed. These are **unpenalized**
regression splines with a chosen number of coefficients, not smooths whose
complexity is estimated by a smoothing penalty. Knots, boundary knots, and
column-centering offsets are retained from training. Plotting and new-surface
prediction reuse that basis rather than selecting new knots. Centering removes
an arbitrary constant conductance multiplier; it does not change the fitted
shape.''')
s = s.replace('B_ns <- ns(x_seq, df = 4)', 'B_ns <- scale(ns(x_seq, df = 4), center = TRUE, scale = FALSE)')
s = s.replace('B_bs <- bs(x_seq, df = 4, degree = 3)', 'B_bs <- scale(bs(x_seq, df = 4, degree = 3), center = TRUE, scale = FALSE)')
s = s.replace('The figure below shows what natural spline', 'The figure below shows centered natural spline')
s = s.replace('compare against simpler models with AICc.', 'compare against simpler predeclared models with held-out prediction.')
s = s.replace('- Consider **AICc** when the number of focal sites is small relative to the number of estimated parameters. `terradish` uses focal sites, not the number of dependent population pairs, as `n` in its AICc calculation.', '''- Use common cross-validation folds for the candidate shapes. AICc is a
  secondary diagnostic under its focal-site sample-size convention; it does
  not make all site pairs independent or calibrate an uncertain Wishart `nu`.''')
s = s.replace('Use AICc (small-sample-corrected AIC) to account for the different number of parameters across models:', '''Predeclare a small candidate set and compare held-out scores on common folds.
AIC and AICc can over-select flexible conductance models, especially when
Wishart `nu` overstates effective information. The following table remains a
useful in-sample diagnostic, but it should not determine spline complexity:''')
anchor = '## Marginal conductance associations'
new = '''### Cross-validation and monotonicity

```{r spline-cv, eval = FALSE}
coords_projected <- terra::project(melip.coords, "EPSG:5070")
folds <- terradish_folds(coords_projected, k = 5, seed = 42)
cv_shape <- terradish_cv_folds(
  data = surface,
  formulas = list(linear = melip.Fst ~ forestcover + altitude,
    spline = melip.Fst ~ forestcover + s(altitude, df = 4)),
  folds = folds, model = mlpe, nuisance = "fixed",
  conductance_model = list(linear = loglinear_conductance,
                           spline = smooth_loglinear_conductance)
)
cv_shape$summary
```

**How to read the output.** Compare paired score differences on common
successful folds. Prefer the simpler shape when the difference is small
relative to its fold-to-fold uncertainty, and inspect failures before ranking.
See `vignette("model-comparison", package = "terradish")` for repeated folds
and the distinction between fixed and reprofiled nuisance parameters.

```{r monotonicity}
summary(fit_ns_alt)$spline_monotonicity
```

The monotonicity table checks the derivative of each fitted smooth over its
focal-site covariate range. It reports direction and derivative sign changes;
it is a shape diagnostic, not a hypothesis test or evidence of an ecological
threshold. Even a predictively useful spline can reflect how a population
process maps onto the graph rather than a direct ecological response to the
covariate. Inspect the conditional curve and its uncertainty.

------------------------------------------------------------------------

'''
s = s.replace(anchor, new + anchor, 1)
s = s.replace('```{r spline-support-clamp, eval = FALSE}', '```{r spline-support-clamp}')
anchor = '## Observed vs. fitted genetic distances'
s = s.replace(anchor, '''Clamping replaces covariate values outside the requested focal-cell
quantiles with the nearest limit. The fitted spline basis and centering remain
fixed. Unclamped `conductance()` and marginal plots use the same centered
prediction basis, so they should agree on the relative conductance scale.
Clamping changes predictions, not the likelihood or the estimated shape.

''' + anchor, 1)
s = s.replace('**Reading the df comparison**: Start from the most constrained model (df = 2) and add flexibility only if the AICc improves meaningfully *and* the shape is ecologically plausible. Oscillating or extreme-at-boundary curves with high df are likely fitting noise.', '''**Reading the df comparison.** Compare these predeclared shapes using the
same held-out folds. An AICc improvement alone does not justify extra
coefficients. Inspect monotonicity, uncertainty, boundary behavior, and the
stability of the selected shape across folds.''')
s = s.replace('Evaluate stability, focal-site support, AICc under its stated convention, and simulation under plausible scenarios.', 'Evaluate held-out prediction, stability, focal-site support, and simulation under plausible scenarios.')
s = s.replace('The spline and linear models give similar AICc; parsimony favors the simpler model.', 'The spline and linear models give similar held-out scores; parsimony favors the simpler model.')
a = s.index('# Step 4: Compare with AICc')
b = s.index('# Step 7: Visualize the selected model', a)
s = s[:a] + '''# Step 4: Compare the predeclared shapes on the same folds.
folds <- terradish_folds(coords_projected, k = 5, seed = 42)
cv <- terradish_cv_folds(
  surface,
  formulas = list(linear = S ~ covar2 + covar1,
                  spline = S ~ covar2 + s(covar1, df = 3)),
  folds = folds, model = mlpe, nuisance = "fixed",
  conductance_model = list(linear = loglinear_conductance,
                           spline = smooth_loglinear_conductance)
)
cv$summary

# Step 5: Inspect the fitted shape and its derivative sign changes.
summary(fit_spline)$spline_monotonicity
plot(fit_spline, type = "marginal", data = surface)

# Step 6: Prefer the simpler model if predictive differences are inconclusive.
# Add larger-df models only as a declared sensitivity analysis, not an AIC search.

''' + s[b:]
s = s.replace('best_fit <- fit_spline  # or whichever won', 'best_fit <- fit_spline  # choose after reviewing CV, convergence, and shape')
p.write_text(s, encoding='utf-8')
