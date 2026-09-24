from pathlib import Path

p = Path('dev/check/phase7-package/vignettes/gaussian-scale-optimization.Rmd')
s = p.read_text(encoding='utf-8')
s = s.replace('$x_{i,\\sigma}$ is the Gaussian-smoothed raster value at cell *i*', '$x_{i,\\sigma}$ is the Gaussian-smoothed raster value at cell *i*,\n  standardized over retained cells after smoothing')
s = s.replace('Large `sigma` means that it is averaged over a broader neighborhood.', '''Large `sigma` means that it is averaged over a broader neighborhood. The
coefficient `beta` describes a one-standard-deviation change in that smoothed
layer, with standardization recomputed at each trial scale. Its magnitude
therefore needs to be interpreted jointly with `sigma`.''')
s = s.replace('spatial scale of a causal ecological process.', '''spatial scale of a causal ecological process. Dispersal can average exposure
over a neighborhood, so a fitted positive scale may reflect that averaging
even when the underlying ecological response is local.''', 1)
s = s.replace('- the upper bound is the diagonal length of the retained raster extent', '''- the upper bound is `min(nrow(raster), ncol(raster)) * mean(res(raster)) / 6`,
  approximately one-sixth of the smaller raster dimension for square cells''')
s = s.replace('you create the conductance model.', '''you create the conductance model. The kernel is truncated to a raster-sized
window. The default upper bound keeps approximately three standard deviations
inside half the smaller dimension; larger requested bounds issue a truncation
warning. Near that upper bound the covariate can become almost constant, so
scale and effect size may be weakly identified. The half-cell lower bound
still smooths and is not equivalent to an unsmoothed raster.''', 1)
s = s.replace('maxit = 12', 'maxit = 100')
s = s.replace('fit_gaussian\n```', 'fit_gaussian\nfit_gaussian$convergence\n```', 1)
s = s.replace('```{r coef}\ncoef(fit_gaussian)', '```{r coef}\ncoef(fit_gaussian)\nsummary(fit_gaussian)$sigma_table\nconfint(fit_gaussian)')
anchor = '`gaussian_scale_summary()` converts'
s = s.replace(anchor, '''**How to read the output.** The scale table includes its estimate, standard
error, and proximity to a fitting bound. Wald intervals from `confint()` are
limited to those bounds; a bound-limited interval does not establish an
interior optimum. Inspect the likelihood profile when scale is important:

```{r sigma-profile, eval = FALSE}
# Refit the other parameters at a grid of fixed forest-cover scales.
profile <- gaussian_scale_profile(fit_gaussian, "forestcover", n = 15)
profile$interval
confint(profile$fit)
```

The returned fit stores the profile interval. Increase grid resolution when
the profile is irregular. Failed profile fits stop the calculation rather
than being silently omitted. Boundary likelihood-ratio calibration remains
approximate; report a bound-limited interval as such.

''' + anchor, 1)
s = s.replace('Gaussian kernel mass, while', 'untruncated Gaussian kernel mass, while')
s = s.replace('the smoothing kernel is relative to the raster resolution.', '''the smoothing kernel is relative to the raster resolution. These ideal
Gaussian radii are descriptive; raster truncation and edge handling can
change the realized weights.''')
s = s.replace('```{r sigma-plot-km, fig.cap', '```{r sigma-plot-km, eval = FALSE, fig.cap')
s = s.replace('in kilometers:', '''in kilometers. This conversion is only valid for a meter-based raster;
do not apply it to longitude/latitude degrees:''')
a = s.index('At this stage, the main question is not simply')
b = s.index('### Step 6:', a)
s = s[:a] + '''Assess whether the scale is identified away from its bounds, whether the
profile supports a useful range, and whether estimating it improves held-out
prediction. With only eight sites this illustration cannot provide a strong
spatial-transfer assessment. On a larger sampling design, use common folds
and a distinct conductance factory for each candidate:

```{r gaussian-cv, eval = FALSE}
coords_projected <- terra::project(coords, "EPSG:5070")
folds <- terradish_folds(coords_projected, k = 2, seed = 42)
cv_scale <- terradish_cv_folds(
  data = surface,
  formulas = list(fixed = melip.Fst ~ forestcover,
                  gaussian = melip.Fst ~ forestcover),
  folds = folds, model = mlpe, nuisance = "fixed",
  conductance_model = list(fixed = loglinear_conductance,
                           gaussian = gaussian_model),
  control = NewtonRaphsonControl(maxit = 100)
)
cv_scale$summary
```

Inspect failures and paired differences on common successful folds. A scale
that improves in-sample fit but not prediction needs cautious interpretation.

''' + s[b:]
s = s.replace('That does **not** mean the model failed. It means the scale parameters are\ntelling you that the data may prefer a simpler spatial structure.', '''Check convergence separately from boundary status. A converged boundary fit
can indicate limited scale information or preference for an extreme allowed
scale. A nonconverged fit requires numerical investigation before ecological
interpretation.''')
anchor = '## Important current limitations'
s = s.replace(anchor, '''## Outer scale searches

`terradish_scale_optim()` remains available for an outer search, while joint
Gaussian fitting provides derivatives and joint scale inference directly.
Outer-search radii are in map units and the default outer kernel is Gaussian.
A Gaussian standard deviation, a uniform-disk radius, and other kernel
parameters represent different weighting functions even when their numeric
values match. Record the kernel, units, search bounds, and estimated scales.
The parameter count includes optimized scales; repeated searching is still
part of the declared model-selection procedure. Prefer projected rasters
when distances are to be interpreted in meters or kilometers.

''' + anchor, 1)
s = s.replace('   units, and defaults to the retained raster extent\'s diagonal when left NULL.', '   units, and defaults to the smaller raster dimension / 6 for square cells.')
s = s.replace('melip.Fst ~ altitude, data = surface', 'melip.Fst ~ forestcover, data = surface')
s = s.replace('| `conductance()` | Extract', '| `gaussian_scale_profile()` | Profile scale uncertainty while reoptimizing other parameters |\n| `terradish_cv_folds()` | Compare fixed and estimated scales on common held-out sites |\n| `conductance()` | Extract')
p.write_text(s, encoding='utf-8')
