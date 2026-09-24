from pathlib import Path

root = Path('dev/check/phase7-package')
p = root / 'vignettes/large-landscapes.Rmd'
s = p.read_text(encoding='utf-8')

def between(start, end, replacement):
    global s
    a, b = s.index(start), s.index(end, s.index(start))
    s = s[:a] + replacement + '\n\n' + s[b:]

between('A note on the numbers below.', '```{r load}', '''The runnable code uses the bundled `melip` data to show the workflow.
Time and memory depend on graph geometry, conductance contrast, the number of
focal sites, and hardware. Measure those costs on your own design before
committing to a large analysis.''')
between('"Large" here means', '---\n\n## Choosing a solver', '''"Large" here means that graph factorization or repeated derivative solves
dominate the analysis. There is no universal cell-count limit: missing cells,
graph connectivity, and conductance contrast can change the practical limit.''')
between('The choices:', 'This applies to every symmetric resistance model:', '''The choices are:

- **`direct`** (default) uses sparse Cholesky factorization. It is a useful
  reference for checking iterative results. The factor can contain many more
  nonzero entries than the graph, so its memory requirements can grow rapidly.
- **`amg`** uses algebraic multigrid as an iterative solver. It can reduce memory
  demands on large grids, but convergence depends on conductance contrast and
  graph structure. Strong barriers can require more iterations or tighter
  solver settings. Check agreement with a direct solve on a manageable case.
- **`auto`** uses `direct` below 1,500,000 graph vertices and `amg` at or above
  that threshold with the default controls. The number of focal sites does not
  override the large-graph switch. `solver_control$auto_amg_min_vertices`
  changes the threshold; the legacy `auto_direct_max_rhs` is ignored.

Use `assess_terradish_settings()` for a local speed and accuracy probe. Its
recommendation concerns the tested settings and graph, not the scientific
adequacy of the measurement model. Record the resolved solver and compare
coefficient estimates and likelihoods when changing solver tolerances.''')
between('On Windows, `cores > 1`', 'Use one level of parallelism', '''Worker reuse is **experimental**. On Windows, `cores > 1` requests PSOCK
workers for supported Hessian and partial-derivative calculations. It does
not parallelize every step of a fit. The pool is reused during optimization
and stopped before return; fitted objects do not retain live workers.
AMG and the cached CHOLMOD backend ignore this setting. On small graphs,
worker startup and data transfer can make multiple cores slower than one.
Benchmark `cores = 1` first, then compare the elapsed time and numerical
results with the requested worker count.''')
between('On the development workstation,', '## Storing large fitted objects', '''**How to read the benchmark.** Compare coefficients and likelihoods first,
then standard errors, Hessian eigenvalues, and elapsed time. Agreement in
point estimates does not guarantee agreement in curvature. A large difference
in uncertainty warrants investigation of residuals, conditioning, and model
misspecification before choosing the faster option.

## Information and effective degrees of freedom

For Wishart fits, `nu` is a positive effective degrees-of-freedom value supplied
by the analyst. It is not estimated from the response and is not automatically
the marker count. Linkage and shared population history can reduce the
information represented by many markers. More graph cells alone do not supply
more independent genetic information.

Holding the response and model fixed, changing `nu` leaves the likelihood
optimum unchanged in theory and changes interior standard errors in proportion
to `1 / sqrt(nu)`. Use `terradish_rescale_nu()` to inspect that sensitivity
without refitting. Boundary inference needs separate care. Choose conductance
complexity using held-out prediction and report inference across plausible
`nu` values; see `vignette("model-comparison", package = "terradish")`.

## Cropping and coarse warm starts

Cropping removes routes outside the retained domain and can increase resistance
between sites even when every site remains inside the crop. The supplied
release audit found increases of 15–23% with a two-cell buffer and 1–2% with a
ten-cell buffer. Those are results for that landscape, not universal buffer
recommendations. Compare progressively wider buffers while retaining the same
response, sites, formula, and covariate scaling.

```{r crop-sensitivity, eval = FALSE}
# Scale on the original domain once, before cropping.
buffers <- c(2, 5, 10) * max(terra::res(covariates))
fits_buffer <- lapply(buffers, function(buffer) {
  cropped_surface <- conductance_surface(
    covariates, melip.coords, directions = 8, crop_buffer = buffer
  )
  terradish(melip.Fst ~ forestcover + altitude, cropped_surface,
            measurement_model = mlpe)
})
lapply(fits_buffer, coef)
lapply(fits_buffer, confint)
```

**How to read the output.** Compare estimates and intervals across buffers,
and inspect fitted resistance distances with `resistance_distances()`. Continue
expanding the domain if the substantive conclusion changes. The small-buffer
message is a diagnostic prompt, not proof that a larger buffer is adequate.

Coarse rasters can provide starting values, but the final fit must be refined
on the original graph. Coarsening changes both routes and covariate summaries.
The supported coarse-raster workflow therefore requires exact refinement:

```{r coarse-warm-start, eval = FALSE}
fit <- terradish(
  melip.Fst ~ forestcover + altitude, surface,
  measurement_model = mlpe,
  approximation = "coarse_raster",
  approximation_control = list(factor = c(4, 2), exact_refine = TRUE)
)
summary(fit)
```

Inspect convergence after the final refinement. A converged coarse fit does
not establish convergence on the original graph. `terradish_grid()` accepts
only `approximation = "none"`; coarse-only likelihoods are not supported
for ranking models.

---''')
s = s.replace('<!-- Phase 7: rewrite -->\n\n', '')
s = s.replace('  directional model, which uses a separate algorithm.\n', '')
s = s.replace('| `slim_terradish()` |', '| `crop_to_focal_buffer()` | Crop with an explicit map-unit buffer; check domain sensitivity |\n| `assess_terradish_settings()` | Compare computational settings on the current graph |\n| `slim_terradish()` |')
p.write_text(s, encoding='utf-8')
