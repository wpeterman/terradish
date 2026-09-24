from pathlib import Path
import json
import re

p = Path('dev/check/phase7-package')
desc = p / 'DESCRIPTION'
s = desc.read_text(encoding='utf-8').replace('Version: 0.0.55', 'Version: 0.0.56')
a, b = s.index('Description:'), s.index('URL:')
description = '''Maximum likelihood estimation of landscape conductance surfaces
  from genetic data on raster graphs, using Newton and quasi-Newton optimization
  with analytic derivatives. Fits log-linear conductance models with optional
  joint estimation of Gaussian smoothing scales and spline response shapes,
  linked to pairwise genetic distances through maximum likelihood population
  effects models or to allele-frequency covariance through a Wishart likelihood
  on site contrasts, with the same pairwise distance and environment covariates
  available in both. Provides exact sparse and algebraic multigrid solvers for
  large rasters, fixed-domain spatial cross-validation for comparing conductance
  formulas, and simulation tools for study design. Works natively with 'terra'
  objects.'''
s = s[:a] + 'Description: ' + description + '\n' + s[b:]
desc.write_text(s, encoding='utf-8')
cff = p / 'CITATION.cff'
s = cff.read_text(encoding='utf-8').replace('version: 0.0.47', 'version: 0.0.56').replace('2026-09-22', '2026-09-24')
cff.write_text(s, encoding='utf-8')
zenodo = p / '.zenodo.json'
z = json.loads(zenodo.read_text(encoding='utf-8'))
z['description'] = ' '.join(description.split())
z['version'] = '0.0.56'
zenodo.write_text(json.dumps(z, indent=2) + '\n', encoding='utf-8')
news = p / 'NEWS.md'
s = news.read_text(encoding='utf-8')
news.write_text('''terradish 0.0.56
----------------
* Rebuilt the eight core guides and README around conditional interpretation,
  effective Wishart information, predictive comparison, and common pairwise inputs.
* Replaced unavailable precomputed power sections with explicit reproducible
  study-design recipes, and removed unverified timing claims.
* Updated package and citation metadata for the supported core scope.

Planned terradish 0.1.0 core release
----------------------------------
* Scope: directed, hierarchical, drift, selected-pair, Kron, block-CG, and legacy
  cross-validation prototypes remain on the experimental branch because their
  identifiability, statistical assumptions, or numerical behavior need further work.
* Breaking: covariance Wishart now uses site contrasts. Environmental Wishart
  kernels use the same pairwise transforms as MLPE, defaulting to absdiff;
  normalize is removed and shared pairwise objects are accepted.
* Breaking: graph neighborhoods default to eight directions. Supported solvers
  are direct, auto, and AMG. Approximate starts require exact refinement;
  terradish_grid accepts only approximation = "none". Legacy CV is removed.
* Breaking: Gaussian outer searches default to the terradish kernel, count
  estimated scales in degrees of freedom, and use map units consistently with
  simulation. Joint Gaussian upper bounds reflect truncated kernel support.
* Breaking: spline prediction retains fitted knots and column centering, which
  can shift absolute conductance by a constant without changing relative shape.
* Fixes: even-grid Gaussian alignment, layerwise scaling, covariance centering,
  retained new-raster transformations, nested tests, boundary inference,
  optimizer stopping, and failed-fold accounting.
* New: covariance and interval methods, explicit convergence records and boundary
  restarts, Wishart nu rescaling, Gaussian scale profiles and bound summaries,
  focal spline monotonicity, stored raster scaling, and IBE:IBR ratio uncertainty.
* New: per-formula conductance and measurement models, repeated fixed-domain
  folds, fixed-nuisance prediction, common-fold paired comparisons, and safe
  checkpoint resumption. Power studies can vary nu_fit separately from nu.
* Documentation: revised effective-information guidance, educational core
  workflows, and limitations of environmental, smoothing, and conductance effects.

''' + s, encoding='utf-8')
(p / 'cran-comments.md').write_text('''## Submission status

This is preparation for the first CRAN release of terradish 0.1.0.
The current development version is 0.0.56. No submission has been made for
this core overhaul; owner release approval remains pending.

## Validation

Windows, R 4.6.1: the .54 robustness candidate passed a complete check with
0 errors, 0 warnings, and 0 notes. The .55 function-documentation candidate
passed examples and vignette rebuilding with the same clean result; its code
was unchanged from .54 and its redundant test run was skipped.
Final .56 source-archive and as-cran checks are pending. Earlier R-hub and
win-builder runs predate this overhaul and do not validate this candidate.

## Dependencies

terradish requires landgraph >= 0.0.3 for the documented covariance metadata
and grouped-data Gower default. This companion version has passed local checks
but its CRAN publication is still a release prerequisite. The private validation
library contains that version. Do not submit terradish before the required
landgraph release is available on CRAN.

## Scope

The core retains symmetric resistance likelihoods, supported conductance
factories, direct/auto/AMG solvers, and fixed-domain cross-validation.
Research prototypes are retained on the experimental branch. Final check
results and submission metadata will be updated after release approval.
''', encoding='utf-8')
# Replace the US-only projection in Brazilian example recipes with a local
# azimuthal equidistant projection. No distance is measured in raw degrees.
projection = '''# Center a local distance projection on the sampled Brazilian sites.
lonlat <- terra::crds(terra::project(melip.coords, "EPSG:4326"))
local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
                     mean(lonlat[, 2]), mean(lonlat[, 1]))
coords_projected <- terra::project(melip.coords, local_crs)'''
for name in ('model-comparison', 'spline-conductance', 'gaussian-scale-optimization'):
    f = p / 'vignettes' / f'{name}.Rmd'
    s = f.read_text(encoding='utf-8')
    s = s.replace('coords_projected <- terra::project(melip.coords, "EPSG:5070")', projection)
    s = s.replace('coords_projected <- terra::project(coords, "EPSG:5070")', projection.replace('melip.coords', 'coords'))
    f.write_text(s, encoding='utf-8')
for name in ('getting-started', 'gaussian-scale-optimization'):
    f = p / 'vignettes' / f'{name}.Rmd'
    s = f.read_text(encoding='utf-8')
    s = s.replace('```{r support-clamp, eval = FALSE}', '```{r support-clamp}')
    s = s.replace('```{r gaussian-support-clamp, eval = FALSE}', '```{r gaussian-support-clamp}')
    f.write_text(s, encoding='utf-8')
for f in (p / 'inst/examples').glob('*.R'):
    s = f.read_text(encoding='utf-8')
    s = re.sub(r'maxit = (8|12),', 'maxit = 100,', s)
    if f.name == 'covariance-response-power-example.R':
        s = s.replace('seed = 2026)', 'seed = 2026, nu_fit = nu)')
        s = s.replace('    nu = nu,', '    nu = nu,\n    nu_fit = nu_fit,')
        s = s.replace('"Use selected_AICc_rate to ask whether the design can distinguish the flexible candidate from simpler alternatives."', '"Treat selected_AICc_rate as model-conditional recovery; check held-out prediction and sensitivity to nu_fit before selecting flexibility."')
    f.write_text(s, encoding='utf-8')
