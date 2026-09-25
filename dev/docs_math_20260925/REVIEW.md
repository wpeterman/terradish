# Documentation math audit, September 25, 2026

The 0.0.58 documentation correction repairs source encoding damage and
mathematical notation. Numerical implementation and public interfaces are
unchanged. Existing paragraph-wrapping edits in model-comparison were preserved.

## Corrections

- Model comparison: replaced corrupted multiplication, subtraction,
  chi-squared, and theta text with proper Markdown math or ASCII labels.
  Defined the likelihood-ratio statistic and its degrees of freedom, stated
  the interior-parameter condition for the ordinary chi-squared reference,
  and matched the table explanation to K, loglik, and Delta_AIC.
- Splines: repaired the minus sign, arrows, inequality, delta-AICc, and grid
  dimensions; corrected the introductory basis-function sentence.
- Getting started and Gaussian scales: standardized theta notation, coefficient
  subscripts, and exponential transformations. The Gaussian guide now uses
  theta consistently for conductance coefficients, reserving beta for the
  measurement-model slope elsewhere in the documentation.
- Function help: kept FST inside a single math expression; typeset the
  Gaussian quantile and radial-distance formulas; set descriptive subscripts
  upright; corrected the MLPE row-vector product from transposed Z to Z gamma;
  identified the omitted parameter-constant terms in the Wishart objective.
  All Rd edits were regenerated from roxygen sources.

## Accuracy checks

Equations were checked against the implementation: MLPE's additive mean and
shared-site correlation, the contrast-Wishart objective and nugget scale,
environmental Gower kernels, squared-distance centering, Gaussian kernel
quantiles/radii, conductance transformations, and information-criterion counts.
The simulation helpers use nugget variance on its natural scale; fitted
Wishart models store its logarithm. That documented distinction was retained.
No citations or numerical claims were added from unverified external sources.

## Validation

- Scanned 179 active documentation/code text files, including README, package
  help sources, and all eight vignettes: no remaining suspicious encoding
  sequences. Historical dev logs and frozen release candidates were excluded.
- Parsed and checked all 62 Rd files successfully; rendered every help page.
- Rebuilt all eight vignettes with their enabled code chunks evaluated.
  Model comparison and Gaussian scales were rebuilt once more after the final
  table-label and coefficient-name refinements. Final rendered inputs match
  the working vignette sources.
- Browser checks: all 61 vignette math expressions and 124 help-page math
  expressions rendered, with no MathJax/KaTeX error nodes. The six vignettes
  containing typeset math were checked for equation overflow; none was found.
  The other two guides contain code-formatted notation, which was inspected
  for corrupted symbols and unparsed math markup.
- Visual checks covered the likelihood-ratio equation and chi-squared
  reference, corrected model-selection table labels, theta subscripts and
  exponentials, Gaussian smoothing formula, Wishart generative equation,
  centering matrix, spline formulas and delta-AICc, and representative help
  equations. This is equation-focused inspection, not a new audit of every
  plot or every page's complete layout.
- All executable R syntax matches the checked 0.0.57 implementation. The
  only executable vignette change is a display label in model comparison.
  No numerical tests or full R CMD check were repeated for this documentation
  correction. The prior .57 archive checks are not presented as a .58 check.

R 4.6.1, RStudio's bundled Pandoc, and the private landgraph 0.0.3 library were
used. Browser inspection used the normal MathJax/KaTeX assets referenced by
the generated HTML. It does not establish offline rendering without those
assets. Source encodings remain UTF-8; math commands themselves use ASCII.

## Evidence and preview

Source inventory: source_inventory.json. R checks: validation.log and
check_rd.log. Full render receipts: render_results.rds; final refinements:
refined.log and refined_results.rds. Rendered hashes:
rendered_fingerprints.json. Scripts in this directory reproduce the checks.

Rendered pages are in dev/check/docs-math-20260925-render, including
model-comparison.html and the help subdirectory. These are local previews;
the old .57 release archive remains preserved unchanged. Rebuild the release
archive from current sources before distribution. No push or master merge
is included in this correction.
