"""Install the reviewed vignette without changing its scientific results."""
from pathlib import Path
import shutil
import re

root = Path(__file__).resolve().parent.parent
bundle = root / "dev/amphibian_example/terrestrial_individual_candidate_20261023/vignette-candidate"
ext = root / "inst/extdata/individual-salamanders"
ext.mkdir(parents=True, exist_ok=True)
text = (bundle / "vignettes/individual-salamanders.Rmd").read_text(encoding="utf-8")
start = text.index("<!-- Development vignette")
end = text.index("# What you will learn")
setup = '''```{r setup, include=FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
  fig.width = 6, fig.height = 4.25)
```

'''
text = text[:start] + setup + text[end:]
text = text.replace('    self_contained: true',
    '    self_contained: true\n    mathjax: null\n    pandoc_args: ["--mathml"]')
# Native MathML supports offline reading; give equations visible labels because
# Pandoc's MathML writer does not display TeX's tag command.
text = re.sub(r'\$\$\n(.*?)\\tag\{(\d+)\}\n\$\$',
    lambda match: f'**Equation {match[2]}.**\n\n$$\n{match[1]}$$', text, flags=re.S)
text = text.replace('The animals are simulated. The terrain and land-cover inputs come from southeastern Ohio, near Logan, but neither actual salamander occurrences nor measured water loss are provided. The life cycle is inspired by terrestrial *Plethodon* salamanders.',
    'These are simulated animals on real terrain near Logan in southeastern Ohio. No actual salamander occurrences or measured water loss are provided. Terrestrial *Plethodon* salamanders inspire the life cycle.')
text = text.replace('2 x 10^-8', '$2 \\times 10^{-8}$')
text = text.replace('10^-8', '$10^{-8}$').replace('10^-4', '$10^{-4}$')
anchor = "# The landscape, animals, and sampling"
load = '''# Load the matched example

No simulation software or download is needed. The fitted models, genomic-block
resampling, and predictive comparisons are bundled as checked results. You can
rebuild the genetic response and inspect predictions quickly; the longer fits
are shown separately and never run during a routine vignette build.

```{r load-example}
library(terradish)
data(individual_salamanders)
files <- system.file("extdata", "individual-salamanders",
                     package = "terradish", mustWork = TRUE)
source(file.path(files, "teaching-helpers.R"), local = TRUE)
example <- individual_salamanders
results <- readRDS(file.path(files, "assessment.rds"))
predictions <- readRDS(file.path(files, "predictions.rds"))
main <- example$data$compensatory_500
m <- example$metadata
record <- example$models$compensatory_500$SNP2500__IBR__0
```

The `data` element contains two scenarios with matched genotypes and coordinates.
The `models` element contains coefficient and diagnostic records, rather than
live fitted objects. The `metadata` element fixes the environmental transformations.
Read `?individual_salamanders` for the data dictionary. The small helper script
formats comparison tables and evaluates conductance curves; genetic responses
below use the package's public functions.

'''
text = text.replace(anchor, load + anchor, 1)
for name in ("landscape-and-samples", "drying-response-and-markers", "predictive-comparisons"):
    text = text.replace(f'file.path(bundle,"figures","{name}.png")',
        f'file.path(files, "{name}.png")')
for name in ("landscape-and-samples", "drying-response-and-markers", "predictive-comparisons"):
    text = text.replace(f'knitr::include_graphics(file.path(files, "{name}.png"))',
        f'knitr::include_graphics(file.path(files, "{name}.png"), rel_path = FALSE)')
old = '''# These helpers make the construction explicit and work without extra packages.
C_snp <- snp_covariance(G)
C_msat <- microsat_covariance(A)'''
new = '''# Keep loci unstandardized: this is the response validated for this example.
# Each biallelic locus contributes two allele-frequency features, hence 1/2.
C_snp <- cov_from_genetic_data(G, scale = FALSE, diagonal = "gower",
                                normalize = "none") / (2 * ncol(G))

# Allele-call features count two gene copies. Divide by 4 to get frequency
# inner products, then average across the 24 loci, including monomorphic loci.
L_msat <- ncol(A) / 2
C_msat <- cov_from_genetic_data(A, input = "allele_calls",
  loci = rep(seq_len(L_msat), each = 2), scale = FALSE,
  diagonal = "gower", normalize = "none") / (4 * L_msat)'''
assert old in text
text = text.replace(old, new)
start = text.index("The public `landgraph` helper gives the same result")
end = text.index("The squared genetic distance implied by a covariance", start)
text = text[:start] + text[end:]
text = text.replace("D_snp <- squared_genetic_distance(C_snp)",
                    "D_snp <- as.matrix(dist_from_cov(C_snp))")
anchor = "# Build a nonlinear drying model"
alternative = '''## Why not use the biallelic helpers here?

`cov_from_biallelic()` is also a valid package entry point. For an individual
dosage matrix, give each row two sampled gene copies. Its locus-specific
denominator uses the pooled allele frequency *p*, so rare and common alleles
receive different weights. Equation 1 deliberately leaves loci unstandardized,
which also gives a common allele-frequency interpretation for matched
microsatellites. A constant rescaling cannot remove locus-specific weights.

```{r biallelic-alternative}
N <- rep(2, nrow(G))
C_weighted <- cov_from_biallelic(G, N)
D_weighted <- as.matrix(dist_from_biallelic(G, N))
stopifnot(max(abs(D_weighted - as.matrix(dist_from_cov(C_weighted)))) < 1e-10)
```

**How to read this construction:** The two biallelic functions give coherent
covariance and squared-distance representations of the same weighted response.
That response has a different scale and need not lie between 0 and 2. It is
not the response used for the fitted results below. If you choose it, refit
every candidate and repeat predictive checks using that same weighted response;
do not attach the stored fits to a newly defined genetic quantity.

'''
text = text.replace(anchor, alternative + anchor, 1)
text = text.replace('terms <- rast(file.path(bundle,"data","drying_terms_10m.tif"))',
                    'terms <- rast(file.path(files, "drying_terms_10m.tif"))')
text = text.replace('The fixed values *m*, *s*1, and *s*2',
                    'The fixed values $m$, $s_1$, and $s_2$')
text = text.replace('quadratic$fit$phi["tau",]',
                    'quadratic$fit$phi["tau", "Estimate"]')
text = text.replace('quadratic$fit$phi["sigma",]',
                    'quadratic$fit$phi["sigma", "Estimate"]')
text = text.replace('contrast_ranges <- assay$contrast_central_95_ranges',
                    'contrast_ranges <- results$block_contrasts')
old = '''uniform_test <- read.csv(file.path(bundle,"results",
  paste0(prefix,"uniform__1__predictions.csv")))
drying_test <- read.csv(file.path(bundle,"results",
  paste0(prefix,"IBR__1__predictions.csv")))'''
new = '''uniform_test <- predictions[[paste0(prefix, "uniform__1")]]
drying_test <- predictions[[paste0(prefix, "IBR__1")]]'''
assert old in text
text = text.replace(old, new)
text = text.replace('comparison_table(assessment$comparisons)', 'comparison_table(results$original)')
text = text.replace('comparison_table(validation$comparisons)', 'comparison_table(results$independent_assay)')
text = text.replace('Use the same graph and sampling locations for both competing models. Adjacent cells connect in eight directions. The uniform model gives all cells equal conductance. The alternative makes log conductance quadratic in drying:',
'''Use the same graph and sampling locations for both competing models. Adjacent
cells connect in eight directions. The uniform model gives all cells equal
conductance, while retaining geographic positions, raster connections, and
domain boundaries. The comparison asks whether drying adds predictive value
beyond that geometry. The alternative makes log conductance quadratic in drying:''')
anchor = "Inspect one fold's predictions before reading the pooled comparison:"
cv_workflow = '''For a new analysis, `terradish_cv_folds()` provides a fixed-landscape
workflow. Set `nuisance = "fixed"` so the held-out responses do not refit graph
scale or nugget variance. The default, `"reprofile"`, answers a different
question by reestimating measurement parameters on the test responses.

```{r optional-fixed-cv, eval=FALSE}
# Long-running: run explicitly, after constructing S and graph above.
cv <- terradish_cv_folds(
  graph, formulas = list(quadratic = S ~ drying1 + drying2),
  folds = samples$fold, model = wishart_covariance, nu = 1,
  nuisance = "fixed", baseline = TRUE,
  conductance_model = loglinear_conductance,
  optimizer = "bfgs", cores = 1, keep_fits = "slim",
  leverage = FALSE,
  control = NewtonRaphsonControl(maxit = 100, ctol = 1e-6, ftol = 1e-8),
  measurement_control = NewtonRaphsonControl(
    maxit = 100, ctol = 1e-7, ftol = 1e-9)
)
print(cv)
```

**How to read this comparison:** Positive `loglik_gain` favors quadratic drying
over the uniform baseline for that fold. In this example, interpret it as a
quasi-likelihood prediction score, not literal genotype log probability. A
failed fold is missing evidence; do not omit it from a favorable comparison.
The bundled results additionally audit three starting points, gradients, and
curvature. The installed refit script retains all those starting points.

'''
assert anchor in text
text = text.replace(anchor, cv_workflow + anchor, 1)
text = text.replace('The microsatellites are matched to the SNP samples, but 24 loci need not yield the same inference as 2,500 SNPs. Keep a mixed result visible. It teaches why the strength and consistency of predictive evidence matter more than a visually appealing response curve.',
'''The primary SNP model reduces held-out RMSE by 10.4%, with positive covariance-score
gains in all four folds. The control SNP model increases RMSE by 9.5%, so the
additional drying terms are not automatically useful in both scenarios.

Matched microsatellites give a different lesson. Their covariance score improves
overall, but squared-distance RMSE worsens by 3.4%, and one fold has a negative
score gain. The two measures assess different aspects of the predictions. These
24 loci therefore give mixed evidence, even though their fitted curve resembles
the SNP curve. Matching individuals removes a sampling difference; it does not
give 24 microsatellites the information in 2,500 SNPs. Do not compare raw RMSE or
objective values between marker types.''')
text = text.replace('**How to read the table:** The new assay asks whether the frozen connectivity model predicts another genomic realization conditional on the same ecological history. It does not provide independent ecological replication. Both successful and unsuccessful marker comparisons remain in the table.',
'''**How to read the table:** The new SNP assay reduces RMSE by 10.8%, with positive
covariance-score gains in three of four folds. The matched microsatellite RMSE
worsens by 3.0% despite an overall covariance-score improvement. The SNP result
therefore reproduces across genomic realizations of these same adults; the
microsatellite evidence remains mixed. Neither result supplies independent
ecological replication.''')
text = text.replace('S <- snp_covariance(G)',
'''S <- cov_from_genetic_data(G, scale = FALSE, diagonal = "gower",
                           normalize = "none") / (2 * ncol(G))''')
text = text.replace('| `cov_from_genetic_data()` | Construct genetic covariance with explicit scaling |',
'''| `cov_from_genetic_data()` | Construct the unweighted individual response with explicit normalization |
| `cov_from_biallelic()` / `dist_from_biallelic()` | Construct the alternative allele-frequency-weighted response |
| `dist_from_cov()` | Convert the chosen covariance to its coherent squared-distance representation |''')
text = text.replace('All previous failed response/model comparisons are retained in the development audit. The final workflow was revised with those results in view, so its exploratory history is part of the reproducible record.',
'''Earlier response and model comparisons failed to provide a clear teaching
example. They remain in the development archive. The final workflow was revised
with those results in view, so the reused spatial folds describe exploratory
prediction, not a prospectively untouched test.

For a longer reproduction, copy the installed `refit.R` and
`teaching-helpers.R` files to a writable directory. The refit script requires an
explicit opt-in, runs serially, saves every starting point, and can resume from
completed records. Ordinary rendering reads the stored results; it does not
silently rerun fitting, cross-validation, or simulation.''')
(root / "inst/vignette-source/individual-salamanders.Rmd").write_text(text, encoding="utf-8")
# This guide executes only response construction and stored-result inspection.
# Keep it executable in the routine vignette, just like wishart-covariance.
(root / "vignettes/individual-salamanders.Rmd").write_text(text, encoding="utf-8")

helpers = '''# Local teaching utilities, not additions to the public package API.
comparison_table <- function(comparisons) {
  out <- lapply(comparisons, function(x) data.frame(
    Comparison = paste(if (x$arm == "compensatory_500") "Compensatory" else "Control",
      if (x$panel == "SNP2500") "SNPs" else "microsatellites"),
    "Error reduction (%)" = 100 * x$RMSE_reduction,
    "Covariance score gain" = x$total_score_gain,
    "Positive folds" = paste(x$positive_folds, "of 4"), check.names = FALSE))
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

curve_from_record <- function(record, D, metadata) {
  theta <- unlist(record$theta)
  exp(theta["drying1"] * (D - metadata$center) / metadata$sd_linear +
      theta["drying2"] * (D - metadata$center)^2 / metadata$sd_quadratic)
}

predict_from_record <- function(record, graph) {
  theta <- unlist(record$theta)
  # Both zero coefficients give uniform conductance on the same graph.
  if (record$model == "uniform") theta <- c(drying1 = 0, drying2 = 0)
  columns <- names(theta)
  E <- terradish_distance(
    matrix(unname(theta), 1, dimnames = list(NULL, columns)),
    reformulate(columns), graph,
    conductance_model = loglinear_conductance, covariance = TRUE
  )$covariance[, , 1]
  phi <- unlist(record$phi)
  phi["tau"] * E + exp(phi["sigma"]) * diag(length(graph$demes))
}
'''
(ext / "teaching-helpers.R").write_text(helpers, encoding="utf-8")
for name in ("landscape-and-samples", "drying-response-and-markers", "predictive-comparisons"):
    shutil.copy2(bundle / f"figures/{name}.png", ext / f"{name}.png")

# Add navigation at stable anchors, preserving all existing development edits.
link = '| How do individual genotypes reveal nonlinear resistance? | `vignette("individual-salamanders", package = "terradish")` | Terrestrial SNP IBR with matched microsatellite comparison |\n'
for name in ("vignettes/getting-started.Rmd", "inst/vignette-source/getting-started.Rmd"):
    path = root / name
    content = path.read_text(encoding="utf-8")
    if link not in content:
        anchor = '| Which genetic response and information assumption should I use?'
        assert anchor in content
        content = content.replace(anchor, link + anchor, 1)
        path.write_text(content, encoding="utf-8")
for name in ("vignettes/ibe-ibr-workflow.Rmd", "inst/vignette-source/ibe-ibr-workflow.Rmd"):
    path = root / name
    content = path.read_text(encoding="utf-8")
    addition = '''For individual sampling, see `vignette("individual-salamanders", package = "terradish")`.
That SLiM example demonstrates nonlinear IBR, retains mixed microsatellite
prediction results, and explains endpoint IBE terms without claiming a separate
simulated IBE effect.

'''
    anchor = '## Three environmental hypotheses that should stay distinct'
    if addition not in content:
        assert anchor in content
        path.write_text(content.replace(anchor, addition + anchor, 1), encoding="utf-8")
path = root / "vignettes/wishart-covariance.Rmd"
content = path.read_text(encoding="utf-8")
addition = '''The worked individual-level example in
`vignette("individual-salamanders", package = "terradish")` uses an unweighted
allele-frequency response for matched SNPs and microsatellites. It contrasts
that response with the weighting in `cov_from_biallelic()` and checks prediction
without treating the SNP count as independent Wishart information.

'''
anchor = '## Preserve sampling order and marker meaning'
if addition not in content:
    assert anchor in content
    path.write_text(content.replace(anchor, addition + anchor, 1), encoding="utf-8")
path = root / "dev/refresh-precomputed-vignettes.R"
content = path.read_text(encoding="utf-8")
if '  "individual-salamanders",' in content:
    path.write_text(content.replace('  "individual-salamanders",\n', '', 1), encoding="utf-8")
path = root / "NEWS.md"
content = path.read_text(encoding="utf-8")
bullet = '''* Added a focused individual-genotype SNP resistance vignette with raw matched
  microsatellites, nonlinear drying, spatial prediction, and explicit IBE limits.
'''
if bullet not in content:
    content = content.replace('----------------\n', '----------------\n' + bullet, 1)
    path.write_text(content, encoding="utf-8")
print("Integrated vignette source, helpers, figures, and navigation.")
