# Individual salamander teaching data

Use `data(individual_salamanders)` and
`vignette("individual-salamanders", package = "terradish")`.

The primary 2,500-SNP analysis is a positive nonlinear IBR demonstration.
The 24 matched microsatellites have mixed predictive evidence and remain an
explicit comparison. No independent IBE mechanism was manipulated.

## Installed files

- `drying_terms_10m.tif`: fixed `drying1` and `drying2` predictors on all
  7,056 cells. Use the saved center and divisors in the dataset metadata.
- `environment_10m.tif`: interpolated drying and forest cover. Interpolation
  from 30 m to 10 m supplies graph support, not new measured information.
- `terrain_30m.tif`: elevation, smoothed terrain wetness, northness, and eastness.
  All rasters use EPSG:26917 with meter coordinates. The domain is 840 by 840 m.
- `assessment.rds`: `original` and `independent_assay` comparisons, each with
  scenario (`arm`), marker panel, validity, uniform and IBR RMSE, proportional
  RMSE reduction, per-fold covariance-score gains, total gain, and number of
  positive folds. `block_contrasts` contains 2.5%, median, and 97.5% rows for
  conductance ratios C(0.65)/C(0.20) and C(0.65)/C(0.95). `independent_scores`
  contains the 40 frozen-model evaluations on the second marker assay.
- `predictions.rds`: 40 data frames, one per scenario, marker, model, and
  partition. Names have the form `compensatory_500__5000__SNP2500__IBR__1`.
  The 5000 records ancestral effective size; final 0 denotes full-data prediction,
  and 1 through 4 denote held-out spatial quadrants. Columns `id1`, `id2`,
  `observed`, and `predicted` identify pairs and squared allele-frequency distances.
- `independent-assay.rds`: two scenarios, each containing `snp` (160 by 2500)
  and `microsat` (160 by 48) matrices matched to the original ordered adults.
  These are new genomic realizations of the same fixed ecological pedigrees.
- `teaching-helpers.R`: local table, conductance-curve, and frozen-prediction
  utilities. They are educational code, not new exported package functions.
- `refit.R`: opt-in serial refits using installed data and existing public APIs.
  Copy to a writable working directory, set `TERRADISH_RUN_EXTENDED_EXAMPLE=true`,
  then source it. It checkpoints all starts into `refit-results`. Check shared
  resources before opting in. Stored results need no refits to read the guide.
- Three PNG figures are the checked illustrations used by the vignette.
- `provenance/`: source download receipts, roster protocol, original analysis
  code hashes, 200 block-resampling estimates, nested-panel sensitivity records,
  and the pre-generation freeze for the second assay.

## Model record dictionary

Records in `individual_salamanders$models` contain these fields:

| Fields | Meaning |
|:--|:--|
| `arm`, `Ne`, `panel`, `model`, `partition` | Scenario, ancestral effective size, marker panel, conductance model, and full-data (0) or held-out (1 through 4) partition |
| `valid`, `attempts`, `selected_start` | Accepted diagnostic status, all attempted starting-point records, and the selected start |
| `n_train`, `n_scored`, `nu` | Numbers of training and scored individuals, and the fixed objective multiplier (1) |
| `theta` | Named linear and quadratic coefficients on log conductance; empty for uniform conductance |
| `phi` | `tau`, the natural nonnegative graph-covariance scale, and `sigma`, log nugget variance |
| `logLik` | Stored objective value; not used here for nominal likelihood inference |
| `convergence` | Optimizer code, message, iterations, maximum absolute projected gradient, stopping criterion, boundary flag, and restart flag |
| `interior`, `within_diagnostic_guard` | Legacy bound-diagnostic placeholder (NULL for the unbounded model) and finite safety-guard check |
| `raw_gradient_tolerance_met` | Whether the ordinary, unprojected gradient meets the recorded tolerance |
| `start_logLik_spread`, `all_starts_same_optimum`, `all_starts_converged` | Objective spread and agreement across all starting points |
| `curvature_positive`, `graph_scale_positive`, `covariance_PD` | Positive objective curvature, positive graph contribution, and positive covariance on contrasts |
| `normalized_covariance_score`, `RMSE` | Held-out or full-data second-moment score and squared-distance root mean squared error; lower is better |
| `minimum_mean`, `maximum_mean`, `means_in_bounds` | Predicted squared-distance range and its check against 0 to 2; no clipping |
| `conductance_D020_D065_D095`, `hump` | Conductance at drying 0.20, 0.65, and 0.95, and the expected intermediate-drying peak check |

## Response and interpretation

The main response is the centered covariance of unstandardized diploid
allele-frequency profiles, averaged over loci. Construct it with
`cov_from_genetic_data()` as shown in the vignette, with explicit frequency
normalization. `dist_from_cov()` supplies its coherent squared distances.
`cov_from_biallelic()` and `dist_from_biallelic()` instead weight loci by pooled
allele frequency. Changing response requires refitting every candidate;
the stored fits cannot be attached to an alternative response.

Use `wishart_covariance` with `nu=1` as a covariance quasi-likelihood here.
Linked loci and diploid frequencies are not literal independent Gaussian
replicates. Do not use nominal AIC, Hessian standard errors, or significance
tests as genetic-sampling evidence. Block ranges condition on this ecological
history and sample. The spatial folds were reused during development; the
second assay uses frozen models but the same pedigrees and individuals.

Samples were selected from the saved ecological history without genetic data:
160 adults in 16 strata, at least 10 m apart, distinct analysis cells, with no
direct or shared recorded parent or grandparent. Distant relationships remain.
There is one ecological history per scenario, not ecological replication.
All earlier unsuccessful response/model comparisons remain in the development
archive; neither ecology nor the original genotypes was changed for this guide.

## Landscape attribution and biological motivation

The real landscape is near Logan in southeastern Ohio. Elevation derives from
the [USGS 3DEP service](https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer).
Forest derives from 2023 Annual NLCD Collection 1.0. Exact requests and source
hashes are in `provenance/landscape_download_manifest.json`. Drying is a
hypothetical terrain-derived index, with nonforest native cells assigned 1
before interpolation. No actual salamander occurrence or water-loss measurements
are provided. The USGS inputs are public source products; the simulated
genotypes are teaching data distributed with terradish under its package license.

Peterman et al. (2014), *Ecological resistance surfaces predict fine-scale
genetic differentiation in a terrestrial woodland salamander*, Molecular Ecology
23, doi:10.1111/mec.12747 motivates the compensatory-movement hypothesis.
It does not calibrate the numerical rates. Both sexes can live up to 20 years;
intermediate drying increases speed and heading persistence, while very dry
terrain adds movement danger. The control holds behavioral drying at 0.5 but
retains forest settlement and demography, so it is not a uniform-habitat control.

IBR describes connectivity through the intervening raster. MLPE endpoint
differences enter a pairwise mean; the Wishart environmental construction adds
centered environmental covariance kernels with nonnegative weights. These are
different statistical constructions. Neither alone demonstrates causal IBE.
