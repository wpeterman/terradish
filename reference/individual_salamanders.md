# Simulated individual salamander genotypes and resistance analyses

Matched SNP and microsatellite genotypes for 160 adult terrestrial
salamanders in each of two simulated movement scenarios on a real
southeastern Ohio landscape. The primary example demonstrates nonlinear
isolation by resistance with individual sampling. The animals and drying
index are hypothetical; neither species occurrence nor field-calibrated
movement is claimed.

## Usage

``` r
data(individual_salamanders)
```

## Format

A list with three elements:

- `data`:

  Two named scenarios, `compensatory_500` and `control_500`. Each
  contains `samples`, `snp_dosage`, `microsat_alleles`, and
  `covariance`, described below.

- `models`:

  Twenty checked model records per scenario: uniform and quadratic
  drying models, two marker types, and partitions 0 (all adults) through
  4 (one spatial quadrant withheld). Names combine marker, model, and
  partition, for example `SNP2500__IBR__0`. Records contain
  log-conductance coefficients `theta`, graph covariance multiplier
  `phi$tau` on its natural nonnegative scale, log nugget variance
  `phi$sigma`, starting-point diagnostics, convergence, curvature,
  covariance score, RMSE, and predicted-mean range. These are plain
  lists, not fitted `terradish` objects.

- `metadata`:

  Fixed drying center `center`, linear divisor `sd_linear`, and
  quadratic divisor `sd_quadratic`; the latter is SD of squared
  deviations, not squared SD. `start_halfwidth` defines additional
  starting points, not coefficient bounds. `epsg`,
  `native_resolution_m`, `graph_resolution_m`, `directions`, `extent`,
  `core`, and `graph_nodes` describe the spatial grid (EPSG:26917,
  meters). `sample_size`, `main_SNPs`, and `matched_microsatellites`
  give panel sizes. `objective`, `biological_scope`, and
  `raw_marker_selection` state interpretation and selection rules.

## Source

Simulated genotypes from recorded individual pedigrees on terrain
derived from USGS 3DEP elevation and 2023 Annual NLCD Collection 1.0.
See
`system.file("extdata", "individual-salamanders", "README.md", package = "terradish")`
for data sources, processing, and limitations.

## Details

Within each scenario, `samples` has 160 rows in genotype order.
`sample_id` and `sample_key` identify the animal; `id`, `mom`, and `dad`
are recorded ecological pedigree identifiers. `stratum` gives one of 16
sampling strata; `quadrant` and `fold` identify the four spatial
validation blocks. `age` is in years, `sex` is 0 for female and 1 for
male, `x_local` and `y_local` are meters from the domain origin, and
`easting` and `northing` are projected UTM coordinates. `drying` is the
interpolated hypothetical index from 0 to 1.

`snp_dosage` contains nested integer matrices named `500`, `1000`, and
`2500`. Entries 0, 1, and 2 count one allele. Locus names identify
synthetic genomic blocks and within-block positions. `microsat_alleles`
is a 160 by 48 integer matrix; adjacent columns hold the two
repeat-count alleles at each of 24 loci. All marker row names match
`samples$sample_id`. Microsatellite allele labels are categorical, not
continuous response values.

`covariance` contains centered matrices named `SNP2500` and
`microsat24`, averaging unstandardized allele-frequency inner products
across loci. Their implied squared distances range from 0 to 2. This is
different from the locus-frequency weighting in
[`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html)
and is not a Dps response.

Each scenario follows 500 years, with maximum age 20 for both sexes. The
compensatory scenario allows faster, straighter movement at intermediate
drying but greater danger in very dry terrain. The control holds
behavioral drying at 0.5 while retaining forest settlement and
demography. It is a control for drying-dependent movement, not every
possible resistance effect. Samples occupy distinct 10 m cells, are at
least 10 m apart, and exclude direct or shared recorded parents and
grandparents. Distant kin can remain.

These analyses use `wishart_covariance` with `nu = 1` as a covariance
quasi-likelihood, not a literal independent-Gaussian-locus sampling
model. Nominal Hessian uncertainty, AIC, and significance tests are not
used as genetic-sampling evidence. Two hundred genomic-block resamples
and a second genomic assay check conditional marker variation; neither
supplies independent ecological replication. Spatial folds were reused
during workflow development. Marker discovery conditions on the full
sampled panel. The matched microsatellites retain mixed predictive
evidence. No separate causal isolation-by-environment effect was
simulated.

Installed files in `extdata/individual-salamanders` provide the 10 m
drying terms and environment, native 30 m terrain, all 40 original
prediction records, comparative scores, the independent assay, and
provenance. The 10 m grid interpolates native 30 m information; it adds
no new measurements. The vignette explains these files and the longer
opt-in reproduction.

## Worked examples

[`vignette("v12-data-generation", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v12-data-generation.md),
[`vignette("v09-individual-salamanders", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v09-individual-salamanders.md).

## See also

[`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html),
[`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html),
[`dist_from_cov`](https://rdrr.io/pkg/landgraph/man/dist_from_cov.html),
[`vignette("v09-individual-salamanders", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v09-individual-salamanders.md).

Other datasets:
[`melip`](https://wpeterman.github.io/terradish/reference/melip.md),
[`ohio_ponds`](https://wpeterman.github.io/terradish/reference/ohio_ponds.md)

## Examples

``` r
data(individual_salamanders)
animals <- individual_salamanders$data$compensatory_500
G <- animals$snp_dosage[["2500"]]
S <- cov_from_genetic_data(G, scale = FALSE, normalize = "none") /
  (2 * ncol(G))
dim(S)                              # 160 individual sampling locations
#> [1] 160 160
range(as.matrix(dist_from_cov(S)))   # squared allele-frequency distance
#> [1] 0.000 0.582
table(animals$samples$fold)           # 40 withheld adults per spatial fold
#> 
#>  1  2  3  4 
#> 40 40 40 40 
```
