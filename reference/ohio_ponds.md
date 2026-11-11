# Simulated Ohio forest-pond SNP and microsatellite data

Matched SNP and microsatellite genotypes for 20 retained diploid animals
in each of 24 hypothetical forest ponds on a real southeastern Ohio
landscape. The neutral SLiM simulation was motivated by a forest
pond-breeding salamander. Pond locations, hydroperiods, demographic
parameters, and genetic data are simulated. Elevation and land cover are
derived from public USGS products.

## Usage

``` r
data(ohio_ponds)
```

## Format

A list with 12 elements:

- `covariates`:

  A packed
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with 102 rows and 102 columns at 100 m resolution in NAD83 / UTM zone
  17N (EPSG:26917). Layers are `forest` (fraction from 0 to 1), `slope`
  (degrees), and `elevation` (meters).

- `forest500`:

  A packed raster on the same grid containing forest fraction smoothed
  with a Gaussian kernel of standard deviation 500 m. This is the
  context layer used to generate pond exchange.

- `ponds`:

  A 24-row data frame in population order. `pond` is the identifier; `x`
  and `y` are UTM coordinates in meters; `cell` is the 100 m raster
  cell; `slope` is in degrees; `edge_m` is distance to the domain edge
  in meters; `forest_buffer_100m`, `forest_buffer_300m`, and
  `forest_buffer_500m` are cell-center forest fractions in circular
  buffers of the indicated radius. `hydroperiod` is a hypothetical pond
  covariate in days with no simulated genetic effect; `nlcd_30m` is the
  native land-cover code at the pond; `displacement_from_trace25_m` and
  `candidate_count_at_selection` record the location-screening
  procedure; `focal` flags a sampled pond; `fold` gives three spatial
  validation blocks.

- `samples`:

  A 480-row data frame in genotype row order. `individual` and `pond`
  are identifiers; `population` is the one-based pond index;
  `pedigree_id`, `parent_1`, and `parent_2` are SLiM pedigree
  identifiers (`-1` denotes an unavailable parent); `x` and `y` are the
  pond coordinates.

- `raw_samples`:

  The same columns for 780 initial and additional captures before
  pedigree screening. Sampling began with 30 per pond, with additional
  captures where necessary to retain 20.

- `snp`:

  A 480 by 2500 integer matrix of alternate-allele dosage (0, 1, or 2),
  with individual row names and locus column names.

- `raw_snp`:

  A 780 by 2500 dosage matrix for the raw captures, using exactly the
  same loci as `snp`.

- `snp_positions`:

  The 2500 SNP positions in base pairs along a 50 Mb simulated
  chromosome, in column order.

- `microsat`:

  A 480 by 48 integer matrix of repeat counts at 24 stepwise-mutation
  loci. Adjacent columns hold the two alleles at one locus. Row names
  match `snp`; values are repeat counts, not dosage.

- `raw_microsat`:

  A 780 by 48 matrix for the raw captures at the same microsatellite
  loci.

- `microsat_loci`:

  A length-48 character vector assigning each microsatellite allele
  column to its locus.

- `simulation`:

  A named list recording simulation and mutation seeds, the 1200
  compressed generations, demographic and marker parameters, the known
  500 m forest context, the raw log-conductance coefficient of 6, the
  maximum pond-exchange probability of 0.025, the null hydroperiod
  effect, and the assumptions used in the teaching analyses.

## Source

USGS 3D Elevation Program and USGS Annual National Land Cover Database
2023 (Collection 1.0), processed to a common projected grid. Source
URLs, hashes, processing, and simulation provenance are recorded in
`system.file("extdata", "ohio-ponds", "PROVENANCE.md", package = "terradish")`.

## Details

Call
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html)
on the two packed raster elements before use. No downloads or simulation
software are needed to analyze the bundled data. All retained
individuals have both marker panels from the same ancestry. Screening
removes known immediate pedigree relationships within ponds; this is not
a claim that the retained animals are unrelated at deeper generations or
that such relationships would be known in field samples.

Ponds were screened for native forest, at least 90 percent forest within
100 m, and at least 60 percent within 300 m. These are hypothetical
design rules, not estimated habitat thresholds. Exchange has strong
reproductive philopatry with a long tail among ponds. The exchange model
uses endpoints on a resistance graph; it does not simulate observed
animal travel paths.

The SNP panel was filtered for pooled retained-sample minor allele
frequency of at least 0.05. Locus spacing does not establish independent
loci. The effective degrees of freedom used in a Wishart example is an
assumption, distinct from the SNP count. Spatial validation conditions
on this globally ascertained panel, rather than validating marker
discovery itself.

The source history was specified before fitting. Other histories and SNP
panel sizes were checked during development; the bundled history was not
selected for the strongest fit. The simulation illustrates an estimable
mechanism, not calibration or validation of a species-specific
demographic model. See the installed provenance file for sources and
processing details.

## Worked examples

[`vignette("v12-data-generation", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v12-data-generation.md),
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md).

## See also

[`cov_from_biallelic`](https://rdrr.io/pkg/landgraph/man/cov_from_biallelic.html),
[`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html),
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`vignette("v01-getting-started", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v01-getting-started.md),
[`vignette("v07-gaussian-scale-optimization", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v07-gaussian-scale-optimization.md).

Other datasets:
[`individual_salamanders`](https://wpeterman.github.io/terradish/reference/individual_salamanders.md),
[`melip`](https://wpeterman.github.io/terradish/reference/melip.md)

## Examples

``` r
data(ohio_ponds)
covariates <- terra::unwrap(ohio_ponds$covariates)
dim(ohio_ponds$snp)       # individuals by SNPs
#> [1]  480 2500
dim(ohio_ponds$microsat)  # individuals by two allele columns per locus
#> [1] 480  48
table(ohio_ponds$samples$pond)
#> 
#> P01 P02 P03 P04 P05 P06 P07 P08 P09 P10 P11 P12 P13 P14 P15 P16 P17 P18 P19 P20 
#>  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20  20 
#> P21 P22 P23 P24 
#>  20  20  20  20 
terra::res(covariates)   # map units are meters
#> [1] 100 100

# Aggregate diploid dosage to counts of alternate alleles in each pond.
groups <- factor(ohio_ponds$samples$pond, levels = ohio_ponds$ponds$pond)
Y <- rowsum(ohio_ponds$snp, groups, reorder = FALSE)
Y <- Y[ohio_ponds$ponds$pond, , drop = FALSE]
N <- 2 * as.numeric(table(groups))
S <- cov_from_biallelic(Y, N) / N[1]
round(S[1:4, 1:4], 3)
#>        P01   P02   P03    P04
#> P01  0.336 0.053 0.025 -0.015
#> P02  0.053 0.260 0.067  0.031
#> P03  0.025 0.067 0.156  0.031
#> P04 -0.015 0.031 0.031  0.293
```
