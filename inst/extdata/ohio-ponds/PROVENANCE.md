# Ohio forest-pond teaching data

The bundled `ohio_ponds` object is an analysis-ready extract from the
prespecified first history, `teaching_pond_v3_79101`. It contains simulated
genotypes and hypothetical ponds on a real landscape. It is not field data,
a map of known breeding ponds, or a fitted demographic model of spotted
salamanders. The parameters are teaching assumptions.

## Landscape inputs

The map uses NAD83 / UTM zone 17N (EPSG:26917). The 100 m analysis grid has
102 rows and 102 columns, with bounds x = 367650 to 377850 m and
y = 4365900 to 4376100 m. It retains all 10,404 cells, including matrix
habitat between ponds. A surrounding 30.6 km input domain was used during
landscape preparation. Elevation was sampled bilinearly at 30 m; land cover
was sampled by nearest neighbor. The continuous `forest` layer aggregates
forest indicator cells onto the 100 m grid. Forest includes NLCD deciduous,
evergreen, mixed forest, and woody wetlands (codes 41, 42, 43, and 90).
Elevation is in meters and slope is in degrees. Forest fractions are from
0 to 1. Forest within a circular screening buffer is a cell-center estimate.

The additional `landcover100.rds` file is a packed factor raster on the same
grid. It uses a modal resample of the native NLCD classes, then groups codes
41, 42, 43, and 90 as forest and all remaining classes as nonforest. This
dominant-class map is distinct from mean forest fraction. The broad
nonforest group is an educational simplification, not a biological claim
that water, roads, crops, and open habitat have the same conductance.

Public source services:

- USGS 3D Elevation Program:
  <https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer>
- USGS Annual NLCD land cover:
  <https://di-nlcd.img.arcgis.com/arcgis/rest/services/USA_NLCD_Annual_LandCover/ImageServer>
  The selected record is `Annual_NLCD_LndCov_2023_CU_C1V0`, year 2023,
  Collection 1.0, service object ID 39.

The original 30 m exports have SHA-256 hashes:

| Export | SHA-256 |
|---|---|
| DEM | `dc90b397dc9fe1a536196ca8e582d35544781258fd4ff9e27aa201c2ccc7bfe9` |
| Annual NLCD 2023 | `22c268caddea04806e5d87ec6d494eb5e175fc9baa54618280d13957b6dc1823` |

Derived run inputs:

| Input | SHA-256 |
|---|---|
| `covariates_100m.tif` | `5e4044b3a1d770b3902652ca03b49289ade2f8f2931895d0391ed64bf7b0593d` |
| `truth_100m.tif` | `ebc560782adf9c113e1421ffdf72cd688f1fa5e5400f1e9ca2a298992bfa5cc4` |
| Original pond table | `b25c1aed7f4371fc8a627071f50f57e45f49e31e6682a113a9963f5b0c2654c9` |

## Biological design and ancestry

All 24 pond points fall in native deciduous, evergreen, or mixed forest.
Screening required at least 90% forest within 100 m and at least 60% within
300 m, slope at most 12 degrees, spacing at least 600 m, and a domain-edge
margin of at least 1500 m. These are design rules, not habitat thresholds
estimated from field observations. Hydroperiod is a hypothetical covariate
with no effect on exchange in this history.

The forest layer is smoothed with Gaussian standard deviation 500 m. The
generating log-conductance is six times the centered, unstandardized
smoothed forest fraction. A symmetric eight-neighbor graph is traced to
the ponds by a Schur complement. Its exchange probabilities have maximum
off-diagonal row sum 0.025. Realized reproductive exchange is assessed from
surviving offspring. It is not a simulated route or an observed animal
movement track. Most reproduction remains at the natal pond, with a long
tail of exchange among more distant ponds.

SLiM 5.2 used a neutral nonWF model, target 120 diploids per pond, and 1200
compressed, nonoverlapping generations. One parent can originate from the
exchange distribution and the other is local. This simplified lifecycle is
not an annual age-structured salamander model. The simulated chromosome is
50 Mb with recombination rate 2e-8 per base per generation. Structured
recapitation uses population size 120 and half the forward exchange matrix
as a backward-time approximation. This is not proof of demographic
equilibrium or exact coalescent equivalence.

Seeds are 79101 for the history, 79102 for recapitation, 79103 for SNP
mutation, and 79104 for microsatellite mutation. Python processing used
tskit 1.0.3, pyslim 1.1.1, msprime 1.4.2, and NumPy 2.3.5.

## Matched marker panels and sampling

Sampling begins with 30 animals per pond in a genotype-independent capture
order. Additional captures are taken when needed to retain 20 without
known immediate pedigree relationships within that pond. This history has
780 raw captures and 480 retained animals. Older shared ancestry remains
possible. A pedigree screen in a simulation is not a replacement for
estimating relatedness and its uncertainty in field data.

The current 2500-SNP panel is drawn from the same ancestry and individuals
as the original demonstration, with JC69 mutation rate 1e-8. It retains
biallelic loci with pooled retained-sample minor allele frequency at least
0.05, uses candidate spacing of at least 10 kb, and deterministically
spreads the panel along the chromosome. Nested 1000- and 500-SNP panels
were evaluated during development but are not additional packaged objects.
The chromosome is linked; spacing does not establish independent loci.
Validation is conditional on a panel ascertained using all retained
animals, including those in held-out ponds.

The matched microsatellite panel uses 24 stepwise-mutation loci, mutation
rate 0.001 per locus per generation, ancestral repeat count 20, and repeat
bounds 5 to 45. Both allele copies are stored explicitly. SNP dosage and
microsatellite repeat counts are different input representations.

The source recapitated tree sequence has SHA-256
`a80d169c9e9f21e47d589c4f45a583a7811e6018711cf9b18afa2893362a14ff`.
The 2500-SNP retained dosage export has SHA-256
`4bc6607aeef55e47fe2c197faa96ecc34a94aa63cb06da7d2e78803a0313b540`.
The raw dosage export has SHA-256
`7ee9b99d11eaef6c11a22c164620fac1dd74245352b045fca096b6f25cdb2d75`.
The position export has SHA-256
`587d8e339b88d5ef1bec6eda1e783bbea577b257b714cb08dca399a64df6d2ff`.

## Analysis and reproducibility scope

The package contains raw and retained calls, pedigree metadata, landscape
values, and fixed spatial blocks. Reproduce the displayed analyses by
rendering copies of the installed `vignette-source` files. Downloads,
SLiM, Python, and tree-sequence processing are not required for those
analyses. Precomputed vignette outputs keep routine CRAN builds modest.

The full ancestry, simulation producers, original exports, and all three
independent histories are retained in the separate development archive.
They are not included in this small R dataset. Recreating the forward
simulation requires that archive and its pinned software environment;
the bundled analysis object alone cannot recreate the genealogy.

The first history was specified before fitting. Development included three
histories and 500-, 1000-, and 2500-SNP panels. The 500 m generating context
remains fixed across histories; estimates vary with history and panel.
These checks do not constitute a power study or calibrated interval coverage.

The teaching Wishart information value, `nu = 40`, is a sensitivity assumption.
It is neither the number of SNPs nor an estimate of independent genetic
features. Report it, examine sensitivity, and avoid comparisons of raw
likelihoods across different panels, marker types, or likelihood families.
