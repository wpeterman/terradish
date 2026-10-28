# Assemble the prespecified first pond history; never select a history by fit.
# Run from the development package root. Ancestry and original exports stay intact.
library(terra)
root <- file.path('dev', 'amphibian_example')
run <- file.path(root, 'runs', 'teaching_pond_v3_79101')
panel <- file.path(root, 'snp_panel_sensitivity_20261012', 'datasets',
                   'teaching_pond_v3_79101', '2500')
read_calls <- function(path) {
  x <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(x) <- 'integer'
  x
}
ponds <- read.csv(file.path(run, 'inputs', 'ponds.csv'))
samples <- read.csv(file.path(run, 'samples.csv'))
raw_samples <- read.csv(file.path(run, 'raw_samples.csv'))
folds <- read.csv(file.path(root, 'snp_panel_sensitivity_20261012', 'analysis',
                            'teaching_pond_v3_79101', 'cv', 'spatial_folds.csv'))
stopifnot(identical(ponds$x, folds$x), identical(ponds$y, folds$y))
ponds$fold <- as.integer(folds$fold)
ohio_ponds <- list(
  covariates = wrap(rast(file.path(run, 'inputs', 'covariates_100m.tif')), proxy = FALSE),
  forest500 = wrap(rast(file.path(run, 'inputs', 'truth_100m.tif'))[['forest500']],
                   proxy = FALSE),
  ponds = ponds,
  samples = samples,
  raw_samples = raw_samples,
  snp = read_calls(file.path(panel, 'snp_dosage.csv')),
  raw_snp = read_calls(file.path(panel, 'raw_snp_dosage.csv')),
  snp_positions = scan(file.path(panel, 'snp_positions.csv'), quiet = TRUE),
  microsat = read_calls(file.path(run, 'microsat_alleles.csv')),
  raw_microsat = read_calls(file.path(run, 'raw_microsat_alleles.csv')),
  microsat_loci = rep(sprintf('MS%02d', 1:24), each = 2),
  simulation = list(
    history_seed = 79101L, recapitation_seed = 79102L,
    snp_seed = 79103L, microsat_seed = 79104L,
    generations = 1200L, target_diploids_per_pond = 120L,
    genome_length_bp = 50000000L, recombination_rate = 2e-8,
    snp_mutation_rate = 1e-8, microsat_mutation_rate = 0.001,
    microsat_repeat_bounds = c(5L, 45L), forest_sigma_m = 500,
    forest_log_conductance_coefficient = 6,
    maximum_exchange_probability = 0.025,
    hydroperiod_effect = 0, degrees_of_freedom_demo = 40,
    scope = paste('Neutral teaching simulation on real terrain and land cover;',
                  'hypothetical ponds and compressed generations;',
                  'not a fitted model of an observed salamander population.'),
    marker_ascertainment = paste('Biallelic JC69 sites; pooled retained-sample MAF >= 0.05;',
      'candidate sites >= 10 kb apart; deterministic genome-wide 2500-site panel.'),
    cv_scope = 'Conditional on the globally ascertained retained marker panel.'
  )
)
stopifnot(identical(rownames(ohio_ponds$snp), samples$individual),
          identical(rownames(ohio_ponds$microsat), samples$individual),
          identical(rownames(ohio_ponds$raw_snp), raw_samples$individual),
          identical(rownames(ohio_ponds$raw_microsat), raw_samples$individual),
          all(table(factor(samples$pond, levels = ponds$pond)) == 20),
          identical(dim(ohio_ponds$snp), c(480L, 2500L)),
          length(ohio_ponds$snp_positions) == 2500L,
          all(ohio_ponds$snp %in% 0:2),
          all(diff(ohio_ponds$snp_positions) >= 10000),
          all(ponds$forest_buffer_100m >= 0.9),
          all(ponds$forest_buffer_300m >= 0.6))
dir.create('data', showWarnings = FALSE)
save(ohio_ponds, file = 'data/ohio_ponds.rda', compress = 'xz', version = 2)
print(c(bytes = file.info('data/ohio_ponds.rda')$size,
        ponds = nrow(ponds), retained = nrow(samples), raw = nrow(raw_samples)))
print(table(raw_samples$pond))
print(tools::md5sum('data/ohio_ponds.rda'))
