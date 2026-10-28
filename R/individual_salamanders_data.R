#' Simulated Individual Salamander Genotypes and Resistance Analyses
#'
#' Matched SNP and microsatellite genotypes for 160 adult terrestrial salamanders
#' in each of two simulated movement scenarios on a real southeastern Ohio
#' landscape. The primary example demonstrates nonlinear isolation by resistance
#' with individual sampling. The animals and drying index are hypothetical;
#' neither species occurrence nor field-calibrated movement is claimed.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{\code{data}}{Two named scenarios, \code{compensatory_500} and
#'     \code{control_500}. Each contains \code{samples}, \code{snp_dosage},
#'     \code{microsat_alleles}, and \code{covariance}, described below.}
#'   \item{\code{models}}{Twenty checked model records per scenario: uniform
#'     and quadratic drying models, two marker types, and partitions 0 (all
#'     adults) through 4 (one spatial quadrant withheld). Names combine marker,
#'     model, and partition, for example \code{SNP2500__IBR__0}. Records contain
#'     log-conductance coefficients \code{theta}, graph covariance multiplier
#'     \code{phi$tau} on its natural nonnegative scale, log nugget variance
#'     \code{phi$sigma}, starting-point diagnostics, convergence, curvature,
#'     covariance score, RMSE, and predicted-mean range. These are plain lists,
#'     not fitted \code{terradish} objects.}
#'   \item{\code{metadata}}{Fixed drying center \code{center}, linear divisor
#'     \code{sd_linear}, and quadratic divisor \code{sd_quadratic}; the latter
#'     is SD of squared deviations, not squared SD. \code{start_halfwidth}
#'     defines additional starting points, not coefficient bounds.
#'     \code{epsg}, \code{native_resolution_m}, \code{graph_resolution_m},
#'     \code{directions}, \code{extent}, \code{core}, and \code{graph_nodes}
#'     describe the spatial grid (EPSG:26917, meters). \code{sample_size},
#'     \code{main_SNPs}, and \code{matched_microsatellites} give panel sizes.
#'     \code{objective}, \code{biological_scope}, and
#'     \code{raw_marker_selection} state interpretation and selection rules.}
#' }
#'
#' @details
#' Within each scenario, \code{samples} has 160 rows in genotype order.
#' \code{sample_id} and \code{sample_key} identify the animal;
#' \code{id}, \code{mom}, and \code{dad} are recorded ecological pedigree
#' identifiers. \code{stratum} gives one of 16 sampling strata;
#' \code{quadrant} and \code{fold} identify the four spatial validation blocks.
#' \code{age} is in years, \code{sex} is 0 for female and 1 for male,
#' \code{x_local} and \code{y_local} are meters from the domain origin, and
#' \code{easting} and \code{northing} are projected UTM coordinates.
#' \code{drying} is the interpolated hypothetical index from 0 to 1.
#'
#' \code{snp_dosage} contains nested integer matrices named \code{500},
#' \code{1000}, and \code{2500}. Entries 0, 1, and 2 count one allele.
#' Locus names identify synthetic genomic blocks and within-block positions.
#' \code{microsat_alleles} is a 160 by 48 integer matrix; adjacent columns hold
#' the two repeat-count alleles at each of 24 loci. All marker row names match
#' \code{samples$sample_id}. Microsatellite allele labels are categorical,
#' not continuous response values.
#'
#' \code{covariance} contains centered matrices named \code{SNP2500} and
#' \code{microsat24}, averaging unstandardized allele-frequency inner products
#' across loci. Their implied squared distances range from 0 to 2. This is
#' different from the locus-frequency weighting in
#' \code{\link{cov_from_biallelic}} and is not a Dps response.
#'
#' Each scenario follows 500 years, with maximum age 20 for both sexes.
#' The compensatory scenario allows faster, straighter movement at intermediate
#' drying but greater danger in very dry terrain. The control holds behavioral
#' drying at 0.5 while retaining forest settlement and demography. It is a
#' control for drying-dependent movement, not every possible resistance effect.
#' Samples occupy distinct 10 m cells, are at least 10 m apart, and exclude
#' direct or shared recorded parents and grandparents. Distant kin can remain.
#'
#' These analyses use \code{wishart_covariance} with \code{nu = 1} as a
#' covariance quasi-likelihood, not a literal independent-Gaussian-locus
#' sampling model. Nominal Hessian uncertainty, AIC, and significance tests
#' are not used as genetic-sampling evidence. Two hundred genomic-block
#' resamples and a second genomic assay check conditional marker variation;
#' neither supplies independent ecological replication. Spatial folds were
#' reused during workflow development. Marker discovery conditions on the
#' full sampled panel. The matched microsatellites retain mixed predictive
#' evidence. No separate causal isolation-by-environment effect was simulated.
#'
#' Installed files in \code{extdata/individual-salamanders} provide the 10 m
#' drying terms and environment, native 30 m terrain, all 40 original prediction
#' records, comparative scores, the independent assay, and provenance. The
#' 10 m grid interpolates native 30 m information; it adds no new measurements.
#' The vignette explains these files and the longer opt-in reproduction.
#'
#' @source Simulated genotypes from recorded individual pedigrees on terrain
#'   derived from USGS 3DEP elevation and 2023 Annual NLCD Collection 1.0.
#'   See \code{system.file("extdata", "individual-salamanders", "README.md",
#'   package = "terradish")} for data sources, processing, and limitations.
#' @section Worked examples:
#' \code{vignette("data-generation", package = "terradish")},
#' \code{vignette("individual-salamanders", package = "terradish")}.
#'
#' @seealso \code{\link{cov_from_genetic_data}},
#'   \code{\link{cov_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{vignette("individual-salamanders", package = "terradish")}.
#' @examples
#' data(individual_salamanders)
#' animals <- individual_salamanders$data$compensatory_500
#' G <- animals$snp_dosage[["2500"]]
#' S <- cov_from_genetic_data(G, scale = FALSE, normalize = "none") /
#'   (2 * ncol(G))
#' dim(S)                              # 160 individual sampling locations
#' range(as.matrix(dist_from_cov(S)))   # squared allele-frequency distance
#' table(animals$samples$fold)           # 40 withheld adults per spatial fold
#' @usage data(individual_salamanders)
#' @name individual_salamanders
#' @docType data
#' @keywords datasets
NULL
