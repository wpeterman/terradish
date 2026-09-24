# Re-exported from landgraph (the shared base package).
#
# The genetic-covariance and distance helpers now live
# in landgraph and are re-exported here so existing terradish workflows and
# documentation links point to landgraph's authoritative help. Grouped covariance
# should use the coherent gower diagonal; within-population diagonals have a
# different scale. Rare-variant standardization and unequal group sizes need
# sensitivity checks. FST ratio estimators are not Wishart responses. pca_dist
# replaces missing entries with means. These caveats also belong in the
# covariance and getting-started vignettes.

#' @importFrom landgraph cov_from_biallelic
#' @export
landgraph::cov_from_biallelic

#' @importFrom landgraph cov_from_genetic_data
#' @export
landgraph::cov_from_genetic_data

#' @importFrom landgraph fst_from_biallelic
#' @export
landgraph::fst_from_biallelic

#' @importFrom landgraph dist_from_cov
#' @export
landgraph::dist_from_cov

#' @importFrom landgraph dist_from_biallelic
#' @export
landgraph::dist_from_biallelic
