# Re-exported from landgraph (the shared base package).
#
# The genetic-covariance, distance, and directional-edge-covariate helpers now live
# in landgraph and are re-exported here so existing terradish workflows and
# documentation links are unchanged. terradish's own functions (resistance, Tier-1/2/3,
# wishart measurement models, simulation) continue to call them as before.

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

#' @importFrom landgraph edge_gradient
#' @export
landgraph::edge_gradient

#' @importFrom landgraph edge_flow
#' @export
landgraph::edge_flow
