#' @section Interpreting pairwise effects:
#' Construct one \code{\link{pairwise_covariates}} object to use the same
#' environmental transforms and geographic distances in both likelihood
#' families. Their fitted conductance surfaces are conditional on these terms,
#' and conductance signs can change when the mean structure changes. Report
#' fits with and without pairwise effects. A uniform \code{~ 1} conductance
#' model with pairwise terms represents IBD plus IBE, not IBE alone.
#'
#' A positive environmental coefficient means differentiation increases with
#' environmental difference conditional on resistance. It does not distinguish
#' assortative mating, dispersal filtering, local adaptation, or site variance
#' that rises at environmental extremes. MLPE coefficients are unconstrained;
#' Wishart kernel weights are nonnegative. Compare their effects in
#' resistance-distance units with \code{\link{terradish_ibe_ratio}}, which
#' divides by the resistance coefficient and propagates joint uncertainty.
#' Ratios are undefined when that denominator is zero.
