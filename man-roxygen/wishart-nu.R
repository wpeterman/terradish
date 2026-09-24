#' @section Effective Wishart information:
#' For Wishart likelihoods, \code{nu} scales the log likelihood. Point
#' estimates do not depend on it, but standard errors, likelihood-ratio tests,
#' and information criteria do. It is a dispersion parameter that you set,
#' rather than the number of markers. Linkage and the population history
#' shared by loci can make effective information much smaller than the SNP
#' count. In forward-time simulations used to validate terradish, standard
#' errors at the SNP count were three to five times smaller than the spread
#' among replicate histories. Those model-based errors understated uncertainty.
#'
#' Use \code{\link{terradish_cv_folds}} to compare conductance formulas, and
#' \code{\link{terradish_rescale_nu}} to report inference across plausible
#' values. For the same successful folds and nuisance settings, changing
#' \code{nu} rescales Wishart CV scores without changing their ranking.
#' Neither procedure identifies the true effective information or repairs
#' an unsuitable response matrix.
