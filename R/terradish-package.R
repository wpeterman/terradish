#' @title terradish: Fast Gradient-Based Optimization of Resistance Surfaces
#'
#' @description
#' Fits associations between raster covariates and population-genetic distance
#' or covariance by maximum likelihood. You supply raster covariates, sampling
#' locations, and an observed genetic matrix. \pkg{terradish} maps covariates to
#' vertex conductance, builds a weighted graph Laplacian, derives its generalized
#' inverse \eqn{E(\theta)}, and links that shared process layer to the response
#' through a selected measurement likelihood.
#'
#' @details
#' \if{html}{
#' \figure{terradish-sticker.png}{options: width=150 alt='Package logo'}
#'}
#'
#' \strong{Where to start.}  \code{\link{conductance_surface}} builds the
#' graph, \code{\link{terradish}} fits the model, and the methods documented in
#' \code{\link{terradish_methods}} (\code{summary()}, \code{coef()},
#' \code{plot()}, \code{anova()}) read the result.  Conductance parameters are
#' on the log scale, so \code{exp(coef(fit))} is the multiplicative change in
#' conductance per unit of a covariate.
#' These coefficients describe conditional associations within the fitted
#' graph and candidate set. They do not, by themselves, identify causal
#' mechanisms, absolute migration or dispersal rates, or habitat suitability.
#'
#' \strong{Main components.}
#' \describe{
#'   \item{Conductance models}{\code{\link{loglinear_conductance}} (the
#'     default), \code{\link{linear_conductance}},
#'     \code{\link{smooth_loglinear_conductance}} for spline responses, and
#'     \code{\link{gaussian_smoothed_loglinear_conductance}} for estimating a
#'     covariate's Gaussian raster-smoothing scale inside the model.}
#'   \item{Measurement models}{\code{\link{leastsquares}} and
#'     \code{\link{mlpe}} for a genetic distance matrix;
#'     \code{\link{generalized_wishart}} and \code{\link{wishart_covariance}}
#'     for a Wishart likelihood, which need the effective degrees of freedom
#'     \code{nu}; \code{\link{check_distance_response}} to verify a
#'     generalized-Wishart squared-distance response; \code{\link{mlpe_covariates}} and
#'     \code{\link{wishart_covariates}} to add isolation-by-environment terms.}
#'   \item{Model extensions}{\code{\link{terradish_hierarchical}} adds a smooth
#'     diagnostic residual field for spatial conductance structure not captured
#'     by supplied covariates. \code{\link{terradish_directed}} fits
#'     antisymmetric edge-rate bias through a symmetric commute-time response.}
#'   \item{Comparison and diagnostics}{\code{\link{aic_table}},
#'     \code{\link{terradish_cv}}, \code{\link{terradish_grid}}, and
#'     \code{\link{terradish_assess_settings}}, which profiles the graph and
#'     recommends solver and optimizer settings.}
#' }
#'
#' \strong{Vignettes.}  Start with
#' \code{vignette("getting-started", package = "terradish")}.  The others cover
#' model comparison, joint isolation by environment and resistance, the
#' covariance-based Wishart workflow, spline conductance, scale optimization,
#' hierarchical and directional models, large landscapes, and simulation-based
#' study design.  List them all with
#' \code{browseVignettes("terradish")}.
#'
#' @seealso \code{\link{terradish}}, \code{\link{conductance_surface}},
#'   \code{\link{terradish_methods}}, \code{\link{melip}}
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom Matrix Cholesky Diagonal forceSymmetric rowSums solve sparseMatrix t update
#' @importFrom methods as new setRefClass
#' @importFrom multiScaleR kernel_scale.raster
#' @importFrom nlme gls
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply stopCluster
#' @importFrom stats AIC D anova as.dist as.formula coef cov2cor delete.response dist
#' @importFrom stats dnorm fft fitted formula lm logLik model.matrix optim optimize pchisq plogis
#' @importFrom stats pnorm printCoefmat prcomp qlogis qnorm reformulate residuals rnorm
#' @importFrom stats rWishart sd setNames sigma simulate terms
#' @importFrom terra adjacent cellFromXY crop geom global is.factor levels ncell nlyr patches rast
#' @importFrom terra ext extract is.lonlat res rowColFromCell unwrap values values<- xyFromCell
#' @importFrom utils globalVariables modifyList write.csv
#' @importFrom grDevices terrain.colors
#' @useDynLib terradish, .registration = TRUE
#' @md
"_PACKAGE"

globalVariables(c(
  "abs_log_rate_ratio", "conductance",
  "distance", "distance_lower", "distance_upper", "est", "label",
  "label_y", "observed", "upper", "weight", "x", "xend", "y", "yend"
))
