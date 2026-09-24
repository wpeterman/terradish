#' @title terradish: Fast Gradient-Based Optimization of Resistance Surfaces
#'
#' @description
#' Estimates landscape conductance surfaces from genetic data by maximum
#' likelihood on a raster graph. It fits log-linear conductance models and can
#' optionally estimate each covariate's Gaussian smoothing scale and spline
#' shape. Models fit pairwise genetic distances through the MLPE likelihood,
#' or allele-frequency covariance and squared distances derived from it through
#' a Wishart likelihood on site contrasts. Both families accept the same
#' pairwise environmental covariates alongside resistance. Exact sparse
#' Cholesky and algebraic multigrid solvers handle large rasters; fixed-domain
#' spatial cross-validation compares conductance formulas.
#'
#' @details
#' \if{html}{\figure{terradish-sticker.png}{options: width=150 alt='terradish logo'}}
#' \if{latex}{\figure{terradish-sticker.png}{options: width=1.5in}}
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
#' For Wishart models, the user-supplied effective degrees of freedom,
#' \code{nu}, sets inferential precision. Using the SNP count substantially
#' overstated information in the package's forward-time validation experiments.
#' Report sensitivity across plausible values with
#' \code{\link{terradish_rescale_nu}}.
#'
#' \strong{Main components.}
#' \describe{
#'   \item{Conductance models}{\code{\link{loglinear_conductance}} (the
#'     default),
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
#'   \item{Pairwise inputs and inference}{\code{\link{pairwise_covariates}}
#'     combines geographic distances and environmental differences;
#'     \code{\link{terradish_ibe_ratio}} expresses their fitted effects in
#'     resistance-distance units. \code{\link{terradish_inference}} provides
#'     covariance matrices and confidence intervals, and
#'     \code{\link{gaussian_scale_profile}} profiles smoothing-scale uncertainty.}
#'   \item{Comparison and diagnostics}{\code{\link{aic_table}},
#'     \code{\link{terradish_cv_folds}}, \code{\link{terradish_grid}}, and
#'     \code{\link{terradish_assess_settings}}, a speed probe for solver and
#'     optimizer settings. It does not establish statistical validity.}
#'   \item{Preparation and design}{\code{\link{scale_covariates}},
#'     \code{\link{crop_to_focal_buffer}},
#'     \code{\link{simulate_covariance_response}}, and
#'     \code{\link{covariance_response_power}} support reproducible inputs and
#'     model-based study design. \code{\link{slim_terradish}} reduces stored
#'     fit size when later prediction is unnecessary.}
#' }
#'
#' \strong{Experimental features.} Directed and hierarchical models,
#' site-specific drift terms, selected-pair likelihoods, and other research
#' prototypes remain on the \code{experimental} branch. Their identifiability
#' and validation limits are described in that branch's \code{EXPERIMENTAL.md}.
#' Install it explicitly with
#' \code{remotes::install_github("wpeterman/terradish@experimental")}.
#'
#' \strong{Vignettes.}  Start with
#' \code{vignette("getting-started", package = "terradish")}.  The others cover
#' model comparison, joint isolation by environment and resistance, the
#' covariance-based Wishart workflow, spline conductance, scale optimization,
#' large landscapes, and simulation-based
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
#' @importFrom terra aggregate adjacent cellFromXY crop geom global is.factor levels ncell nlyr patches rast
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
