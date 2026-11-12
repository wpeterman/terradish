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
#' \strong{Parameter names and their scales.} Several names recur with
#' different meanings. Read each one on the scale of the model that reports
#' it.
#' \describe{
#'   \item{\code{theta}}{Conductance parameters, returned by \code{coef()}.
#'     Log-linear and spline coefficients are on the log-conductance scale, so
#'     \code{exp(theta)} is a conductance multiplier per unit of the
#'     covariate. Gaussian \code{sigma.<layer>} entries of \code{theta} are
#'     smoothing scales in map units.}
#'   \item{\code{phi}}{Nuisance parameters of the measurement model, reported
#'     by \code{summary()}. Their names and scales depend on the model:
#'     \code{alpha}, \code{beta}, \code{tau}, \code{rho}, \code{sigma}, and
#'     environmental coefficients.}
#'   \item{\code{tau}}{In \code{\link{mlpe}} and \code{\link{leastsquares}},
#'     the log residual precision: the residual variance is \code{exp(-tau)}.
#'     In \code{\link{wishart_covariance}}, \code{\link{generalized_wishart}},
#'     and \code{\link{simulate_covariance_response}}, the nonnegative
#'     multiplier on the graph covariance, on its natural scale.}
#'   \item{\code{sigma}}{In fitted Wishart models, the log nugget variance:
#'     the nugget is \code{exp(sigma)}. In
#'     \code{\link{simulate_covariance_response}} and
#'     \code{\link{covariance_response_power}}, the nugget variance itself, on
#'     its natural scale. In
#'     \code{\link{gaussian_smoothed_loglinear_conductance}}, a Gaussian
#'     smoothing scale in map units (\code{sigma.<layer>}). In
#'     \code{\link{HagerZhangControl}}, a line-search constant unrelated to
#'     any model.}
#'   \item{\code{rho}}{In \code{\link{mlpe}}, the shared-site correlation on a
#'     logit scale: the correlation is \code{plogis(rho) / 2}, between 0 and
#'     0.5. \code{rho_start} in \code{\link{mlpe_covariates}} is on the
#'     natural scale. In \code{\link{HagerZhangControl}}, an unrelated
#'     line-search constant.}
#'   \item{\code{nu}}{Effective Wishart degrees of freedom, supplied by the
#'     user and never estimated. It scales the log-likelihood, so standard
#'     errors scale as \code{1 / sqrt(nu)} while point estimates do not
#'     change. It is not the marker count. Unused by the distance models.}
#'   \item{\code{alpha}, \code{beta}}{In the distance models, the intercept
#'     and the slope of genetic distance on resistance distance, in the units
#'     of the response.}
#'   \item{\code{lambda_<covariate>}}{In \code{\link{wishart_covariates}}, a
#'     nonnegative weight on an environmental covariance kernel.}
#' }
#'
#' \strong{Moving from radish.} The original \code{radish*} function names
#' still work but are deprecated; see \code{\link{radish_compatibility}}.
#'
#' \strong{Experimental features.} Directed and hierarchical models,
#' site-specific drift terms, selected-pair likelihoods, and other research
#' prototypes remain on the \code{experimental} branch. Their identifiability
#' and validation limits are described in that branch's \code{EXPERIMENTAL.md}.
#' Install it explicitly with
#' \code{remotes::install_github("wpeterman/terradish@experimental")}.
#'
#' \strong{Vignettes.}  Start with
#' \code{vignette("v01-getting-started", package = "terradish")}, or with
#' \code{vignette("v03-mlpe-pairwise-distance", package = "terradish")} if your
#' data are a pairwise genetic distance matrix.  The others cover
#' model comparison, joint isolation by environment and resistance, the
#' covariance-based Wishart workflow, spline conductance, scale optimization,
#' large landscapes, and simulation-based
#' study design.  List them all with
#' \code{browseVignettes("terradish")}.
#'
#' @section Worked examples:
#' \code{vignette("v12-data-generation", package = "terradish")},
#' \code{vignette("v01-getting-started", package = "terradish")},
#' \code{vignette("v09-individual-salamanders", package = "terradish")},
#' \code{vignette("v03-mlpe-pairwise-distance", package = "terradish")}.
#'
#' @references
#' Peterman WE, Pope NS (2021). The use and misuse of regression models in
#' landscape genetic analyses. Molecular Ecology, 30, 37--47.
#' \doi{10.1111/mec.15716}.
#'
#' Clarke RT, Rothery P, Raybould AF (2002). Confidence limits for regression
#' relationships between distance matrices: estimating gene flow with distance.
#' Journal of Agricultural, Biological, and Environmental Statistics, 7,
#' 361--372. \doi{10.1198/108571102320}.
#'
#' Hanks EM, Hooten MB (2013). Circuit theory and model-based inference for
#' landscape connectivity. Journal of the American Statistical Association,
#' 108, 22--33. \doi{10.1080/01621459.2012.724647}.
#'
#' @seealso \code{\link{terradish}}, \code{\link{conductance_surface}},
#'   \code{\link{terradish_methods}}, \code{\link{melip}},
#'   \code{\link{radish_compatibility}}
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom Matrix Cholesky Diagonal forceSymmetric rowSums solve sparseMatrix t update
#' @importFrom methods as new setRefClass
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
