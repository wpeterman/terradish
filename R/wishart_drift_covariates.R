#' Wishart measurement models with site-specific diagonal variance
#'
#' Creates a Wishart measurement model whose fitted covariance replaces the
#' single scalar nugget of \code{\link{wishart_covariance}} or
#' \code{\link{generalized_wishart}} with a \strong{per-site diagonal}
#' parameterized by site-level covariates. This models heterogeneity in diagonal genetic
#' variance. Under an explicit population-genetic model and compatible sampling
#' design, that heterogeneity may be consistent with differences in drift or
#' effective population size, but the fitted diagonal is not a direct estimate
#' of either quantity. It is the diagonal analogue of
#' \code{\link{wishart_covariates}}, which adds covariance kernels.
#'
#' @param x Site-level covariates that parameterize the diagonal variance surface. Supported
#'   inputs are the same as \code{\link{wishart_covariates}}: a numeric vector,
#'   matrix, data frame, or \code{terra::SpatRaster}.  \code{NULL} (the default)
#'   gives an intercept-only model, which reproduces the scalar-nugget
#'   \code{\link{wishart_covariance}} / \code{\link{generalized_wishart}}.
#' @param coords Required when \code{x} is a raster.  Focal-point coordinates in
#'   the same projection as \code{x}; accepts the same inputs as
#'   \code{\link{conductance_surface}}.
#' @param model Which Wishart likelihood to use.
#'   \code{"wishart_covariance"} is appropriate when the response \code{S} is a
#'   \strong{covariance} matrix (e.g. from \code{\link{cov_from_genetic_data}});
#'   \code{"generalized_wishart"} is appropriate only when \code{S} is an
#'   admissible squared-distance representation coherently related to a centered
#'   positive semidefinite covariance matrix. An arbitrary distance matrix,
#'   including an unchecked F\eqn{_{ST}} matrix, is not sufficient.
#' @param scale Logical.  Standardize site-level covariates to zero mean and
#'   unit variance before building the design matrix?  Recommended when
#'   covariates use different units or scales.
#'
#' @details
#' The fitted covariance is
#'
#' \deqn{\Sigma = \tau\, E(\theta) + \mathrm{diag}\!\big(\exp(Z\gamma)\big),}
#'
#' where \eqn{E(\theta)} is the resistance-implied covariance (generalized
#' inverse of the graph Laplacian at conductance parameters \eqn{\theta}),
#' \eqn{\tau \ge 0} is the resistance weight (IBR signal), and \eqn{Z} is the
#' site design matrix consisting of an intercept column plus one (mean-centered)
#' column per diagonal-variance covariate. The per-site nugget
#' \eqn{n_i = \exp((Z\gamma)_i)} is strictly positive for any real
#' \eqn{\gamma}, so \eqn{\Sigma} stays positive definite whenever \eqn{\tau \ge 0}.
#'
#' The nuisance parameters estimated alongside \eqn{\theta} are:
#' \describe{
#'   \item{\code{tau}}{Nonnegative scale on the resistance-implied covariance
#'     \eqn{E}.  A value near zero indicates no detectable IBR signal.}
#'   \item{\code{sigma}}{Intercept of the log-nugget: the baseline diagonal
#'     variance at mean covariate values. With no diagonal covariates this
#'     is exactly the scalar nugget of \code{\link{wishart_covariance}}.}
#'   \item{\code{gamma_<covariate>}}{Slope of the log-nugget on each supplied
#'     site covariate. \eqn{\gamma_j > 0} means that higher covariate values are
#'     associated with greater fitted diagonal variance, conditional on the
#'     graph covariance and other diagonal terms.}
#' }
#'
#' \strong{When to use this.} Use the extension when the response and study
#' design support a hypothesis about site-specific diagonal variance. It can
#' reduce confounding between a common nugget and structured diagonal
#' heterogeneity, but it does not guarantee separation of landscape
#' connectivity, drift, sampling variation, density, or effective size.
#'
#' \strong{Scope.} This is a deme-structured covariance parameterization: the
#' diagonal is a per-focal-site variance. Reading that term as genetic drift
#' requires additional assumptions about the genetic summary, marker sampling,
#' deme definition, and demographic model. Do not report it as local
#' \eqn{N_e}, abundance, carrying capacity, or density without an external
#' calibration that identifies that relationship.
#'
#' For the \code{"generalized_wishart"} (squared-distance) model the per-site
#' diagonal is a site-specific diagonal contribution to the admissible distance
#' representation. Its biological meaning remains conditional on the same
#' assumptions described above.
#'
#' @return A function of class \code{"terradish_measurement_model"} suitable for
#'   the \code{measurement_model} argument of \code{\link{terradish}} and
#'   \code{\link{terradish_grid}}.  The nuisance parameter vector \eqn{\phi}
#'   contains \code{tau}, \code{sigma} (the log-nugget intercept), and one
#'   \code{gamma_<covariate>} per drift covariate, in that order.  The function
#'   stores the diagonal covariates in attribute \code{"drift_covariates"} and
#'   supports site-subsetting for cross-validation through its \code{"subsetter"}
#'   attribute.
#'
#' @seealso \code{\link{check_distance_response}},
#'   \code{\link{wishart_covariance}}, \code{\link{generalized_wishart}},
#'   \code{\link{wishart_covariates}}, \code{\link{terradish}}
#'
#' @references
#' McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica
#' Sinica 19:631-649.
#'
#' @examples
#' # Intercept-only reduces to the scalar-nugget covariance Wishart.
#' E <- matrix(c(1.2, 0.3, 0.2,
#'               0.3, 1.5, 0.4,
#'               0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
#' Sigma <- 0.8 * E + 0.2 * diag(3)
#' S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25
#'
#' g0 <- wishart_drift_covariates(model = "wishart_covariance")
#' start <- g0(E, S, nu = 25)
#' names(start$phi)        # "tau"   "sigma"
#'
#' # With a site covariate associated with diagonal variance.
#' g <- wishart_drift_covariates(data.frame(local_variance = c(-1, 0, 1)),
#'                               model = "wishart_covariance")
#' fit <- g(E, S, phi = c(tau = 0.8, sigma = log(0.2),
#'                         gamma_local_variance = 0.1),
#'          nu = 25)
#' fit$objective
#'
#' @export
wishart_drift_covariates <- function(x = NULL,
                                     coords = NULL,
                                     model = c("wishart_covariance",
                                               "generalized_wishart"),
                                     scale = FALSE)
{
  model <- match.arg(model)
  site_covariates <- if (is.null(x))
    NULL
  else
    .pairwise_site_covariates(x, coords = coords, scale = scale)
  covariance <- identical(model, "wishart_covariance")

  g <- function(E, S, phi, nu,
                gradient = TRUE,
                hessian = TRUE,
                partial = TRUE,
                nonnegative = TRUE,
                validate = FALSE)
  {
    .wishart_drift_fit(E = E, S = S,
                       phi = if (missing(phi)) NULL else phi,
                       nu = nu,
                       site_covariates = site_covariates,
                       covariance = covariance,
                       gradient = gradient,
                       hessian = hessian,
                       partial = partial,
                       nonnegative = nonnegative,
                       validate = validate)
  }

  attr(g, "base_model") <- model
  attr(g, "drift_covariates") <- site_covariates
  attr(g, "subsetter") <- function(index)
    wishart_drift_covariates(
      if (is.null(site_covariates)) NULL else
        site_covariates[index, , drop = FALSE],
      model = model,
      scale = FALSE)
  class(g) <- unique(c("terradish_wishart_drift_model",
                       "terradish_measurement_model",
                       "radish_measurement_model"))
  g
}

# Build the per-site design matrix Z = [intercept, centered drift covariates].
# The intercept column is named "sigma" so that an intercept-only model has phi
# = (tau, sigma), identical to wishart_covariance.
.wishart_drift_design <- function(site_covariates, n)
{
  if (is.null(site_covariates))
  {
    Z <- matrix(1, n, 1L)
    colnames(Z) <- "sigma"
    return(Z)
  }

  site_covariates <- as.matrix(site_covariates)
  if (!is.numeric(site_covariates))
    stop("drift covariates must be numeric", call. = FALSE)
  if (anyNA(site_covariates))
    stop("missing values are not supported in drift covariates", call. = FALSE)
  if (nrow(site_covariates) != n)
    stop("drift covariates have ", nrow(site_covariates),
         " rows but the response matrix has ", n, " sites.", call. = FALSE)
  if (is.null(colnames(site_covariates)))
    colnames(site_covariates) <- paste0("var", seq_len(ncol(site_covariates)))

  centered <- scale(site_covariates, center = TRUE, scale = FALSE)
  Z <- cbind(1, centered)
  colnames(Z) <- c("sigma", paste0("gamma_", colnames(site_covariates)))
  Z
}

.wishart_drift_fit <- function(E, S, phi, nu, site_covariates, covariance,
                               gradient = TRUE, hessian = TRUE,
                               partial = TRUE, nonnegative = TRUE,
                               validate = FALSE)
{
  if (!(is.matrix(E) && is.matrix(S) && all(dim(E) == dim(S))))
    stop("invalid inputs", call. = FALSE)
  if (anyNA(E) || anyNA(S))
    stop("missing values are not supported", call. = FALSE)

  n <- nrow(E)
  Z <- .wishart_drift_design(site_covariates, n)
  diag_names <- colnames(Z)              # "sigma" plus "gamma_<cov>"
  phi_names <- c("tau", diag_names)

  E <- .pair_subset_symm(E)
  if (!isTRUE(covariance))
    S <- .prepare_gw_response(S)
  else
    S <- .pair_subset_symm(S)

  if (is.null(phi))
  {
    if (isTRUE(covariance))
    {
      X <- cbind(c(E), c(diag(n)))
      coef0 <- tryCatch(qr.solve(X, c(S)),
                        error = function(e) c(1, mean(diag(S))))
      tau0 <- max(as.numeric(coef0[1]), 1e-6)
      sigma0 <- max(as.numeric(coef0[2]), 1e-6)
      phi <- c(tau0, log(sigma0), rep(0, ncol(Z) - 1L))
    }
    else
    {
      phi <- c(1, 0, rep(0, ncol(Z) - 1L))
    }
    names(phi) <- phi_names
    return(list(phi = phi,
                lower = c(0, rep(-Inf, ncol(Z))),
                upper = rep(Inf, length(phi_names))))
  }

  if (!(is.numeric(phi) && length(phi) == length(phi_names)))
    stop("invalid inputs", call. = FALSE)
  if (anyNA(phi))
    stop("missing values are not supported", call. = FALSE)
  stopifnot(nu > 0)

  names(phi) <- phi_names
  tau <- phi["tau"]
  gamma <- phi[diag_names]
  stopifnot(tau >= 0)
  nonnegative <- TRUE

  n_vec <- as.vector(exp(Z %*% gamma))
  I_diag <- diag(n_vec)
  Sigma <- tau * E + I_diag

  # Per-parameter basis matrices B[[nm]] = dSigma/dphi_nm.
  #   B[[tau]]      = E
  #   B[[diag d]]   = diag(n_vec * Z[, d])   (since dn_i/dgamma_d = n_i Z_id)
  B <- c(list(tau = E),
         lapply(seq_len(ncol(Z)),
                function(d) diag(n_vec * Z[, d])))
  names(B) <- phi_names

  # Second-order ("curvature") term: contribution of d2Sigma/dphi_i dphi_j.
  # Only diagonal parameters are curved, with
  #   d2Sigma/dphi_a dphi_b = diag(n_vec * Z[,a] * Z[,b]),
  # so the curvature block is Z' diag(w) Z with w = n_vec * diag(grad_Sigma).
  # This reduces to the scalar-nugget special case when Z is intercept-only.
  curvature <- function(grad_Sigma)
  {
    w <- n_vec * diag(grad_Sigma)
    C <- matrix(0, length(phi_names), length(phi_names),
                dimnames = list(phi_names, phi_names))
    C[diag_names, diag_names] <- crossprod(Z, w * Z)
    C
  }

  if (isTRUE(covariance))
    .wishart_covariate_covariance(E = E, S = S, phi = phi, nu = nu,
                                  Sigma = Sigma, B = B, tau = tau,
                                  gradient = gradient, hessian = hessian,
                                  partial = partial,
                                  nonnegative = nonnegative,
                                  validate = validate,
                                  curvature = curvature)
  else
    .wishart_covariate_generalized(E = E, S = S, phi = phi, nu = nu,
                                   Sigma = Sigma, B = B, tau = tau,
                                   sigma = exp(phi["sigma"]),
                                   gradient = gradient, hessian = hessian,
                                   partial = partial,
                                   nonnegative = nonnegative,
                                   validate = validate,
                                   curvature = curvature)
}
