#' Wishart measurement models with site-level covariance kernels
#'
#' Creates a Wishart measurement model whose covariance structure includes the
#' resistance-implied covariance plus one or more fixed positive semidefinite
#' kernels built from site-level environmental covariates.
#'
#' @param x Site-level covariates.  Supported inputs are the same as
#'   \code{\link{pairwise_endpoint_covariates}}: a numeric vector, matrix, data
#'   frame, or \code{terra::SpatRaster}.
#' @param coords Required when \code{x} is a raster.  Focal-point coordinates
#'   in the same projection as \code{x}; accepts the same inputs as
#'   \code{\link{conductance_surface}}.
#' @param model Which Wishart likelihood should be used?  Use
#'   \code{"wishart_covariance"} for covariance-matrix responses and
#'   \code{"generalized_wishart"} for distance-matrix responses.
#' @param scale Logical.  Standardize site-level covariates before constructing
#'   kernels?  Recommended when covariates use different units.
#' @param normalize Logical.  If \code{TRUE}, rescale each kernel by its mean
#'   diagonal so the estimated kernel weights are on comparable scales.
#'
#' @details
#' \code{wishart_covariates()} is the Wishart analogue of adding endpoint
#' covariates to an MLPE mean model, but the covariates enter differently.  The
#' Wishart models are full-matrix likelihoods, so site-level covariates are
#' represented as covariance components rather than pairwise regression terms.
#' For each site covariate \eqn{z_k}, this helper centers the vector and builds
#' the positive semidefinite kernel \eqn{K_k = z_k z_k^\top}.  The fitted
#' covariance becomes
#'
#' \deqn{\Sigma = \tau E + \sum_k \lambda_k K_k + \exp(\sigma) I}
#'
#' where \eqn{E} is the resistance-implied covariance, \eqn{\tau} is the
#' resistance weight, \eqn{\lambda_k} are nonnegative kernel weights, and
#' \eqn{\exp(\sigma)} is the nugget variance.  With
#' \code{model = "generalized_wishart"}, the same \eqn{\Sigma} is projected to
#' the implied distance matrix by the generalized Wishart likelihood.
#'
#' This design keeps \eqn{\Sigma} positive definite when the nugget is positive
#' and all covariance-component weights are nonnegative.  It is intentionally
#' different from \code{\link{pairwise_endpoint_covariates}}, whose pairwise
#' dissimilarities are regression covariates for \code{\link{mlpe_covariates}}
#' rather than covariance kernels.
#'
#' @return A function of class \code{"terradish_measurement_model"} suitable
#'   for the \code{measurement_model} argument of \code{\link{terradish}} and
#'   \code{\link{terradish_grid}}.  The returned function stores the kernel
#'   array in attribute \code{"kernel_covariates"} and supports site-subsetting
#'   through \code{\link{terradish_cv}}.
#'
#' @seealso \code{\link{wishart_covariance}},
#'   \code{\link{generalized_wishart}}, \code{\link{mlpe_covariates}}
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' g <- wishart_covariates(melip.altitude, coords = melip.coords,
#'                         model = "generalized_wishart", scale = TRUE)
#' inherits(g, "terradish_measurement_model")
#'
#' @export
wishart_covariates <- function(x,
                               coords = NULL,
                               model = c("wishart_covariance",
                                         "generalized_wishart"),
                               scale = FALSE,
                               normalize = TRUE)
{
  model <- match.arg(model)
  site_covariates <- .pairwise_site_covariates(x, coords = coords, scale = scale)
  kernels <- .make_wishart_kernel_covariates(site_covariates,
                                             normalize = normalize)

  g <- switch(
    model,
    wishart_covariance = .wishart_covariate_model(kernels, covariance = TRUE),
    generalized_wishart = .wishart_covariate_model(kernels, covariance = FALSE)
  )

  attr(g, "base_model") <- model
  attr(g, "kernel_covariates") <- kernels
  attr(g, "subsetter") <- function(index)
    wishart_covariates(attr(kernels, "site_covariates")[index, , drop = FALSE],
                       model = model,
                       scale = FALSE,
                       normalize = normalize)
  class(g) <- unique(c("terradish_wishart_covariate_model",
                       "terradish_measurement_model",
                       "radish_measurement_model",
                       class(g)))
  g
}

.make_wishart_kernel_covariates <- function(site_covariates, normalize = TRUE)
{
  site_covariates <- as.matrix(site_covariates)
  if (!is.numeric(site_covariates))
    stop("site-level covariates must be numeric")
  if (anyNA(site_covariates))
    stop("missing values are not supported in Wishart kernel covariates")
  if (nrow(site_covariates) < 2L)
    stop("need at least two focal points to construct Wishart kernel covariates")
  if (is.null(colnames(site_covariates)))
    colnames(site_covariates) <- paste0("var", seq_len(ncol(site_covariates)))

  centered <- scale(site_covariates, center = TRUE, scale = FALSE)
  kernels <- array(0, dim = c(nrow(centered), nrow(centered), ncol(centered)),
                   dimnames = list(rownames(centered), rownames(centered),
                                   paste0("kernel_", colnames(centered))))
  scales <- rep(NA_real_, ncol(centered))
  names(scales) <- dimnames(kernels)[[3]]

  for (j in seq_len(ncol(centered)))
  {
    z <- centered[, j]
    K <- tcrossprod(z)
    if (isTRUE(normalize))
    {
      scale_j <- mean(diag(K))
      if (!is.finite(scale_j) || scale_j <= .Machine$double.eps)
        stop("Wishart kernel covariate `", colnames(centered)[j],
             "` has no variation.", call. = FALSE)
      K <- K / scale_j
      scales[j] <- scale_j
    }
    kernels[, , j] <- K
  }

  structure(kernels,
            site_covariates = site_covariates,
            kernel = "linear",
            normalized = isTRUE(normalize),
            kernel_scale = scales,
            class = unique(c("terradish_wishart_kernel_covariates",
                             class(kernels))))
}

.wishart_covariate_model <- function(kernels, covariance)
{
  force(kernels)
  force(covariance)

  function(E, S, phi, nu,
           gradient = TRUE,
           hessian = TRUE,
           partial = TRUE,
           nonnegative = TRUE,
           validate = FALSE)
  {
    .wishart_covariate_fit(E = E, S = S,
                           phi = if (missing(phi)) NULL else phi,
                           nu = nu, kernels = kernels,
                           covariance = covariance,
                           gradient = gradient,
                           hessian = hessian,
                           partial = partial,
                           nonnegative = nonnegative,
                           validate = validate)
  }
}

.wishart_covariate_fit <- function(E, S, phi, nu, kernels, covariance,
                                   gradient = TRUE, hessian = TRUE,
                                   partial = TRUE, nonnegative = TRUE,
                                   validate = FALSE)
{
  if (!(is.matrix(E) && is.matrix(S) && all(dim(E) == dim(S))))
    stop("invalid inputs", call. = FALSE)
  if (anyNA(E) || anyNA(S))
    stop("missing values are not supported", call. = FALSE)
  if (!is.array(kernels) || length(dim(kernels)) != 3L ||
      !all(dim(kernels)[1:2] == dim(E)))
    stop("Wishart kernel covariates do not match the response matrix.",
         call. = FALSE)

  E <- .pair_subset_symm(E)
  S <- .pair_subset_symm(S)
  if (!isTRUE(covariance))
  {
    if (any(diag(S) != 0))
      warning("Ignoring non-zero diagonal entries in `S`.")
    diag(S) <- 0
    if (any(S < 0))
      warning("Some distances are negative after symmetrization.")
  }

  kernel_names <- dimnames(kernels)[[3]]
  if (is.null(kernel_names))
    kernel_names <- paste0("kernel_", seq_len(dim(kernels)[3]))
  lambda_names <- paste0("lambda_", sub("^kernel_", "", kernel_names))
  phi_names <- c("tau", lambda_names, "sigma")

  if (is.null(phi))
  {
    if (isTRUE(covariance))
    {
      X <- cbind(c(E),
                 do.call(cbind, lapply(seq_len(dim(kernels)[3]),
                                       function(k) c(kernels[, , k]))),
                 c(diag(nrow(E))))
      coef0 <- tryCatch(qr.solve(X, c(S)),
                        error = function(e) rep(1e-6, length(phi_names)))
      coef0 <- pmax(as.numeric(coef0), 1e-6)
      phi <- c(coef0[-length(coef0)], log(coef0[length(coef0)]))
    }
    else
    {
      phi <- c(1, rep(1e-6, length(lambda_names)), 0)
    }
    names(phi) <- phi_names
    return(list(phi = phi,
                lower = c(0, rep(0, length(lambda_names)), -Inf),
                upper = rep(Inf, length(phi_names))))
  }

  if (!(is.numeric(phi) && length(phi) == length(phi_names)))
    stop("invalid inputs", call. = FALSE)
  if (anyNA(phi))
    stop("missing values are not supported", call. = FALSE)
  stopifnot(nu > 0)

  names(phi) <- phi_names
  tau <- phi["tau"]
  lambdas <- phi[lambda_names]
  sigma <- exp(phi["sigma"])
  stopifnot(tau >= 0, all(lambdas >= 0))
  nonnegative <- TRUE

  I <- diag(nrow(E))
  kernel_sum <- matrix(0, nrow(E), ncol(E))
  for (k in seq_along(lambdas))
    kernel_sum <- kernel_sum + lambdas[k] * kernels[, , k]
  Sigma <- tau * E + kernel_sum + sigma * I

  B <- c(list(tau = E),
         lapply(seq_along(lambdas), function(k) kernels[, , k]),
         list(sigma = sigma * I))
  names(B) <- phi_names

  if (isTRUE(covariance))
    .wishart_covariate_covariance(E = E, S = S, phi = phi, nu = nu,
                                  Sigma = Sigma, B = B, tau = tau,
                                  gradient = gradient, hessian = hessian,
                                  partial = partial,
                                  nonnegative = nonnegative,
                                  validate = validate)
  else
    .wishart_covariate_generalized(E = E, S = S, phi = phi, nu = nu,
                                   Sigma = Sigma, B = B, tau = tau,
                                   sigma = sigma,
                                   gradient = gradient, hessian = hessian,
                                   partial = partial,
                                   nonnegative = nonnegative,
                                   validate = validate)
}

.wishart_covariate_covariance <- function(E, S, phi, nu, Sigma, B, tau,
                                          gradient, hessian, partial,
                                          nonnegative, validate)
{
  symm <- function(X) (X + t(X)) / 2
  A <- solve(Sigma)
  ASA <- A %*% S %*% A
  objective <- nu / 2 * (as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) +
                           sum(diag(A %*% S)))
  grad_Sigma <- nu / 2 * (A - ASA)

  dgrad_from_dSigma <- function(U)
  {
    U <- symm(U)
    out <- nu / 2 * (-A %*% U %*% A +
                       A %*% U %*% ASA +
                       ASA %*% U %*% A)
    symm(out)
  }

  .wishart_covariate_finish(objective = c(objective),
                            fitted = Sigma,
                            grad_Sigma = grad_Sigma,
                            dgrad_from_dSigma = dgrad_from_dSigma,
                            E = E, S = S, phi = phi, B = B, tau = tau,
                            nu = nu, sigma_sign = 1,
                            gradient = gradient, hessian = hessian,
                            partial = partial, nonnegative = nonnegative,
                            validate = validate,
                            jacobian_S_factor = -nu / 2,
                            A_left = A, A_right = A,
                            jacobian_S_left = A,
                            jacobian_S_right = A)
}

.wishart_covariate_generalized <- function(E, S, phi, nu, Sigma, B, tau,
                                           sigma, gradient, hessian, partial,
                                           nonnegative, validate)
{
  symm <- function(X) (X + t(X)) / 2
  ones <- matrix(1, nrow(S), 1)
  I <- diag(nrow(S))
  SigmaInv <- solve(Sigma)
  SigInvOne <- SigmaInv %*% ones
  W <- I - ones %*% solve(t(ones) %*% SigInvOne) %*% t(SigInvOne)
  SigInvW <- SigmaInv %*% W
  eigSigW <- eigen(SigInvW)
  P <- eigSigW$vectors[, -nrow(Sigma)]
  D <- diag(eigSigW$values[-nrow(Sigma)])

  loglik <- nu / 4 * sum(diag(SigInvW %*% S)) + nu / 2 * sum(log(diag(D)))
  fitted <- diag(Sigma) %*% t(ones) + ones %*% t(diag(Sigma)) - 2 * Sigma
  grad_Sigma <- -nu / 2 * SigInvW - nu / 4 * SigInvW %*% S %*% t(SigInvW)

  dgrad_from_dSigma <- function(U)
  {
    U <- symm(U)
    dSigInvW <- -SigmaInv %*% U %*% SigInvW
    denom <- solve(t(ones) %*% SigmaInv %*% ones)
    dW_part <- ones %*% denom %*% t(ones) %*% SigmaInv %*% U %*% SigmaInv -
      c(t(ones) %*% SigmaInv %*% U %*% SigmaInv %*% ones) *
      c(denom^2) * ones %*% t(ones) %*% SigmaInv

    out <- -0.5 * nu * dSigInvW %*% W -
      0.25 * nu * dSigInvW %*% S %*% t(SigInvW) -
      0.5 * nu * W %*% dSigInvW -
      0.25 * nu * SigInvW %*% S %*% t(dSigInvW) +
      0.5 * nu * W %*% dSigInvW %*% W
    out <- -0.5 * nu * SigmaInv %*% dW_part -
      nu / 4 * SigmaInv %*% dW_part %*% S %*% t(W) %*% SigmaInv -
      nu / 4 * SigmaInv %*% W %*% S %*% t(dW_part) %*% SigmaInv + out
    out
  }

  .wishart_covariate_finish(objective = -c(loglik),
                            fitted = fitted,
                            grad_Sigma = grad_Sigma,
                            dgrad_from_dSigma = dgrad_from_dSigma,
                            E = E, S = S, phi = phi, B = B, tau = tau,
                            nu = nu, sigma_sign = -1,
                            gradient = gradient, hessian = hessian,
                            partial = partial, nonnegative = nonnegative,
                            validate = validate,
                            jacobian_S_factor = -nu / 4,
                            A_left = SigInvW, A_right = t(SigInvW),
                            jacobian_S_left = t(SigInvW),
                            jacobian_S_right = SigInvW)
}

.wishart_covariate_finish <- function(objective, fitted, grad_Sigma,
                                      dgrad_from_dSigma, E, S, phi, B, tau, nu,
                                      sigma_sign, gradient, hessian, partial,
                                      nonnegative, validate, jacobian_S_factor,
                                      A_left, A_right,
                                      jacobian_S_left, jacobian_S_right)
{
  p <- length(phi)
  dPhi <- matrix(0, p, 1, dimnames = list(names(phi), NULL))
  ddPhi <- matrix(0, p, p, dimnames = list(names(phi), names(phi)))
  ddEdPhi <- matrix(0, length(E), p,
                    dimnames = list(NULL, names(phi)))
  ddPhidS <- matrix(0, p, sum(lower.tri(S, diag = sigma_sign > 0)),
                    dimnames = list(names(phi), NULL))

  if (gradient || hessian || partial)
  {
    for (nm in names(phi))
      dPhi[nm, ] <- sum(B[[nm]] * grad_Sigma)

    dgrad <- NULL
    if (hessian || partial)
    {
      dgrad <- lapply(B, dgrad_from_dSigma)
      for (i in names(phi))
        for (j in names(phi))
          ddPhi[i, j] <- sum(B[[i]] * dgrad[[j]])
      ddPhi["sigma", "sigma"] <- ddPhi["sigma", "sigma"] + dPhi["sigma", ]
      ddPhi <- (ddPhi + t(ddPhi)) / 2

      if (partial)
      {
        dE <- tau * grad_Sigma
        for (nm in names(phi))
          ddEdPhi[, nm] <- c(tau * dgrad[[nm]])
        ddEdPhi[, "tau"] <- c(grad_Sigma + tau * dgrad[["tau"]])

        for (nm in names(phi))
        {
          dparam_dS <- jacobian_S_factor * A_left %*% B[[nm]] %*% A_right
          if (sigma_sign > 0)
            ddPhidS[nm, ] <- dparam_dS[lower.tri(dparam_dS, diag = TRUE)]
          else
            ddPhidS[nm, ] <- dparam_dS[lower.tri(dparam_dS)]
        }

        jacobian_E <- function(dotdotE)
        {
          U <- .pair_subset_symm(dotdotE)
          sigma_sign * tau^2 * dgrad_from_dSigma(U)
        }

        jacobian_S <- function(dotdotE)
        {
          U <- .pair_subset_symm(dotdotE)
          sigma_sign * jacobian_S_factor * tau *
            jacobian_S_left %*% U %*% jacobian_S_right
        }
      }
    }
  }

  list(objective = objective,
       fitted = fitted,
       boundary = nonnegative && tau == 0,
       gradient = if (!gradient) NULL else sigma_sign * dPhi,
       hessian = if (!hessian) NULL else sigma_sign * ddPhi,
       gradient_E = if (!partial) NULL else sigma_sign * dE,
       partial_E = if (!partial) NULL else sigma_sign * ddEdPhi,
       partial_S = if (!partial) NULL else sigma_sign * ddPhidS,
       jacobian_E = if (!partial) NULL else jacobian_E,
       jacobian_S = if (!partial) NULL else jacobian_S,
       num_gradient = NULL,
       num_hessian = NULL,
       num_gradient_E = NULL,
       num_partial_E = NULL,
       num_partial_S = NULL,
       num_jacobian_E = NULL,
       num_jacobian_S = NULL)
}
