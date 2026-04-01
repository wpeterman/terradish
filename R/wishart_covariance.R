#' Wishart covariance regression
#'
#' A function of class \code{measurement_model} that calculates likelihood,
#' gradient, hessian, and partial derivatives of nuisance parameters and the
#' Laplacian generalized inverse, using a Wishart model for an observed
#' covariance matrix.
#'
#' @param E A submatrix of the generalized inverse of the graph Laplacian (e.g.
#'   a covariance matrix implied by the resistance surface).
#' @param S A matrix of observed genetic covariance.
#' @param phi Nuisance parameters (see details).
#' @param nu Number of genetic markers.
#' @param gradient Compute gradient of negative loglikelihood with regard to
#'   \code{phi}?
#' @param hessian Compute Hessian matrix of negative loglikelihood with regard
#'   to \code{phi}?
#' @param partial Compute second partial derivatives of negative loglikelihood
#'   with regard to \code{phi}, \code{E}, \code{S}?
#' @param nonnegative Unused.
#' @param validate Numerical validation via package \code{numDeriv} (very slow,
#'   use for debugging small examples).
#'
#' @details The nuisance parameters are the scaling of the generalized inverse
#' of the graph Laplacian ("tau"; can be zero) and a log scalar multiple of the
#' identity matrix that is added to the generalized inverse ("sigma"). The
#' fitted covariance is \code{Sigma = tau * E + exp(sigma) * I}, and the model
#' evaluates a Wishart negative log-likelihood for the observed covariance
#' matrix \code{S}, treating \code{nu} as the effective number of loci.
#'
#' @seealso \code{\link{radish_measurement_model}},
#'   \code{\link{generalized_wishart}}
#'
#' @return A list containing:
#'  \item{objective}{the negative loglikelihood}
#'  \item{fitted}{the fitted covariance response matrix}
#'  \item{boundary}{is the MLE on the boundary (e.g. no spatial genetic structure)?}
#'  \item{gradient}{gradient of negative loglikelihood with respect to phi}
#'  \item{hessian}{Hessian matrix of negative loglikelihood with respect to phi}
#'  \item{gradient_E}{gradient with respect to the generalized inverse of the graph Laplacian}
#'  \item{partial_E}{Jacobian of \code{gradient_E} with respect to phi}
#'  \item{partial_S}{Jacobian of \code{gradient} with respect to the strict lower triangle of S}
#'  \item{jacobian_E}{a function used for reverse algorithmic differentiation}
#'  \item{jacobian_S}{a function used for reverse algorithmic differentiation}
#'
#' @examples
#' E <- matrix(c(1.2, 0.3, 0.2,
#'               0.3, 1.5, 0.4,
#'               0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
#' Sigma <- 0.8 * E + 0.2 * diag(3)
#' S <- rWishart(1, df = 20, Sigma = Sigma)[,,1] / 20
#' wishart_covariance(E, S, phi = c(0.8, log(0.2)), nu = 20)
#'
#' @export
wishart_covariance <- function(E, S, phi, nu,
                               gradient = TRUE,
                               hessian = TRUE,
                               partial = TRUE,
                               nonnegative = TRUE,
                               validate = FALSE)
{
  symm <- function(X) (X + t(X)) / 2

  if (missing(phi))
  {
    if (!(is.matrix(E) && is.matrix(S) && all(dim(E) == dim(S))))
      stop("invalid inputs")
    if (anyNA(E) || anyNA(S))
      stop("missing values are not supported")

    E <- symm(E)
    S <- symm(S)
    I <- diag(nrow(E))
    X <- cbind(c(E), c(I))
    y <- c(S)
    coef0 <- tryCatch(qr.solve(X, y),
                      error = function(e) c(1, mean(diag(S))))
    tau0 <- max(as.numeric(coef0[1]), 1e-6)
    sigma0 <- max(as.numeric(coef0[2]), 1e-6)

    return(list(phi = c("tau" = tau0, "sigma" = log(sigma0)),
                lower = c(0, -Inf),
                upper = c(Inf, Inf)))
  }
  else if (!(is.matrix(E) &&
             is.matrix(S) &&
             all(dim(E) == dim(S)) &&
             is.numeric(phi) &&
             length(phi) == 2))
    stop("invalid inputs")

  if (anyNA(E) || anyNA(S) || anyNA(phi))
    stop("missing values are not supported")
  stopifnot(nu > 0)

  E <- symm(E)
  S <- symm(S)

  names(phi) <- c("tau", "sigma")
  tau <- phi["tau"]
  sigma <- exp(phi["sigma"])

  stopifnot(tau >= 0)
  nonnegative <- TRUE

  I <- diag(nrow(E))
  Sigma <- tau * E + sigma * I
  SigmaInv <- solve(Sigma)
  A <- SigmaInv
  ASA <- A %*% S %*% A
  objective <- nu / 2 * (as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) +
                           sum(diag(A %*% S)))
  fitted <- Sigma

  if (gradient || hessian || partial)
  {
    dPhi <- matrix(0, length(phi), 1)
    ddPhi <- matrix(0, length(phi), length(phi))
    ddEdPhi <- matrix(0, length(E), length(phi))
    ddPhidS <- matrix(0, length(phi), sum(lower.tri(S, diag = TRUE)))
    rownames(dPhi) <- colnames(ddPhi) <-
      rownames(ddPhi) <- colnames(ddEdPhi) <-
        rownames(ddPhidS) <- names(phi)

    grad_Sigma <- nu / 2 * (A - ASA)

    dPhi["tau", ] <- sum(E * grad_Sigma)
    dPhi["sigma", ] <- sum(diag(grad_Sigma)) * sigma

    if (hessian || partial)
    {
      dgrad_from_dSigma <- function(B)
      {
        B <- symm(B)
        out <- nu / 2 * (-A %*% B %*% A +
                           A %*% B %*% ASA +
                           ASA %*% B %*% A)
        symm(out)
      }

      dgrad_dtau <- dgrad_from_dSigma(E)
      dgrad_dsigma <- dgrad_from_dSigma(sigma * I)

      ddPhi["tau", "tau"] <- sum(E * dgrad_dtau)
      ddPhi["tau", "sigma"] <- sum(E * dgrad_dsigma)
      ddPhi["sigma", "sigma"] <- sigma * sum(diag(dgrad_dsigma)) + dPhi["sigma", ]
      ddPhi <- ddPhi + t(ddPhi)
      diag(ddPhi) <- diag(ddPhi) / 2

      if (partial)
      {
        dE <- tau * grad_Sigma

        ddEdPhi[, "tau"] <- c(grad_Sigma + tau * dgrad_dtau)
        ddEdPhi[, "sigma"] <- c(tau * dgrad_dsigma)

        ddtaudS <- -nu / 2 * A %*% E %*% A
        ddsigmadS <- -nu / 2 * (A %*% A) * sigma
        ddPhidS["tau", ] <- ddtaudS[lower.tri(ddtaudS, diag = TRUE)]
        ddPhidS["sigma", ] <- ddsigmadS[lower.tri(ddsigmadS, diag = TRUE)]

        jacobian_E <- function(dotdotE)
        {
          U <- symm(dotdotE)
          out <- nu / 2 * tau^2 * (-A %*% U %*% A +
                                     ASA %*% U %*% A +
                                     A %*% U %*% ASA)
          symm(out)
        }

        jacobian_S <- function(dotdotE)
        {
          U <- symm(dotdotE)
          symm(-nu / 2 * tau * A %*% U %*% A)
        }
      }
    }
  }

  if (validate)
  {
    num_gradient <- .numderiv_grad(function(x)
      wishart_covariance(E = E, phi = x, nu = nu, S = S)$objective,
    phi)

    num_hessian <- .numderiv_hessian(function(x)
      wishart_covariance(E = E, phi = x, nu = nu, S = S)$objective,
    phi)

    num_gradient_E <- symm(matrix(.numderiv_grad(function(x)
      wishart_covariance(E = x, phi = phi, nu = nu, S = S)$objective,
    E), nrow(E), ncol(E)))

    num_partial_E <- .numderiv_jacobian(function(x)
      wishart_covariance(E = E, phi = x, nu = nu, S = S)$gradient_E,
    phi)

    num_partial_S <- .numderiv_jacobian(function(x)
      wishart_covariance(E = E, phi = phi, nu = nu, S = x)$gradient,
    S)[, lower.tri(S, diag = TRUE)]

    num_jacobian_E <- function(X)
      matrix(c(X) %*% .numderiv_jacobian(function(x)
        wishart_covariance(E = x, phi = phi, nu = nu, S = S)$gradient_E,
      E),
      nrow(X), ncol(X))

    num_jacobian_S <- function(X)
      matrix(c(X) %*% .numderiv_jacobian(function(x)
        wishart_covariance(E = E, phi = phi, nu = nu, S = x)$gradient_E,
      S),
      nrow(X), ncol(X))
  }

  list(objective = c(objective),
       fitted = fitted,
       boundary = nonnegative && tau == 0,
       gradient = if (!gradient) NULL else dPhi,
       hessian = if (!hessian) NULL else ddPhi,
       gradient_E = if (!partial) NULL else dE,
       partial_E = if (!partial) NULL else ddEdPhi,
       partial_S = if (!partial) NULL else ddPhidS,
       jacobian_E = if (!partial) NULL else jacobian_E,
       jacobian_S = if (!partial) NULL else jacobian_S,
       num_gradient = if (!validate) NULL else num_gradient,
       num_hessian = if (!validate) NULL else num_hessian,
       num_gradient_E = if (!validate) NULL else num_gradient_E,
       num_partial_E = if (!validate) NULL else num_partial_E,
       num_partial_S = if (!validate) NULL else num_partial_S,
       num_jacobian_E = if (!validate) NULL else num_jacobian_E,
       num_jacobian_S = if (!validate) NULL else num_jacobian_S)
}
class(wishart_covariance) <- "radish_measurement_model"
