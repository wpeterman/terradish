#' Least-squares measurement model
#'
#' A function of class \code{"terradish_measurement_model"} that evaluates a
#' Gaussian likelihood with independent errors and a linear mean structure
#' relating observed genetic distances to resistance distances.  This is the
#' fastest measurement model and a useful first-pass choice; prefer
#' \code{\link{mlpe}} for inference because it correctly accounts for the
#' non-independence of pairwise measurements.
#'
#' @param E Conductance-implied covariance matrix: the generalized inverse of
#'   the graph Laplacian, evaluated at the current conductance parameters.
#'   Passed automatically by the optimizer; users normally do not call this
#'   function directly.
#' @param S Square, symmetric matrix of observed pairwise genetic distances
#'   (e.g. F\eqn{_{ST}}). Must have the same dimensions as \code{E}.
#' @param phi Named numeric vector of nuisance parameters \code{(alpha, beta,
#'   tau)}.  Omit to obtain maximum-likelihood starting values from an
#'   \code{nlme::gls} fit.
#' @param nu Unused; present for a common interface with Wishart measurement
#'   models.
#' @param gradient Logical. Compute gradient of the negative log-likelihood
#'   with respect to \code{phi}?
#' @param hessian Logical. Compute Hessian of the negative log-likelihood with
#'   respect to \code{phi}?
#' @param partial Logical. Compute second partial derivatives with respect to
#'   \code{phi}, \code{E}, and \code{S}? Required by the optimizer; set
#'   \code{FALSE} only for standalone likelihood evaluations.
#' @param nonnegative Logical. Constrain the IBR slope \code{beta} to be
#'   nonnegative? Default \code{TRUE} prevents a nonsensical negative
#'   resistance-distance effect.
#' @param validate Logical. Numerically validate gradients and Hessians via
#'   \pkg{numDeriv}? Very slow; intended for debugging small examples only.
#'
#' @details
#' The nuisance parameters are:
#' \describe{
#'   \item{\code{alpha}}{Intercept of the regression of genetic distance on
#'     resistance distance.}
#'   \item{\code{beta}}{Slope (IBR effect); constrained \eqn{\geq 0} when
#'     \code{nonnegative = TRUE}.}
#'   \item{\code{tau}}{Log-precision of the Gaussian errors: residual variance
#'     is \eqn{\exp(-\tau)}.}
#' }
#' The mean structure is \eqn{S_{ij} = \alpha + \beta R_{ij} + e_{ij}}, where
#' \eqn{R_{ij}} is the resistance distance between sites \eqn{i} and \eqn{j}
#' (derived from \code{E}) and \eqn{e_{ij}} are independent Gaussian errors
#' with common precision \eqn{\exp(\tau)}.
#'
#' Pairwise genetic distances are not independent: any two pairs sharing a
#' sampling site are correlated.  \code{leastsquares} ignores this, which
#' underestimates standard errors.  For inferential purposes, \code{\link{mlpe}}
#' is strongly preferred.  Use \code{leastsquares} when speed matters more than
#' precision, or for initial parameter exploration.
#'
#' @seealso \code{\link{mlpe}}, \code{\link{generalized_wishart}},
#'   \code{\link{wishart_covariance}}, \code{\link{terradish}}
#'
#' @return When \code{phi} is missing, a list with elements \code{phi}
#'   (starting values), \code{lower}, and \code{upper} (parameter bounds).
#'   Otherwise a list containing:
#'  \item{objective}{Negative log-likelihood.}
#'  \item{fitted}{Matrix of expected genetic distances (same dimensions as \code{S}).}
#'  \item{boundary}{Logical; \code{TRUE} if the MLE is on the boundary (\code{beta = 0}), indicating no detectable IBR signal.}
#'  \item{gradient}{Gradient of the negative log-likelihood with respect to \code{phi} (if \code{gradient = TRUE}).}
#'  \item{hessian}{Hessian matrix with respect to \code{phi} (if \code{hessian = TRUE}).}
#'  \item{gradient_E}{Gradient with respect to \code{E} (if \code{partial = TRUE}).}
#'  \item{partial_E}{Jacobian of \code{gradient_E} with respect to \code{phi} (if \code{partial = TRUE}).}
#'  \item{partial_S}{Jacobian of \code{gradient} with respect to the lower triangle of \code{S} (if \code{partial = TRUE}).}
#'  \item{jacobian_E}{Function for reverse-mode algorithmic differentiation through \code{E} (if \code{partial = TRUE}).}
#'  \item{jacobian_S}{Function for reverse-mode algorithmic differentiation through \code{S} (if \code{partial = TRUE}).}
#'
#' @examples
#'
#' library(terra)
#' 
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#' 
#' covariates <- c(melip.altitude, melip.forestcover)
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' # inverse of graph Laplacian at null model (IBD) 
#' laplacian_inv <- terradish_distance(theta = matrix(0, 1, 2), 
#'                                  formula = ~forestcover + altitude,
#'                                  data = surface,
#'                                  terradish::loglinear_conductance, 
#'                                  covariance = TRUE)$covariance[,,1]
#' 
#' leastsquares(laplacian_inv, melip.Fst) #without 'phi': return MLE of phi
#' leastsquares(laplacian_inv, melip.Fst, phi = c(0., 0.5, -0.1))
#'
#' @export

leastsquares <- function(E, S, phi, nu = NULL, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  symm <- function(X) (X + t(X))/2

  if (missing(phi)) #return starting values and boundaries for optimization
  {
    ones <- matrix(1, nrow(E), 1)
    Ed   <- diag(E)
    R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
    Rl   <- R[lower.tri(R)]
    Sl   <- S[lower.tri(S)]
    fit  <- gls(Sl ~ Rl, method = "ML")

    if (!nonnegative || coef(fit)[2] > 0) 
    {
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = phi[2], "tau" = -2 * log(sigma(fit)))
    }
    else
    {
      fit <- gls(Sl ~ 1, method = "ML")
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = 0, "tau" = -2 * log(sigma(fit)))
    }

    return(list(phi = phi, 
                lower = if(nonnegative) c(-Inf, 0, -Inf) else c(-Inf, -Inf, -Inf), 
                upper = c(Inf, Inf, Inf)))
  }
  else if (!(is.matrix(E)    & 
             is.matrix(S)    & 
             all(dim(E)  == dim(S)) &
             is.numeric(phi) & 
             length(phi) == 3 ))
    stop ("invalid inputs")

  names(phi) <- c("alpha", "beta", "tau")

  alpha <- phi["alpha"]
  tau   <- exp(phi["tau"])
  beta  <- phi["beta"]

  ones <- matrix(1, nrow(E), 1)
  Ed   <- diag(E)
  R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
  Rl   <- R[lower.tri(R)]
  Sl   <- S[lower.tri(S)]

  unos   <- matrix(1, length(Sl), 1)
  e      <- Sl - alpha * unos - beta * Rl
  loglik <- -0.5 * tau * t(e) %*% e + 0.5 * nrow(e) * log(tau)

  distance <- matrix(0, nrow(S), ncol(S))
  distance[lower.tri(distance)] <- Rl
  distance <- distance + t(distance)
  fitted <- alpha + beta * distance

  # gradients, hessians, mixed partial derivatives
  if (gradient || hessian || partial)
  {
    dPhi    <- matrix(0, length(phi), 1)
    ddPhi   <- matrix(0, length(phi), length(phi))
    ddEdPhi <- matrix(0, length(Rl),  length(phi))
    ddPhidS <- matrix(0, length(phi), length(Sl))
    rownames(dPhi) <- colnames(ddPhi) <- 
      rownames(ddPhi) <- colnames(ddEdPhi) <- 
        rownames(ddPhidS) <- names(phi)

    # gradient, phi
    dPhi["alpha",] <- t(unos) %*% e * tau
    dPhi["beta",]  <- t(Rl) %*% e * tau
    dPhi["tau",]   <- -0.5 * tau * t(e) %*% e + 0.5 * length(e)

    if (hessian || partial)
    {
      # hessian, phi x phi
      ddPhi["alpha", "alpha"] <- -t(unos) %*% unos * tau
      ddPhi["alpha",  "beta"] <- -tau * t(unos) %*% Rl
      ddPhi["alpha",   "tau"] <- tau * t(unos) %*% e
      ddPhi[ "beta",  "beta"] <- -t(Rl) %*% Rl * tau
      ddPhi[ "beta",   "tau"] <- tau * t(Rl) %*% e
      ddPhi[  "tau",   "tau"] <- -0.5 * tau * t(e) %*% e
      ddPhi                   <- ddPhi + t(ddPhi)
      diag(ddPhi)             <- diag(ddPhi)/2

      if (partial)
      {
        # gradient wrt E
        dR <- matrix(0, nrow(R), ncol(R))
        dR[lower.tri(dR)] <- 2 * beta * tau * as.vector(e)
        dR <- symm(dR)
        dE <- diag(nrow(R)) * (dR %*% ones %*% t(ones)) - dR

        # hessian offdiagonal, E x phi
        ddEdPhi[, "alpha"] <- -2 * beta * tau * unos
        ddEdPhi[,  "beta"] <- 2 * tau * (e - beta * Rl)
        ddEdPhi[,   "tau"] <- 2 * beta * tau * e
        ddEdPhi            <- apply(ddEdPhi, 2, function(x) { X <- matrix(0,nrow(E),ncol(E)); X[lower.tri(X)] <- x; X <- symm(X); diag(nrow(E)) * (X %*% ones %*% t(ones)) - X })

        # hessian offdiagonal, S x phi
        ddPhidS["alpha",] <- unos * tau
        ddPhidS["beta",]  <- Rl * tau
        ddPhidS["tau",]   <- -tau * e

        # jacobian products (label these properly)
        jacobian_E <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE <- -beta^2 * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          -ddEdE
        }

        jacobian_S <- function(dE)
        {
          ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
          ddEdE <- -beta * tau * ddEdE
          ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
          diag(ddEdE) <- 0
          -ddEdE
        }
      }
    }
  }

  if (validate)
  {
    num_gradient    <- .numderiv_grad(function(x) leastsquares(E = E, phi = x, S = S)$objective, phi)
    num_hessian     <- .numderiv_hessian(function(x) leastsquares(E = E, phi = x, S = S)$objective, phi)
    num_gradient_E  <- symm(matrix(.numderiv_grad(function(x) leastsquares(E = x, phi = phi, S = S)$objective, E), nrow(E), ncol(E)))
    num_partial_E   <- .numderiv_jacobian(function(x) leastsquares(E = E, phi = x, S = S)$gradient_E, phi)
    num_partial_S   <- .numderiv_jacobian(function(x) leastsquares(E = E, phi = phi, S = x)$gradient, S)[,lower.tri(S)]
    num_jacobian_E  <- function(X) matrix(c(X) %*% .numderiv_jacobian(function(x) leastsquares(E = x, phi = phi, S = S)$gradient_E, E), nrow(X), ncol(X))
    num_jacobian_S  <- function(X) matrix(c(X) %*% .numderiv_jacobian(function(x) leastsquares(E = E, phi = phi, S = x)$gradient_E, S), nrow(X), ncol(X))
  }

  list(objective        = -c(loglik), 
       fitted           = fitted,
       boundary         = nonnegative && beta == 0,
       gradient         = if(!gradient) NULL else -dPhi,
       hessian          = if(!hessian)  NULL else -ddPhi,
       gradient_E       = if(!partial)  NULL else -dE, 
       partial_E        = if(!partial)  NULL else -ddEdPhi,   # partial_E[i,k] is d(dl/dE_i)/dPhi_k where i is linearized matrix index
       partial_S        = if(!partial)  NULL else -ddPhidS,   # partial_S[i,k] is d(dl/dPhi_i)/dS_k where k is linearized matrix index
       jacobian_E       = if(!partial)  NULL else jacobian_E, # function mapping vectorized dg/dE to d(dg/dE)/dE
       jacobian_S       = if(!partial)  NULL else jacobian_S, # function mapping vectorized dg/dE to d(dg/dE)/dS
       num_gradient     = if(!validate) NULL else num_gradient,
       num_hessian      = if(!validate) NULL else num_hessian,
       num_gradient_E   = if(!validate) NULL else num_gradient_E,
       num_partial_E    = if(!validate) NULL else num_partial_E,
       num_partial_S    = if(!validate) NULL else num_partial_S,
       num_jacobian_E   = if(!validate) NULL else num_jacobian_E,
       num_jacobian_S   = if(!validate) NULL else num_jacobian_S)
}
class(leastsquares) <- c("terradish_measurement_model",
                         "radish_measurement_model")
