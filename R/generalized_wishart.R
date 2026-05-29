#' Generalized Wishart distance regression
#'
#' A function of class \code{"terradish_measurement_model"} that evaluates the
#' generalized Wishart likelihood for an observed pairwise genetic distance
#' matrix.  This is the principled Wishart-based model when data are available
#' as distances (e.g. F\eqn{_{ST}}) and the effective marker degrees of freedom
#' are known.  For covariance-matrix data use \code{\link{wishart_covariance}}.
#'
#' @param E Conductance-implied covariance matrix: the generalized inverse of
#'   the graph Laplacian at the current conductance parameters.  Passed
#'   automatically by the optimizer.
#' @param S Square, symmetric matrix of observed pairwise genetic distances
#'   (e.g. F\eqn{_{ST}}). Must have the same dimensions as \code{E}.
#' @param phi Named numeric vector of nuisance parameters \code{(tau, sigma)}.
#'   Omit to obtain default starting values \code{c(1, 0)}.
#' @param nu Positive integer.  Effective Wishart degrees of freedom for
#'   \code{S}.  For biallelic SNPs this is usually the number of retained SNPs.
#'   For microsatellites, use the independent allele-frequency count,
#'   approximately \eqn{\sum_l (K_l - 1)} where \eqn{K_l} is the number of
#'   observed alleles at locus \eqn{l}; this is usually larger than the number
#'   of microsatellite loci and smaller than the total expanded allele-column
#'   count.  Must be supplied; it is not estimated.  Pass it via the \code{nu}
#'   argument of \code{\link{terradish}}.
#' @param gradient Logical. Compute gradient of the negative log-likelihood
#'   with respect to \code{phi}?
#' @param hessian Logical. Compute Hessian with respect to \code{phi}?
#' @param partial Logical. Compute second partial derivatives with respect to
#'   \code{phi}, \code{E}, and \code{S}?
#' @param nonnegative Unused; present for interface consistency.
#' @param validate Logical. Numerically validate gradients and Hessians via
#'   \pkg{numDeriv}? Very slow; for debugging small examples only.
#'
#' @details
#' The nuisance parameters are:
#' \describe{
#'   \item{\code{tau}}{Nonnegative scale applied to the conductance-implied
#'     covariance \code{E}. A \code{tau} near zero signals no detectable
#'     resistance effect.}
#'   \item{\code{sigma}}{Log-scale identity component: the nugget variance
#'     added to the diagonal is \eqn{\exp(\sigma)}. It absorbs genetic
#'     variation not explained by landscape resistance.}
#' }
#'
#' The fitted covariance is \eqn{\Sigma = \tau E + \exp(\sigma) I}.  The model
#' evaluates the generalized Wishart log-likelihood for the observed
#' squared-distance matrix \code{S} after projecting out the grand mean, as
#' described in McCullagh (2009) and Peterson et al. (2019).
#'
#' \code{generalized_wishart} and \code{\link{wishart_covariance}} share the
#' same \eqn{\Sigma} parameterization but differ in what \code{S} represents:
#' this function takes a \strong{distance} matrix; \code{wishart_covariance}
#' takes a \strong{covariance} matrix.  When the covariance and its implied
#' distance matrix carry the same information (i.e. no information is lost in
#' the conversion), the two models give identical likelihoods up to a constant.
#'
#' @references
#' McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica
#' Sinica 19:23-41.
#'
#' Peterson EK, Peterman WE, Pope NS. 2019. resist_ga: An R package for
#' landscape genetic resistance surface optimization using genetic algorithms.
#' Methods in Ecology and Evolution 10(9):1502-1509.
#'
#' @seealso \code{\link{wishart_covariance}}, \code{\link{mlpe}},
#'   \code{\link{terradish}}
#'
#' @return When \code{phi} is missing, a list with elements \code{phi}
#'   (starting values \code{c(tau = 1, sigma = 0)}), \code{lower}
#'   (\code{c(0, -Inf)}), and \code{upper}.  Otherwise a list containing:
#'  \item{objective}{Negative log-likelihood.}
#'  \item{fitted}{Matrix of expected genetic distances.}
#'  \item{boundary}{Logical; \code{TRUE} if \code{tau = 0}, indicating no detectable IBR signal.}
#'  \item{gradient}{Gradient with respect to \code{phi} (if \code{gradient = TRUE}).}
#'  \item{hessian}{Hessian matrix with respect to \code{phi} (if \code{hessian = TRUE}).}
#'  \item{gradient_E}{Gradient with respect to \code{E} (if \code{partial = TRUE}).}
#'  \item{partial_E}{Jacobian of \code{gradient_E} with respect to \code{phi} (if \code{partial = TRUE}).}
#'  \item{partial_S}{Jacobian of \code{gradient} with respect to the lower triangle of \code{S} (if \code{partial = TRUE}).}
#'  \item{jacobian_E}{Function for reverse-mode AD through \code{E} (if \code{partial = TRUE}).}
#'  \item{jacobian_S}{Function for reverse-mode AD through \code{S} (if \code{partial = TRUE}).}
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
#' generalized_wishart(laplacian_inv, melip.Fst, nu = 1000, phi = c(0.1, -0.1))
#'
#' @export

generalized_wishart <- function(E, S, phi, nu, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE)
{
  symm <- function(X) (X + t(X))/2

  if (missing(phi)) #return starting values and boundaries for optimization of phi
  {
    return(list(phi = c(1, 0), lower = c(0, -Inf), upper = c(Inf, Inf)))
  }
  else if (!(is.matrix(E)    & 
             is.matrix(S)    & 
             all(dim(E)  == dim(S)) &
             is.numeric(phi) & 
             length(phi) == 2 ))
    stop ("invalid inputs")
  if (anyNA(E) || anyNA(S) || anyNA(phi))
    stop("missing values are not supported")

  stopifnot(nu > 0)

  # Symmetrize the observed distances and ignore any diagonal input.
  S <- symm(S)
  if (any(diag(S) != 0))
    warning("Ignoring non-zero diagonal entries in `S`.")
  diag(S) <- 0
  if (any(S < 0))
    warning("Some distances are negative after symmetrization.")

  names(phi) <- c("tau", "sigma")
  tau   <- phi["tau"]
  sigma <- exp(phi["sigma"])

  # density is undefined if tau is negative
  stopifnot(tau >= 0)
  nonnegative <- TRUE

  ones      <- matrix(1, nrow(S), 1)
  I         <- diag(nrow(S))
  Sigma     <- tau * E + sigma * I
  SigInvOne <- solve(Sigma, ones)

  W        <- I - ones %*% solve(t(ones) %*% SigInvOne) %*% t(SigInvOne)
  SigInvW  <- solve(Sigma, W)
  eigSigW  <- eigen(SigInvW)
  P        <- eigSigW$vectors[,-nrow(Sigma)]
  D        <- diag(eigSigW$values[-nrow(Sigma)])
  ginvSigW <- P %*% solve(D) %*% t(P)

  loglik <- nu/4 * sum(diag(SigInvW %*% S)) + nu/2 * sum(log(diag(D)))

  fitted <- diag(Sigma) %*% t(ones) + ones %*% t(diag(Sigma)) - 2 * Sigma

  if (gradient || hessian || partial)
  {
    dPhi    <- matrix(0, length(phi), 1)
    ddPhi   <- matrix(0, length(phi), length(phi))
    ddEdPhi <- matrix(0, length(E),   length(phi))
    ddPhidS <- matrix(0, length(phi), sum(lower.tri(S)))
    rownames(dPhi) <- colnames(ddPhi) <- 
      rownames(ddPhi) <- colnames(ddEdPhi) <- 
        rownames(ddPhidS) <- names(phi)

    grad_Sigma <- -nu/2 * SigInvW - nu/4 * SigInvW %*% S %*% t(SigInvW)

    # gradient, phi
    dPhi["tau",]   <- sum(E * grad_Sigma) 
    dPhi["sigma",] <- sum(diag(grad_Sigma)) * sigma

    if (hessian || partial)
    {
      # see Golub GH, Pereyra V. 1973. The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate. SIAM Journal on Numerical Analysis 10(2): 413-432
      # ^^actually unnecessary
      dSigInvW_dtau <- -solve(Sigma, E %*% SigInvW)
      dSigInvW_dsigma <- -solve(Sigma, SigInvW)

      fuckoff_tau <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% E %*% solve(Sigma) +
        -c(t(ones) %*% solve(Sigma) %*% E %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
      #dgrad_dtau <- -0.5 * nu * dSigInvW_dtau - 0.25 * nu * dSigInvW_dtau %*% S %*% t(SigInvW) -
      #               0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dtau)
      dgrad_dtau <- -0.5 * nu * dSigInvW_dtau %*% W - 0.25 * nu * dSigInvW_dtau %*% S %*% t(SigInvW) -
                       0.5 * nu * W %*% dSigInvW_dtau - 0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dtau) +
                       0.5 * nu * W %*% dSigInvW_dtau %*% W
      dgrad_dtau <- -0.5 * nu * solve(Sigma) %*% fuckoff_tau - 
        nu/4 * solve(Sigma) %*% fuckoff_tau %*% S %*% t(W) %*% solve(Sigma) -
        nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_tau) %*% solve(Sigma) + dgrad_dtau
      
      fuckoff_sigma <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% solve(Sigma) +
        -c(t(ones) %*% solve(Sigma) %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
      dgrad_dsigma <- -0.5 * nu * dSigInvW_dsigma %*% W - 0.25 * nu * dSigInvW_dsigma %*% S %*% t(SigInvW) -
                       0.5 * nu * W %*% dSigInvW_dsigma - 0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dsigma) +
                       0.5 * nu * W %*% dSigInvW_dsigma %*% W
      dgrad_dsigma <- -0.5 * nu * solve(Sigma) %*% fuckoff_sigma - 
        nu/4 * solve(Sigma) %*% fuckoff_sigma %*% S %*% t(W) %*% solve(Sigma) -
        nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_sigma) %*% solve(Sigma) + dgrad_dsigma
      dgrad_dsigma <- dgrad_dsigma * sigma

      # hessian, phi x phi
      ddPhi["tau","tau"] <- sum(E * dgrad_dtau)
      ddPhi["tau","sigma"] <- sum(E * dgrad_dsigma)
      ddPhi["sigma","sigma"] <- sigma * sum(diag(dgrad_dsigma)) + dPhi["sigma",]
      ddPhi <- ddPhi + t(ddPhi)
      diag(ddPhi) <- diag(ddPhi)/2

      if(partial)
      {
        # gradient wrt E
        dE <- tau * grad_Sigma

        # hessian offdiagonal, E x phi
        ddEdPhi[,"tau"] <- grad_Sigma + tau * dgrad_dtau
        ddEdPhi[,"sigma"] <- dgrad_dsigma * tau

        # hessian offdiagonal, S x phi
        ddtaudS <- -nu/4 * SigInvW %*% E %*% t(SigInvW)
        ddsigmadS <- -nu/4 * SigInvW %*% t(SigInvW) * sigma
        ddPhidS["tau",] <- ddtaudS[lower.tri(ddtaudS)]
        ddPhidS["sigma",] <- ddsigmadS[lower.tri(ddsigmadS)]

        # jacobian products (label these properly)
        jacobian_E <- function(dE)
        {
          dE <- symm(dE)
          dSigInvW_dE <- -solve(Sigma, dE %*% SigInvW)
          fuckoff_E <- ones %*% solve(t(ones) %*% solve(Sigma) %*% ones) %*% t(ones) %*% solve(Sigma) %*% dE %*% solve(Sigma) +
            -c(t(ones) %*% solve(Sigma) %*% dE %*% solve(Sigma) %*% ones) * c(solve(t(ones) %*% solve(Sigma) %*% ones)^2) * ones %*% t(ones) %*% solve(Sigma)
          dgrad_dE <- -0.5 * nu * dSigInvW_dE - 0.25 * nu * dSigInvW_dE %*% S %*% t(SigInvW) -
            0.25 * nu * SigInvW %*% S %*% t(dSigInvW_dE)
          dgrad_dE <- -0.5 * nu * solve(Sigma) %*% fuckoff_E - 
            nu/4 * solve(Sigma) %*% fuckoff_E %*% S %*% t(W) %*% solve(Sigma) -
            nu/4 * solve(Sigma) %*% W %*% S %*% t(fuckoff_E) %*% solve(Sigma) + dgrad_dE
          -dgrad_dE * tau^2
        } #jesus christ

        jacobian_S <- function(dE)
        {
          dE <- symm(dE)
          dEdS <- -nu/4 * t(SigInvW) %*% dE %*% SigInvW * tau
          diag(dEdS) <- 0
          -dEdS
        } 
      }
    }
  }

  if (validate)
  {
    arg <- num_gradient <- .numderiv_jacobian(function(x) 
                                   generalized_wishart(E = E, 
                                        phi = x, 
                                        nu = nu,
                                        S = S)$grSig, 
                                   phi)
    num_gradient <- .numderiv_grad(function(x) 
                                   generalized_wishart(E = E, 
                                        phi = x, 
                                        nu = nu,
                                        S = S)$objective, 
                                   phi)

    num_hessian <- .numderiv_jacobian(function(x) 
                                     generalized_wishart(E = E, 
                                          phi = x, 
                                          nu = nu,
                                          S = S)$gradient, 
                                     phi)

    num_gradient_E <- symm(matrix(.numderiv_grad(function(x) 
                                                 generalized_wishart(E = x, 
                                                      phi = phi, 
                                                      nu = nu,
                                                      S = S)$objective, 
                                                 E), 
                                  nrow(E), ncol(E)))

    num_partial_E <- .numderiv_jacobian(function(x) 
                                        generalized_wishart(E = E, 
                                             phi = x, 
                                             nu = nu,
                                             S = S)$gradient_E, 
                                        phi)

    num_partial_S <- .numderiv_jacobian(function(x) 
                                        generalized_wishart(E = E, 
                                             phi = phi, 
                                             nu = nu,
                                             S = x)$gradient, 
                                        S)[,lower.tri(S)]

    num_jacobian_E <- function(X) 
      matrix(c(X) %*% .numderiv_jacobian(function(x) 
                                         generalized_wishart(E = x, 
                                              phi = phi, 
                                              nu = nu,
                                              S = S)$gradient_E, 
                                         E), 
             nrow(X), ncol(X))

    num_jacobian_S <- function(X) 
      matrix(c(X) %*% .numderiv_jacobian(function(x) 
                                         generalized_wishart(E = E, 
                                              phi = phi, 
                                              nu = nu,
                                              S = x)$gradient_E, 
                                         S), 
             nrow(X), ncol(X))
  }

  list(objective        = -c(loglik), 
       fitted           = fitted,
       boundary         = nonnegative && tau == 0,
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
class(generalized_wishart) <- c("terradish_measurement_model",
                                "radish_measurement_model")
