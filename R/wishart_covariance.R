#' Wishart covariance measurement model
#'
#' A function of class \code{"terradish_measurement_model"} that evaluates a
#' Wishart likelihood for an observed genetic \strong{covariance} matrix.
#' This is the natural model when raw genotype data are summarized with
#' \code{\link{cov_from_genetic_data}} or \code{\link{cov_from_biallelic}}
#' rather than as pairwise distances.  For distance-matrix data use
#' \code{\link{generalized_wishart}}.
#'
#' @param E Conductance-implied covariance matrix: the generalized inverse of
#'   the graph Laplacian at the current conductance parameters.  Passed
#'   automatically by the optimizer; users normally do not call this function
#'   directly.
#' @param S Square, symmetric matrix of observed population-level genetic
#'   covariance.  Typically the output of \code{\link{cov_from_genetic_data}}
#'   or \code{\link{cov_from_biallelic}}.  Must have the same dimensions as
#'   \code{E} and should be positive (semi-)definite.
#' @param phi Named numeric vector of nuisance parameters \code{(tau, sigma)}.
#'   Omit to obtain least-squares starting values.
#' @param nu Positive number.  Effective Wishart degrees of freedom for the
#'   genetic covariance \code{S}: the number of (approximately) independent
#'   genetic markers that went into computing \code{S}.  Must be supplied; it
#'   is not estimated.  Pass via the \code{nu} argument of
#'   \code{\link{terradish}}.  See the \dQuote{The role of \code{nu}} section
#'   below; choosing \code{nu} is consequential and is explained in detail
#'   there.
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
#'     covariance \code{E}.  A value near zero signals no detectable resistance
#'     effect.}
#'   \item{\code{sigma}}{Log-scale nugget: the identity component added to the
#'     model covariance is \eqn{\exp(\sigma)}.  It absorbs genetic drift,
#'     population-size differences, and any variation not explained by landscape
#'     resistance.}
#' }
#'
#' The fitted covariance is \eqn{\Sigma = \tau E + \exp(\sigma) I}.  The
#' negative log-likelihood is:
#'
#' \deqn{-\ell = \frac{\nu}{2} \left[ \log|\Sigma| + \mathrm{tr}(\Sigma^{-1} S) \right]}
#'
#' treating \code{S} as a sample covariance with \eqn{\nu} degrees of freedom
#' (i.e. \eqn{\nu S} follows a Wishart distribution with \eqn{\nu} degrees of
#' freedom and scale matrix \eqn{\Sigma}).
#'
#' \strong{The role of \code{nu} (read this before choosing a value).}
#' \code{nu} is the single most consequential setting in a covariance-based
#' fit, yet it behaves in a way that surprises many users.  Because the entire
#' negative log-likelihood above is multiplied by \eqn{\nu/2}, \code{nu}
#' enters as a global multiplier of the objective surface.  This has three
#' distinct consequences:
#' \itemize{
#'   \item \strong{Point estimates do not depend on \code{nu}.}  The maximizing
#'     conductance parameters \eqn{\theta} and nuisance parameters
#'     \eqn{(\tau, \sigma)} are completely invariant to \code{nu}: scaling the
#'     objective by a positive constant does not move its optimum.  You can
#'     change \code{nu} by orders of magnitude and the estimated surface is
#'     identical.
#'   \item \strong{Standard errors and confidence intervals scale as}
#'     \eqn{1/\sqrt{\nu}}.  The Hessian scales with \code{nu}, so the
#'     asymptotic covariance scales with \eqn{1/\nu}.  Doubling \code{nu}
#'     shrinks every standard error by a factor of \eqn{\sqrt{2}}.  This is the
#'     correct behavior: more markers carry more information and should yield
#'     tighter intervals.
#'   \item \strong{Model selection and likelihood-ratio tests depend strongly
#'     on \code{nu}.}  The log-likelihood scales linearly with \code{nu}, but
#'     the AIC/BIC complexity penalty (and the \eqn{\chi^2} reference
#'     distribution for an LRT) does not.  With large \code{nu} the penalty
#'     becomes negligible and extra covariates are almost always
#'     \dQuote{selected}; with small \code{nu} parsimony dominates.  The same
#'     data can therefore favor a simpler or a more complex model purely
#'     through the choice of \code{nu}.
#' }
#' In short, \code{nu} is an \emph{effective sample size}: it does not change
#' what the data say about the shape of the conductance surface, but it
#' determines how strongly the data speak.  Choosing \code{nu} requires care:
#' \itemize{
#'   \item \emph{Biallelic SNPs:} use the number of retained polymorphic SNPs.
#'     Reduce for linkage disequilibrium if markers are not independent.
#'   \item \emph{Microsatellites:} use the number of loci \eqn{L} as the
#'     conservative default.  Allele frequencies within a locus are correlated
#'     (they sum to a constant), so the locus is the natural unit of
#'     information.  \eqn{\sum_l (K_l - 1)}, where \eqn{K_l} is the number of
#'     observed alleles at locus \eqn{l}, is an upper bound that would only be
#'     correct if alleles within a locus were independent; they are not.
#'     Using it can produce confidence intervals that are too narrow and
#'     model-selection statistics that are too large.  The true effective
#'     \eqn{\nu} lies between \eqn{L} and \eqn{\sum_l (K_l - 1)}.  Report
#'     the value used and conduct a sensitivity analysis across that range.
#' }
#'
#' \strong{Typical workflow:}
#' \enumerate{
#'   \item Compute \code{S} from raw genotypes:
#'     \code{S <- cov_from_genetic_data(dosage_matrix, groups = pop_vector)}.
#'   \item Set \code{nu}: for biallelic SNPs use the retained SNP count; for
#'     microsatellites use the number of loci as the conservative starting point.
#'   \item Fit:
#'     \code{terradish(S ~ ..., measurement_model = wishart_covariance, nu = nu)}.
#' }
#'
#' If \code{S} has a non-positive-definite covariance (which can happen with
#' the \code{diagonal = "within"} option of \code{cov_from_genetic_data}),
#' inspect eigenvalues before fitting.
#'
#' @seealso \code{\link{cov_from_genetic_data}}, \code{\link{cov_from_biallelic}},
#'   \code{\link{generalized_wishart}}, \code{\link{wishart_covariates}},
#'   \code{\link{wishart_drift_covariates}},
#'   \code{\link{simulate_covariance_response}}, \code{\link{terradish}}
#'
#' @return When \code{phi} is missing, a list with elements \code{phi}
#'   (least-squares starting values for \code{(tau, sigma)}), \code{lower}
#'   (\code{c(0, -Inf)}), and \code{upper}.  Otherwise a list containing:
#'  \item{objective}{Negative log-likelihood.}
#'  \item{fitted}{Fitted covariance matrix \eqn{\Sigma} (same dimensions as \code{S}).}
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
class(wishart_covariance) <- c("terradish_measurement_model",
                               "radish_measurement_model")
