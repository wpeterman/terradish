.safe_invert <- function(x)
{
  tryCatch(solve(x), error = function(e) ginv(x))
}

# Which nuisance parameters are free to move at the profiled optimum? A
# parameter sitting on an active box constraint is pinned: perturbing E (or S)
# does not move it, so its implicit derivative is zero.
.free_parameter_index <- function(phi, lower = NULL, upper = NULL)
{
  phi <- as.numeric(phi)
  p <- length(phi)
  if (is.null(lower))
    lower <- rep(-Inf, p)
  if (is.null(upper))
    upper <- rep(Inf, p)
  lower <- rep_len(as.numeric(lower), p)
  upper <- rep_len(as.numeric(upper), p)

  tol <- sqrt(.Machine$double.eps) * pmax(1, abs(phi))
  at_lower <- is.finite(lower) & phi <= lower + tol
  at_upper <- is.finite(upper) & phi >= upper - tol
  !(at_lower | at_upper)
}

# Invert the free-parameter block of the nuisance Hessian and pad the pinned
# rows and columns with zeros, so the profile-likelihood correction contributes
# nothing along directions the optimizer cannot move.
.constrained_inverse_hessian <- function(hessian, free)
{
  hessian <- .as_base_matrix(hessian)
  p <- nrow(hessian)
  out <- matrix(0, p, p, dimnames = dimnames(hessian))
  if (any(free))
    out[free, free] <- .safe_invert(hessian[free, free, drop = FALSE])
  out
}

.as_base_matrix <- function(x)
{
  if (is.null(x))
    return(NULL)
  as.matrix(x)
}

.as_base_vector <- function(x)
{
  as.vector(.as_base_matrix(x))
}

.response_partial_format <- function(partial_S, S)
{
  partial_S <- .as_base_matrix(partial_S)
  n_response <- ncol(partial_S)
  n_lower <- sum(lower.tri(S))
  n_lower_diag <- sum(lower.tri(S, diag = TRUE))
  n_full <- length(S)

  if (n_response == n_lower)
    return("lower")
  if (n_response == n_lower_diag)
    return("lower_diag")
  if (n_response == n_full)
    return("full")

  stop("Unsupported response-derivative dimension in `measurement_model`.", call. = FALSE)
}

.expand_response_partials <- function(x, S, format)
{
  x <- .as_base_vector(x)
  out <- matrix(0, nrow(S), ncol(S))

  if (identical(format, "lower"))
  {
    out[lower.tri(out)] <- x
    return(out + t(out))
  }

  if (identical(format, "lower_diag"))
  {
    out[lower.tri(out, diag = TRUE)] <- x
    diag_vals <- diag(out)
    out <- out + t(out)
    diag(out) <- diag_vals
    return(out)
  }

  if (identical(format, "full"))
    return(matrix(x, nrow(S), ncol(S)))

  stop("Unknown response-derivative format.", call. = FALSE)
}

radish_subproblem <- function(g, E, S, nu, phi = NULL, nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10, verbose = TRUE))
{
  phi_default <- g(E = E, S = S, nonnegative = nonnegative)
  if (is.null(phi))
    phi <- phi_default

  # use Newton-Raphson to profile out nuisance parameters
  fit_subproblem <- function(phi_start)
  {
    BoxConstrainedNewton(phi_start$phi,
                         function(par, gradient, hessian)
                           g(E = E, S = S, nu = nu, phi = c(par),
                             gradient = gradient,
                             hessian = hessian,
                             partial = FALSE,
                             nonnegative = nonnegative),
                         lower = phi_start$lower,
                         upper = phi_start$upper,
                         control = control)
  }
  phi_start  <- phi
  subproblem <- tryCatch(fit_subproblem(phi_start),
                         error = function(e) {
                           if (identical(phi_start, phi_default))
                             stop(e)
                           phi_start <<- phi_default
                           fit_subproblem(phi_default)
                         })

  # refit, computing partial derivatives
  phi         <- subproblem$par
  fit         <- g(E = E, S = S, nu = nu, phi = c(phi), partial = TRUE, nonnegative = nonnegative)
  gradient_E  <- fit$gradient_E

  # for hessian, need to get d(dg/dE)/dE via adjoint method,
  #   dg/dE = \partial (dg/dE)/\partial E + \partial (dg/dE)/\partial \hat{phi} \times \partial \hat{\phi}/\partial E
  # where
  #   dg/dphi = 0 ==> d(dg/dphi)/dE = 0 ==>
  #       \partial (dg/dphi)/\partial phi \times dphi/dE = 0 ==>
  #       dphi/dE = -[\partial (dg/dphi)/\partial phi]^-1 \partial (dg/dphi)/\partial E
  #
  # That stationarity argument only holds for nuisance parameters that are free
  # to move. A parameter pinned at an active box constraint (for example a
  # non-negative kernel coefficient estimated at exactly 0) has dphi/dE = 0, so
  # it must be dropped from the correction; including it inflates the implied
  # curvature and biases the standard errors. `.free_parameter_index()` finds
  # the inactive set and `.constrained_inverse_hessian()` inverts only that
  # block, padding the rest with zeros.
  partial_E   <- .as_base_matrix(fit$partial_E)
  bounds      <- if (is.list(phi_start)) phi_start else phi_default
  free_phi    <- .free_parameter_index(phi, bounds$lower, bounds$upper)
  invhess     <- .constrained_inverse_hessian(fit$hessian, free_phi)
  jacobian_E  <- function(dotdotE)
  {
    dotdotE_matrix <- .as_base_matrix(dotdotE)
    dphi_dE        <- -matrix(.as_base_vector(dotdotE_matrix) %*% partial_E %*% invhess %*% t(partial_E),
                              nrow(dotdotE_matrix), ncol(dotdotE_matrix))
    return (fit$jacobian_E(dotdotE_matrix) + dphi_dE)
  }

  # likewise, to get the leverage dtheta/dy, 
  #    dl/dtheta = 0 ==> d(dl/dtheta)/dy = 0 ==>
  #       \partial (dl/dtheta)/\partial y + \partial (dl/dtheta)/\partial theta \times dtheta/dy = 0 ==>
  #       dtheta/dy = -[\partial (dl/dtheta)/\partial theta]^{-1} \partial (dl/dtheta)/partial y
  # so, need the change in the gradient with y
  partial_S   <- .as_base_matrix(fit$partial_S)
  partial_S_format <- .response_partial_format(partial_S, S)
  jacobian_S  <- function(dotdotE)
  {
    dotdotE_matrix <- .as_base_matrix(dotdotE)
    dphi_dS <- .expand_response_partials(
      -.as_base_vector(dotdotE_matrix) %*% partial_E %*% invhess %*% partial_S,
      S = S,
      format = partial_S_format
    )
    return (fit$jacobian_S(dotdotE_matrix) + dphi_dS)
  }

  # numerical validation
  if (validate)
  {
    silence <- function(control) { control$verbose = FALSE; control }

    num_jacobian_E <- function(X) 
      matrix(c(X) %*% .numderiv_jacobian(function(x) 
                                         radish_subproblem(g = g, 
                                                           E = x, 
                                                           S = S, 
                                                           phi = phi, 
                                                           nonnegative = nonnegative,
                                                           control = silence(control))$gradient, 
                                         E), 
             nrow(E), ncol(E))

    num_jacobian_S <- function(X) 
      matrix(c(X) %*% .numderiv_jacobian(function(x) 
                                         radish_subproblem(g = g, 
                                                           E = E, 
                                                           S = x, 
                                                           phi = phi, 
                                                           control = silence(control))$gradient, 
                                         S), 
             nrow(S), ncol(S))
  }

  list(fit            = fit,
       loglikelihood  = fit$objective,
       boundary       = fit$boundary,
       phi            = phi,
       gradient       = gradient_E,
       jacobian_E     = jacobian_E,
       jacobian_S     = jacobian_S,
       num_jacobian_E = if(!validate) NULL else num_jacobian_E,
       num_jacobian_S = if(!validate) NULL else num_jacobian_S)
}
