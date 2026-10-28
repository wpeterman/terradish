.safe_invert <- function(x)
{
  # A change of parameter units can make an otherwise regular Hessian look
  # singular. Equilibrate before either inversion or the rank-based fallback.
  scale <- sqrt(abs(diag(x)))
  scale[!is.finite(scale) | scale == 0] <- 1
  scaled <- x / outer(scale, scale)
  inverse <- tryCatch(solve(scaled), error = function(e) ginv(scaled))
  inverse / outer(scale, scale)
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

  lower_scale <- pmax(abs(phi), ifelse(is.finite(lower), abs(lower), 0))
  upper_scale <- pmax(abs(phi), ifelse(is.finite(upper), abs(upper), 0))
  at_lower <- is.finite(lower) & phi <= lower + sqrt(.Machine$double.eps) * lower_scale
  at_upper <- is.finite(upper) & phi >= upper - sqrt(.Machine$double.eps) * upper_scale
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

.retry_small_wishart_variance <- function(g, E, S, nu, answer, start, fit_once,
                                          nonnegative, control)
{
  supported <- identical(g, wishart_covariance) || identical(g, generalized_wishart) ||
    inherits(g, "terradish_wishart_covariate_model")
  index <- match("sigma", names(start$phi))
  if (is.na(index) && identical(g, generalized_wishart)) index <- 2L
  unchanged <- list(answer = answer, start = start, zero_limit = FALSE)
  if (!supported || is.na(index) || start$lower[index] == start$upper[index]) return(unchanged)
  parameters <- c(answer$par)
  distance_response <- identical(g, generalized_wishart) ||
    identical(attr(g, "base_model"), "generalized_wishart")
  response_covariance <- S
  if (distance_response) {
    centering <- diag(nrow(S)) - matrix(1 / nrow(S), nrow(S), nrow(S))
    response_covariance <- -.5 * centering %*% S %*% centering
  }
  reference <- max(mean(abs(diag(response_covariance))), .Machine$double.xmin)
  variance <- exp(parameters[index])
  small_variance <- function(value, convergence) is.finite(value) && value <= 1e-4 * reference
  if (!small_variance(variance, answer$convergence)) return(unchanged)
  evaluate <- function(values, gradient = FALSE) tryCatch(
    g(E = E, S = S, phi = values, nu = nu, gradient = gradient,
      hessian = FALSE, partial = FALSE, nonnegative = nonnegative), error = function(e) NULL)
  current <- evaluate(parameters)
  if (is.null(current) || !is.finite(current$objective)) return(unchanged)
  tolerance <- 10 * .Machine$double.eps * max(1, abs(current$objective))
  # Probe finite improving variances before accepting a tiny log-coordinate
  # gradient. Every retry retains the caller's full set of box constraints.
  seeds <- log(reference) + log(c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2))
  for (seed in seeds) {
    retry <- start; retry$phi <- parameters
    retry$phi[index] <- min(max(seed, start$lower[index]), start$upper[index])
    proposed <- evaluate(retry$phi)
    if (is.null(proposed) || proposed$objective >= current$objective - tolerance) next
    candidate <- tryCatch(fit_once(retry), error = function(e) NULL)
    if (is.null(candidate)) next
    evaluated <- evaluate(c(candidate$par))
    if (!is.null(evaluated) && is.finite(evaluated$objective) && evaluated$objective < current$objective - tolerance) {
      answer <- candidate; start <- retry; parameters <- c(candidate$par); current <- evaluated
    }
  }
  # A zero variance is an open endpoint on the public log scale. Represent its
  # supported limit with a variance below floating-point covariance precision,
  # and verify the one-sided variance gradient with the unchanged tolerance.
  zero_limit <- FALSE
  if (!is.finite(start$lower[index]) && small_variance(exp(parameters[index]), answer$convergence)) {
    floor_variance <- .Machine$double.xmin * max(1, reference)
    floor_log <- log(floor_variance)
    if (floor_log <= start$upper[index]) {
      retry <- start; retry$phi <- parameters; retry$phi[index] <- floor_log
      retry$lower[index] <- retry$upper[index] <- floor_log
      candidate <- tryCatch(fit_once(retry), error = function(e) NULL)
      if (!is.null(candidate) && isTRUE(candidate$convergence == 0)) {
        evaluated <- evaluate(c(candidate$par), gradient = TRUE)
        derivative <- if (is.null(evaluated)) NA_real_ else as.numeric(evaluated$gradient)[index] / floor_variance
        if (is.finite(derivative) && derivative >= -control$ctol &&
            is.finite(evaluated$objective) && evaluated$objective <= current$objective + tolerance) {
          answer <- candidate; start <- retry; zero_limit <- TRUE
        }
      }
    }
  }
  if (!zero_limit && exp(parameters[index]) <= sqrt(.Machine$double.eps) * reference &&
      !is.finite(start$lower[index])) {
    evaluated <- evaluate(parameters, gradient = TRUE)
    derivative <- if (is.null(evaluated)) NA_real_ else
      as.numeric(evaluated$gradient)[index] / exp(parameters[index])
    if ((!is.finite(derivative) || derivative < -control$ctol) && isTRUE(answer$convergence == 0L))
      answer$convergence <- 2L
  }
  list(answer = answer, start = start, zero_limit = zero_limit)
}

radish_subproblem <- function(g, E, S, nu, phi = NULL, nonnegative = TRUE, validate = FALSE, control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10, verbose = TRUE))
{
  phi_default <- g(E = E, S = S, nonnegative = nonnegative)
  if (is.null(phi))
    phi <- phi_default
  else if (!is.list(phi))
    phi <- list(phi = setNames(as.numeric(phi), names(phi_default$phi)), lower = phi_default$lower,
                 upper = phi_default$upper)

  # At an extreme logit, a negligible derivative can falsely pin rho when the
  # landscape changes. Reinitialize that coordinate, retaining the other warm
  # values. Equal bounds denote genuinely fixed predictive parameters.
  fixed <- !is.null(phi$lower) && !is.null(phi$upper) && all(phi$lower == phi$upper)
  parameter_scale <- rep(1, length(phi$phi))
  mlpe_bounds <- identical(as.numeric(phi$lower),
                           if (nonnegative) c(-Inf, 0, -Inf, -Inf) else rep(-Inf, 4)) &&
    identical(as.numeric(phi$upper), rep(Inf, 4))
  profiled_mlpe <- !fixed && mlpe_bounds && identical(g, mlpe)
  if (profiled_mlpe)
  {
    profile <- .mlpe_profile_start(E, S, nonnegative)
    phi$phi <- profile$phi
    parameter_scale <- profile$parameter_scale
    if (profile$upper_endpoint)
      stop("MLPE correlation profile reaches its upper numerical endpoint; inspect residual variation and nuisance identification.",
           call. = FALSE)
  }
  # The ordinary two-component Wishart models admit a deterministic scalar
  # variance profile. Use it for the standard bounds, including when a warm
  # log nugget has become spuriously stationary near minus infinity. Keep
  # fixed predictive parameters and user-supplied custom bounds unchanged.
  standard_bounds <- identical(as.numeric(phi$lower), c(0, -Inf)) &&
    identical(as.numeric(phi$upper), c(Inf, Inf))
  if (!fixed && standard_bounds &&
      (identical(g, wishart_covariance) || identical(g, generalized_wishart)))
  {
    profile_S <- if (identical(g, generalized_wishart))
      -0.5 * .prepare_gw_response(S) else S
    profile_phi <- .wishart_profile_start(E, profile_S)
    if (!is.null(profile_phi))
      phi$phi <- profile_phi
  }
  rho <- match("rho", names(phi$phi))
  if (!fixed && !profiled_mlpe && !is.na(rho) && phi$lower[rho] != phi$upper[rho] && abs(phi$phi[rho]) > 8)
    phi$phi[rho] <- min(max(phi_default$phi[rho], phi$lower[rho]), phi$upper[rho])

  # use Newton-Raphson to profile out nuisance parameters
  fit_subproblem <- function(phi_start)
  {
    if (all(phi_start$lower == phi_start$upper)) {
      evaluation <- g(E = E, S = S, nu = nu, phi = c(phi_start$phi),
                      gradient = FALSE, hessian = FALSE, partial = FALSE,
                      nonnegative = nonnegative)
      if (!is.finite(evaluation$objective)) stop("Fixed measurement parameters give a non-finite objective.")
      return(list(par = matrix(phi_start$phi, ncol = 1), convergence = 0L, iters = 0L))
    }
    result <- BoxConstrainedNewton(phi_start$phi / parameter_scale,
                         function(par, gradient, hessian) {
                           out <- g(E = E, S = S, nu = nu, phi = c(par) * parameter_scale,
                             gradient = gradient,
                             hessian = hessian,
                             partial = FALSE,
                             nonnegative = nonnegative)
                           if (gradient) out$gradient <- out$gradient * parameter_scale
                           if (hessian) out$hessian <- out$hessian * outer(parameter_scale, parameter_scale)
                           out
                         },
                         lower = phi_start$lower / parameter_scale,
                         upper = phi_start$upper / parameter_scale,
                         control = control)
    result$par <- result$par * parameter_scale
    result
  }
  phi_start  <- phi
  subproblem <- tryCatch(fit_subproblem(phi_start),
                         error = function(e) {
                           if (identical(phi_start, phi_default))
                             stop(e)
                           if (isTRUE(control$verbose))
                             message("Nuisance warm start failed; retrying the default start: ",
                                     conditionMessage(e))
                           retry <- phi_start
                           retry$phi <- pmin(pmax(phi_default$phi, retry$lower), retry$upper)
                           phi_start <<- retry
                           fit_subproblem(retry)
                         })

  # A warm start can return a finite but unconverged nuisance fit rather than
  # throw an error. Retry the default start within the caller's bounds before
  # passing that failed profile to the outer optimizer. Keep failure visible
  # if the retry also fails; a second attempt is not evidence of convergence.
  if (!isTRUE(subproblem$convergence == 0L)) {
    retry <- phi_start
    retry$phi <- pmin(pmax(phi_default$phi, retry$lower), retry$upper)
    if (!identical(as.numeric(retry$phi), as.numeric(phi_start$phi))) {
      candidate <- tryCatch(fit_subproblem(retry), error = function(e) NULL)
      if (!is.null(candidate) && isTRUE(candidate$convergence == 0L)) {
        subproblem <- candidate
        phi_start <- retry
      }
    }
  }

  recovered <- .retry_small_wishart_variance(g, E, S, nu, subproblem, phi_start,
                                             fit_subproblem, nonnegative, control)
  subproblem <- recovered$answer
  phi_start <- recovered$start

  # refit, computing partial derivatives
  phi         <- subproblem$par
  rownames(phi) <- names(phi_default$phi)
  fit         <- g(E = E, S = S, nu = nu, phi = c(phi), partial = TRUE, nonnegative = nonnegative)
  fit$no_structure_boundary <- if (is.null(fit$no_structure_boundary)) fit$boundary else fit$no_structure_boundary
  fit$residual_variance_zero_limit <- recovered$zero_limit
  fit$boundary <- isTRUE(fit$boundary) || recovered$zero_limit
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
  jacobian_phi <- function(dotdotE)
    -invhess %*% crossprod(partial_E, c(dotdotE))
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
       no_structure_boundary = if (is.null(fit$no_structure_boundary))
         fit$boundary else fit$no_structure_boundary,
       phi            = phi,
       convergence    = subproblem$convergence,
       iters          = subproblem$iters,
       gradient       = gradient_E,
       jacobian_E     = jacobian_E,
       jacobian_phi   = jacobian_phi,
       jacobian_S     = jacobian_S,
       num_jacobian_E = if(!validate) NULL else num_jacobian_E,
       num_jacobian_S = if(!validate) NULL else num_jacobian_S)
}
