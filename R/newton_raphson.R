#' Control settings for Newton-like optimizers
#'
#' Tuning parameters shared by Newton-Raphson and BFGS quasi-Newton
#' optimizers inside \code{\link{terradish}}.
#'
#' @param maxit Maximum number of Newton or quasi-Newton steps.  Increase if
#'   the optimizer reports a convergence warning on large or difficult problems.
#' @param ctol Gradient convergence tolerance.  Optimization stops when
#'   the largest absolute projected gradient is below \code{ctol}.
#' @param ftol Objective convergence tolerance.  Optimization also stops when
#'   the absolute change in the objective is below \code{ftol} and the largest
#'   absolute projected gradient is below \code{sqrt(ctol)}.
#' @param etol Eigenvalue threshold used to detect near-singular Hessians.
#'   Absolute eigenvalues smaller than \code{etol} times the largest absolute
#'   Hessian entry are replaced by one when forming the Newton step.
#'   Increasing this value regularizes the step
#'   at the cost of slower convergence near flat regions.
#' @param verbose Logical.  If \code{TRUE}, print the iteration count,
#'   objective value, and gradient norm at each step.
#' @param eps Box-constraint tolerance.  A parameter is considered to be
#'   "on the boundary" when it is within \code{eps * abs(par)} of the bound.
#'   Active-constraint gradient components are zeroed before the step is taken.
#' @param del Typical initial step length for quasi-Newton (BFGS) methods.
#'   Used to scale the initial approximate Hessian before any updates are
#'   available.
#' @param ls.control Control object for the line search, created by
#'   \code{\link{HagerZhangControl}} (default) or \code{\link{ArmijoControl}}.
#'   \code{ArmijoControl} is only valid when \code{optimizer = "bfgs"} is
#'   passed to \code{\link{terradish}}.
#'
#' @details
#' The projected gradient sets outward-pointing components at active bounds
#' to zero. A sufficiently small projected gradient gives convergence code 0.
#' A small objective change also gives code 0 if the projected gradient is
#' below \code{sqrt(ctol)}. Otherwise the optimizer warns that it has stalled
#' and returns code 2. Reaching the iteration limit gives code 1.
#'
#' For most resistance-surface models the defaults converge reliably.  Consider
#' changing them when:
#'
#' \itemize{
#'   \item The optimizer hits \code{maxit} before converging: increase
#'     \code{maxit}.
#'   \item You want a quick exploratory fit: reduce \code{maxit} (for example
#'     5 to 10).
#'   \item Fitting is slow because of oscillation: tighten \code{ctol} or
#'     switch to a more conservative line search via \code{ls.control}.
#' }
#'
#' @return A named list of class \code{"NewtonRaphsonControl"} holding the
#'   settings passed to the optimizer.  The elements are the arguments above,
#'   validated and stored unchanged: \code{maxit}, \code{ctol}, \code{ftol},
#'   \code{etol}, \code{verbose}, \code{eps}, \code{del}, and
#'   \code{ls.control} (itself a \code{HagerZhangControl} or
#'   \code{ArmijoControl} list).  Pass the object to the \code{control}
#'   argument of \code{\link{terradish}}; inspect it with \code{str()}.  The
#'   list carries no fitted quantities, so nothing in it needs
#'   back-transforming.
#'
#' @seealso \code{\link{HagerZhangControl}}, \code{\link{ArmijoControl}},
#'   \code{\link{terradish}}
#'
#' @examples
#' ctrl <- NewtonRaphsonControl(maxit = 25, verbose = FALSE)
#' str(ctrl)
#'
#' # Relaxed tolerances for a quick exploratory run:
#' ctrl_quick <- NewtonRaphsonControl(maxit = 10, ctol = 1e-3, ftol = 1e-3)
#'
#' \donttest{
#' data(melip)
#' covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#' covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
#'                                                 na.rm = TRUE))
#' surface <- conductance_surface(covariates, terra::unwrap(melip.coords),
#'                                directions = 8)
#'
#' # A tighter, longer-running control than the default
#' fit <- terradish(melip.Fst ~ forestcover + altitude, data = surface,
#'                  conductance_model = loglinear_conductance,
#'                  measurement_model = mlpe,
#'                  control = NewtonRaphsonControl(maxit = 200, ctol = 1e-8))
#' coef(fit)
#'
#' # The same fit, reporting each optimizer step. `verbose` is the simple
#' # switch; it overrides whatever `control` says.
#' fit_loud <- terradish(melip.Fst ~ forestcover + altitude, data = surface,
#'                       conductance_model = loglinear_conductance,
#'                       measurement_model = mlpe,
#'                       control = ctrl_quick, verbose = TRUE)
#' }
#'
#' @export
NewtonRaphsonControl <- function(maxit = 100, 
                                 ctol = sqrt(.Machine$double.eps), 
                                 ftol = sqrt(.Machine$double.eps), 
                                 etol = 10*.Machine$double.eps, 
                                 verbose = FALSE, 
                                 eps = 1e-8, 
                                 del = 1,
                                 ls.control = HagerZhangControl())
  list(maxit = maxit, ctol = ctol, etol = etol, 
       ftol = ftol, verbose = verbose, eps = eps, 
       del = del, ls.control = ls.control)

# Both optimizers iterate `seq_len(maxit)` and then read the loop variable, so
# a non-positive or non-finite `maxit` has to be rejected up front rather than
# silently producing a zero-length or reversed sequence.
.terradish_validate_maxit <- function(maxit)
{
  maxit <- suppressWarnings(as.integer(maxit)[1])
  if (is.na(maxit) || maxit < 1L)
    stop("`maxit` must be a single integer of at least 1.", call. = FALSE)
  maxit
}

BoxConstrainedNewton <- function(par, fn, lower = rep(-Inf, length(par)), upper = rep(Inf, length(par)), control = NewtonRaphsonControl())
{
  BoxConstrainedNewtonNaN <- function()
  {
    list(objective = NaN,
         gradient  = matrix(NaN, length(par), 1),
         hessian   = matrix(NaN, length(par), length(par)))
  }
  
  prettify <- function(x)
    formatC(x, digits=3, width=5, format="e")

  zero_bounded_variables <- function(gradient, par, lower, upper, eps = 1e-8)
  {
    # set gradient to 0 for active constraints
    tol <- eps * abs(par)
    gradient <- ifelse(upper - tol <= par & gradient < 0, 0, gradient)
    gradient <- ifelse(lower + tol >= par & gradient > 0, 0, gradient)
    gradient
  }

  gap_step_bounded_variables <- function(desc, par, gradient, lower, upper, eps = 1e-8)
  {
    # modify search direction so that at alpha == 1, actively constrained variables are set to the boundary
    tol <- eps * abs(par)
    desc <- ifelse(upper - tol <= par & gradient < 0, upper - par, desc)
    desc <- ifelse(lower + tol >= par & gradient > 0, lower - par, desc)
    desc
  }

  project <- function(x, lower, upper)
    pmin(pmax(x, lower), upper)

  stopifnot(lower <= upper)

  maxit <- control$maxit
  ctol <- control$ctol
  etol <- control$etol
  ftol <- control$ftol
  verbose <- control$verbose
  eps <- control$eps
  del <- control$del
  ls.control <- control$ls.control
  etol <- etol * length(par)
  if (.terradish_is_armijo_control(ls.control))
    stop("`ArmijoControl()` is currently supported for BFGS optimization only; use `optimizer = \"bfgs\"` or the default Hager-Zhang line search for Newton.",
         call. = FALSE)

  if (verbose)
    message("Projected Newton-Raphson with Hager-Zhang line search")

  # `maxit` must be at least one step: the loop below defines `fit` and `i`,
  # and the convergence check after it reads both.
  maxit <- .terradish_validate_maxit(maxit)

  convergence <- 0
  criterion <- "iteration_limit"
  line_search_failed <- FALSE
  par <- as.matrix(par)

  for (i in seq_len(maxit))
  {
    fit   <- fn(par, gradient = TRUE, hessian = TRUE)
    delta <- if (i > 1) abs(oldfit$objective - fit$objective) else 0

    if (verbose)
      message(paste0("[", i, "]"),
              " f(x) = ", prettify(-fit$objective),
              "  |f(x)-fold(x)| = ", prettify(delta),
              "  max|f'(x)| = ", prettify(max(abs(fit$gradient))),
              "  |f''(x)| = ", prettify(-det(fit$hessian)))

    gradient     <- fit$gradient
    gradient_box <- zero_bounded_variables(gradient, par, lower, upper, eps)
    projected_norm <- max(abs(gradient_box))
    if (projected_norm < ctol) {
      criterion <- "projected_gradient"
      break
    }
    if (i > 1 && delta < ftol) {
      if (projected_norm < sqrt(ctol)) {
        criterion <- "objective_and_projected_gradient"
      } else {
        convergence <- 2L
        criterion <- "stalled"
        warning("Optimizer stalled: objective change is small but the projected gradient remains large.",
                call. = FALSE)
      }
      break
    }
    # Numerical Hessians can drift slightly away from exact symmetry on larger
    # problems. Symmetrizing keeps the Newton step real-valued and avoids
    # complex eigendecompositions in the line search.
    # A parameter pinned at a bound cannot move to offset a free-coordinate
    # step. Invert only the free Hessian block, rather than using the free
    # rows of the unconstrained inverse (a different quadratic problem).
    bound_tol <- eps * abs(par)
    active <- (upper - bound_tol <= par & gradient < 0) |
      (lower + bound_tol >= par & gradient > 0) | lower == upper
    free <- which(!active)
    hessian_sym <- fit$hessian[free, free, drop = FALSE]
    hessian_sym <- (hessian_sym + t(hessian_sym)) / 2
    ehess        <- eigen(hessian_sym, symmetric = TRUE)
    ehess$values <- abs(ehess$values)
    ehess$values <- ifelse(ehess$values < max(abs(hessian_sym)) * etol, 1, ehess$values)
    ihess <- ehess$vectors %*% solve(diag(ehess$values, nrow=length(free))) %*% t(ehess$vectors)
    desc <- matrix(0, length(par), 1L)
    desc[free, ] <- -ihess %*% gradient_box[free, , drop = FALSE]
    desc <- gap_step_bounded_variables(desc, par, gradient, lower, upper, eps)
    phi0         <- fit$objective
    dphi0        <- c(t(desc) %*% gradient_box)

    line_cache <- new.env(parent = emptyenv())
    dphi_fn <- function(alpha) 
    {
      cache_key <- .terradish_line_search_cache_key(alpha)
      cached <- line_cache[[cache_key]]
      if (!is.null(cached))
        return(cached)

      .terradish_record_line_search_trial(control, gradient = TRUE)
      tryCatch({
        phi <- fn(project(par + alpha*desc, lower, upper), gradient = TRUE, hessian = FALSE)
        grb <- zero_bounded_variables(phi$gradient, par + alpha*desc, lower, upper, eps)
        value <- list(objective = phi$objective, gradient = c(t(desc) %*% grb))
        line_cache[[cache_key]] <- value
        value
      }, error = function(e) {
        BoxConstrainedNewtonNaN()
      })
    }

    #alpha <- HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    alpha <- tryCatch({
      HagerZhang(dphi_fn, phi0, dphi0, control = ls.control)
    }, error = function(err) {
      if (verbose)
        message("Hager-Zhang line search failed; switching to bounded backtracking.")
      Backtracking(dphi_fn, phi0, dphi0, control = ls.control)
    })
    if (!is.finite(alpha) || alpha <= 0 ||
        identical(attr(alpha, "line_search_status"), "failed"))
    {
      convergence <- 2
      line_search_failed <- TRUE
      criterion <- "line_search_failed"
      warning("Failed to find a usable line-search step; returning the current parameter values.",
              call. = FALSE, immediate. = TRUE)
      break
    }
    par <- project(par + alpha*desc, lower, upper)

    oldfit <- fit
  }

  boundary_fit <- any(par == lower | par == upper)
  if (verbose)
    message("Solution on ",
            if (boundary_fit) "boundary" else "interior",
            " with `max(abs(gradient))` == ", max(abs(fit$gradient)),
            " and `diff(f)` == ", delta)

  if (identical(criterion, "iteration_limit"))
  {
    fit <- fn(par, gradient = TRUE, hessian = TRUE)
    gradient_box <- zero_bounded_variables(fit$gradient, par, lower, upper, eps)
    warning("`maxit` reached for Newton steps", immediate. = TRUE)
    convergence = 1
  } 

  list(par = par,
       gradient = fit$gradient,
       hessian = fit$hessian,
       value = fit$objective,
       fit = fit,
       iters = i,
       boundary = boundary_fit,
       criterion = criterion,
       max_abs_projected_gradient = max(abs(gradient_box)),
       convergence = convergence)
}
