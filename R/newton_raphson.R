#' Control settings for Newton-like optimizers
#'
#' Tuning parameters shared by Newton-Raphson and BFGS quasi-Newton
#' optimizers inside \code{\link{terradish}}.
#'
#' @param maxit Maximum number of Newton or quasi-Newton steps.  Increase if
#'   the optimizer reports a convergence warning on large or difficult problems.
#' @param ctol Gradient convergence tolerance.  Optimization stops when
#'   \code{max(abs(gradient)) < ctol}.
#' @param ftol Objective convergence tolerance.  Optimization also stops when
#'   the improvement in the objective over one step falls below \code{ftol}.
#' @param etol Eigenvalue threshold used to detect near-singular Hessians.
#'   Eigenvalues smaller than \code{etol * length(par)} are treated as zero
#'   when forming the Newton step.  Increasing this value regularizes the step
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
#' Both convergence criteria (\code{ctol} and \code{ftol}) must be satisfied
#' simultaneously before the optimizer declares convergence.  If you want
#' gradient-only convergence, set \code{ftol} to a very small number; for
#' objective-only, set \code{ctol} small.
#'
#' For most resistance-surface models the defaults converge reliably.  Consider
#' changing them when:
#'
#' \itemize{
#'   \item The optimizer hits \code{maxit} before converging — increase
#'     \code{maxit}.
#'   \item You want a quick exploratory fit — reduce \code{maxit} (e.g. 5–10).
#'   \item Fitting is slow because of oscillation — tighten \code{ctol} or
#'     switch to a more conservative line search via \code{ls.control}.
#' }
#'
#' @seealso \code{\link{HagerZhangControl}}, \code{\link{ArmijoControl}},
#'   \code{\link{terradish}}
#'
#' @examples
#' ctrl <- NewtonRaphsonControl(maxit = 25, verbose = FALSE)
#' str(ctrl)
#'
#' # Verbose fit with relaxed tolerances for a quick exploratory run:
#' ctrl_quick <- NewtonRaphsonControl(maxit = 10, verbose = TRUE,
#'                                    ctol = 1e-3, ftol = 1e-3)
#'
#' # Pass to terradish():
#' \dontrun{
#' fit <- terradish(melip.Fst ~ forestcover + altitude, data = surface,
#'                  conductance_model = loglinear_conductance,
#'                  measurement_model = mlpe,
#'                  control = NewtonRaphsonControl(maxit = 200))
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

  stopifnot(lower < upper)

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
    cat("Projected Newton-Raphson with Hager-Zhang line search\n")

  convergence <- 0
  line_search_failed <- FALSE
  par <- as.matrix(par)

  for (i in 1:maxit)
  {
    fit   <- fn(par, gradient = TRUE, hessian = TRUE)
    delta <- if (i > 1) abs(oldfit$objective - fit$objective) else 0

    if (verbose)
      cat(paste0("[", i, "]"), 
          "f(x) =", prettify(-fit$objective),
          "  |f(x)-fold(x)| =", prettify(delta),
          "  max|f'(x)| =", prettify(max(abs(fit$gradient))),
          "  |f''(x)| =", prettify(-det(fit$hessian)),
          "\n")

    if (max(abs(fit$gradient)) < ctol || (i > 1 && delta < ftol))
      break

    gradient     <- fit$gradient
    gradient_box <- zero_bounded_variables(gradient, par, lower, upper, eps)
    # Numerical Hessians can drift slightly away from exact symmetry on larger
    # problems. Symmetrizing keeps the Newton step real-valued and avoids
    # complex eigendecompositions in the line search.
    hessian_sym  <- (fit$hessian + t(fit$hessian)) / 2
    ehess        <- eigen(hessian_sym, symmetric = TRUE)
    ehess$values <- abs(ehess$values)
    ehess$values <- ifelse(ehess$values < max(abs(hessian_sym)) * etol, 1, ehess$values)
    ihess        <- ehess$vectors %*% solve(diag(ehess$values, nrow=length(par))) %*% t(ehess$vectors)
    desc         <- gap_step_bounded_variables(-ihess %*% gradient_box, par, gradient, lower, upper, eps)
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
      message("Hager-Zhang line search failed; switching to bounded backtracking.")
      Backtracking(dphi_fn, phi0, dphi0, control = ls.control)
    })
    if (!is.finite(alpha) || alpha <= 0 ||
        identical(attr(alpha, "line_search_status"), "failed"))
    {
      convergence <- 2
      line_search_failed <- TRUE
      warning("Failed to find a usable line-search step; returning the current parameter values.",
              call. = FALSE, immediate. = TRUE)
      break
    }
    par <- project(par + alpha*desc, lower, upper)

    oldfit <- fit
  }

  boundary_fit <- any(par == lower | par == upper)
  if (verbose)
    if (boundary_fit)
      cat ("Solution on boundary with `max(abs(gradient))` ==", max(abs(fit$gradient)), "and `diff(f)` ==", delta, "\n")
    else
      cat ("Solution on interior with `max(abs(gradient))` ==", max(abs(fit$gradient)), "and `diff(f)` ==", delta, "\n")

  if (!line_search_failed && i == maxit)
  {
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
       convergence = convergence)
}
