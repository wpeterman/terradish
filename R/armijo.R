#' Control settings for Armijo backtracking line search
#'
#' Tuning parameters for an objective-only Armijo line search. This line search
#' is currently intended for \code{optimizer = "bfgs"} fits where trial
#' objective evaluations are cheaper than repeated full objective-plus-gradient
#' evaluations during Hager-Zhang line search.
#'
#' @param initial Initial trial step length.
#' @param contraction Multiplicative shrinkage applied after rejected trial
#'   steps. Must be strictly between 0 and 1.
#' @param sufficient_decrease Armijo sufficient-decrease constant.
#' @param maxit Maximum number of backtracking steps.
#' @param min_alpha Smallest step length to try before stopping.
#' @param verbose Should line-search progress be printed?
#'
#' @details
#' \code{ArmijoControl()} is opt-in. The default optimizer control continues to
#' use \code{\link{HagerZhangControl}}. Armijo trial steps only evaluate the
#' objective, so rejected steps avoid the gradient/backpropagation work required
#' by Hager-Zhang. After a step is accepted, BFGS still evaluates the full
#' gradient once at the accepted point so the quasi-Newton update remains
#' well-defined.
#'
#' @examples
#' ctrl <- NewtonRaphsonControl(
#'   maxit = 25,
#'   verbose = FALSE,
#'   ls.control = ArmijoControl()
#' )
#' str(ctrl$ls.control)
#'
#' @export
ArmijoControl <- function(initial = 1.0,
                          contraction = 0.5,
                          sufficient_decrease = 1e-4,
                          maxit = 25L,
                          min_alpha = sqrt(.Machine$double.eps),
                          verbose = FALSE)
{
  if (!is.finite(initial) || initial <= 0)
    stop("`initial` must be a positive finite step length", call. = FALSE)
  if (!is.finite(contraction) || contraction <= 0 || contraction >= 1)
    stop("`contraction` must be strictly between 0 and 1", call. = FALSE)
  if (!is.finite(sufficient_decrease) ||
      sufficient_decrease <= 0 ||
      sufficient_decrease >= 1)
    stop("`sufficient_decrease` must be strictly between 0 and 1", call. = FALSE)

  maxit <- as.integer(maxit)[1]
  if (is.na(maxit) || maxit < 1L)
    stop("`maxit` must be a positive integer", call. = FALSE)
  if (!is.finite(min_alpha) || min_alpha < 0)
    stop("`min_alpha` must be a non-negative finite value", call. = FALSE)

  structure(
    list(
      type = "armijo",
      initial = initial,
      contraction = contraction,
      sufficient_decrease = sufficient_decrease,
      maxit = maxit,
      min_alpha = min_alpha,
      verbose = isTRUE(verbose)
    ),
    class = c("terradish_armijo_control", "list")
  )
}

.terradish_is_armijo_control <- function(control)
{
  inherits(control, "terradish_armijo_control") ||
    identical(control$type, "armijo")
}

Armijo <- function(phifn, phi_0, dphi_0, control = ArmijoControl())
{
  if (!is.finite(phi_0) || !is.finite(dphi_0))
    stop("Value and slope at step length = 0 must be finite.", call. = FALSE)
  if (dphi_0 >= 0)
    return(0)

  alpha <- control$initial
  best_alpha <- 0
  best_objective <- phi_0

  for (iter in seq_len(control$maxit))
  {
    fit <- phifn(alpha)
    phi_alpha <- if (is.list(fit)) fit$objective else fit
    phi_alpha <- c(phi_alpha)[1]

    if (is.finite(phi_alpha))
    {
      if (phi_alpha < best_objective)
      {
        best_objective <- phi_alpha
        best_alpha <- alpha
      }

      if (phi_alpha <= phi_0 + control$sufficient_decrease * alpha * dphi_0)
        return(alpha)
    }

    if (isTRUE(control$verbose))
      cat("Armijo backtracking: alpha =", alpha, "objective =", phi_alpha, "\n")

    alpha <- alpha * control$contraction
    if (alpha < control$min_alpha)
      break
  }

  warning("Armijo line search did not satisfy sufficient decrease; using best finite trial step",
          call. = FALSE)
  best_alpha
}
