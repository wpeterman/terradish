# simple bounded backtracking to find a finite improving function value
Backtracking <- function(dphifn, phi_0, dphi_0, control = NULL)
{
  line_search_step <- function(alpha, status, iter, objective = NA_real_)
  {
    alpha <- as.numeric(alpha)[1]
    attr(alpha, "line_search_status") <- status
    attr(alpha, "line_search_iterations") <- iter
    attr(alpha, "line_search_objective") <- objective
    alpha
  }

  maxit <- if (!is.null(control$linesearchmax))
    min(as.integer(control$linesearchmax)[1], 12L)
  else
    10L
  if (is.na(maxit) || maxit < 1L)
    maxit <- 10L

  rho <- 1e-4
  contraction <- if (!is.null(control$psi3) &&
                     is.finite(control$psi3) &&
                     control$psi3 > 0 &&
                     control$psi3 < 1)
    control$psi3
  else
    0.2
  alpha <- if (!is.null(control$c) && is.finite(control$c) && control$c > 0)
    control$c
  else
    1.0
  if (!is.null(control$alphamax) && is.finite(control$alphamax))
    alpha <- min(alpha, control$alphamax)

  best_alpha <- 0
  best_objective <- phi_0
  for (iter in seq_len(maxit))
  {
    if (!is.finite(alpha) || alpha <= .Machine$double.eps)
      break
    ev <- dphifn(alpha)
    val <- ev$objective
    gra <- ev$gradient

    if (is.finite(val) && is.finite(gra))
    {
      if (val < best_objective)
      {
        best_alpha <- alpha
        best_objective <- val
      }
      if (is.finite(dphi_0) &&
          dphi_0 < 0 &&
          val <= phi_0 + alpha * rho * dphi_0)
        return(line_search_step(alpha, "accepted", iter, val))
    }
    alpha <- alpha * contraction
  }

  if (best_alpha > 0)
  {
    warning("Backtracking did not satisfy sufficient decrease; using the best finite improving step.",
            call. = FALSE, immediate. = TRUE)
    return(line_search_step(best_alpha, "best", maxit, best_objective))
  }

  warning("Backtracking failed to find a finite improving step; returning alpha = 0.",
          call. = FALSE, immediate. = TRUE)
  line_search_step(0, "failed", maxit)
}
