# Profile the two covariance components on site contrasts. Optimizing a log
# nugget directly can falsely converge after exp(sigma) becomes negligible.
# The overall variance has a closed-form optimum for every mixture weight.
# A deterministic search over that weight also removes warm-start history
# from the objective and its derivatives. This is a finite search, not a
# proof that an arbitrary multimodal variance profile has a unique optimum.
.wishart_profile_start <- function(E, S)
{
  L <- qr.Q(qr(stats::contr.helmert(nrow(E))))
  Ec <- crossprod(L, ((E + t(E)) / 2) %*% L)
  Sc <- crossprod(L, ((S + t(S)) / 2) %*% L)
  eig <- eigen(Ec, symmetric = TRUE)
  d <- eig$values
  scale <- mean(d)
  if (!is.finite(scale) || scale <= 0 ||
      min(d) < -sqrt(.Machine$double.eps) * max(abs(d)))
    return(NULL)
  d <- pmax(d, 0) / scale
  y <- diag(crossprod(eig$vectors, Sc %*% eig$vectors))
  if (any(!is.finite(y)) || any(y < 0) || mean(y) <= 0)
    return(NULL)
  m <- length(d)
  dv <- d - 1
  # Keep the limiting zero-nugget covariance representable in the existing
  # log-nugget interface. This is within roundoff of the endpoint, rather
  # than imposing a positive scientific nugget bound.
  upper <- 1 - .Machine$double.eps
  objective <- function(w)
  {
    v <- 1 + w * dv
    a <- mean(y / v)
    m * log(a) + sum(log(v)) + m
  }
  score <- function(w)
  {
    v <- 1 + w * dv
    sum(dv / v) - m * sum(y * dv / v^2) / sum(y / v)
  }
  grid <- sort(unique(c(seq(0, upper, length.out = 33L),
                        stats::plogis(seq(-36, 36, by = 0.5)), upper)))
  values <- vapply(grid, objective, numeric(1))
  slopes <- vapply(grid, score, numeric(1))
  candidates <- grid[c(1L, length(grid), which.min(values))]
  for (i in seq_len(length(grid) - 1L))
  {
    if (is.finite(slopes[i]) && is.finite(slopes[i + 1L]) &&
        slopes[i] < 0 && slopes[i + 1L] > 0)
      candidates <- c(candidates, stats::uniroot(score, grid[c(i, i + 1L)],
                                          tol = 1e-13)$root)
  }
  w <- candidates[which.min(vapply(candidates, objective, numeric(1)))]
  a <- mean(y / (1 + w * dv))
  c(tau = a * w / scale, sigma = log(a * (1 - w)))
}
