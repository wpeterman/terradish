# Center and scale the resistance predictor before regression. The MLPE
# likelihood must not change when resistance is expressed in different units.
.mlpe_pair_data <- function(E, S)
{
  R <- outer(diag(E), diag(E), "+") - E - t(E)
  x <- R[lower.tri(R)]
  y <- S[lower.tri(S)]
  x_mean <- mean(x)
  y_mean <- mean(y)
  x_scale <- max(abs(x - x_mean))
  if (!is.finite(x_scale))
    stop("MLPE resistance values are outside floating-point range.", call. = FALSE)
  if (x_scale == 0)
    x_scale <- 1
  list(n = nrow(E), x = (x - x_mean) / x_scale, y = y - y_mean,
       x_mean = x_mean, y_mean = y_mean, x_scale = x_scale)
}

.mlpe_initial_values <- function(E, S, nonnegative)
{
  data <- .mlpe_pair_data(E, S)
  denominator <- sum(data$x^2)
  beta <- if (denominator > 0) sum(data$x * data$y) / denominator else 0
  if (nonnegative)
    beta <- max(0, beta)
  variance <- mean((data$y - beta * data$x)^2)
  if (!is.finite(variance) || variance <= 0)
    stop("MLPE requires positive residual variance.", call. = FALSE)
  c(alpha = data$y_mean - beta * data$x_mean / data$x_scale,
    beta = beta / data$x_scale, tau = -log(variance),
    rho = stats::qlogis(0.2))
}

# Conditional GLS coefficients and ML residual precision are analytic at a
# fixed shared-site correlation. Search only the correlation, using the same
# complete-pair covariance as mlpe(). No REML term is introduced.
.mlpe_profile_start <- function(E, S, nonnegative)
{
  data <- .mlpe_pair_data(E, S)
  m <- length(data$x)
  ind <- which(lower.tri(E), arr.ind = TRUE)
  U <- Matrix::sparseMatrix(i = rep(seq_len(m), 2), j = c(ind), x = 1)
  eig <- .get_mlpe_eigen(data$n)
  keep <- eig$values > max(eig$values) * .Machine$double.eps * data$n
  d <- eig$values[keep]
  rotation <- eig$vectors[, keep, drop = FALSE]
  # Apply the orthonormal projection through sparse U rather than creating
  # a dense number-of-pairs by number-of-sites basis.
  xp <- as.numeric(crossprod(rotation, as.numeric(crossprod(U, data$x)))) / sqrt(d)
  yp <- as.numeric(crossprod(rotation, as.numeric(crossprod(U, data$y)))) / sqrt(d)
  xr <- as.numeric(data$x - U %*% (rotation %*% (xp / sqrt(d))))
  yr <- as.numeric(data$y - U %*% (rotation %*% (yp / sqrt(d))))
  evaluate <- function(logit, detail = FALSE)
  {
    lambda <- 0.5 * exp(logit)
    v <- 1 + lambda * d
    denominator <- sum(xr^2) + sum(xp^2 / v)
    beta <- if (denominator > 0)
      (sum(xr * yr) + sum(xp * yp / v)) / denominator else 0
    if (nonnegative)
      beta <- max(0, beta)
    er <- yr - beta * xr
    ep <- yp - beta * xp
    quadratic <- sum(er^2) + sum(ep^2 / v)
    if (!is.finite(quadratic) || quadratic <= 0)
      stop("MLPE requires positive residual variance.", call. = FALSE)
    objective <- 0.5 * m * (1 + log(quadratic / m)) +
      0.5 * sum(log1p(lambda * d))
    score <- -0.5 * m * sum(lambda * d * ep^2 / v^2) / quadratic +
      0.5 * sum(lambda * d / v)
    if (!detail)
      return(c(objective = objective, score = score))
    list(phi = c(alpha = data$y_mean - beta * data$x_mean / data$x_scale,
                 beta = beta / data$x_scale,
                 tau = log(m / quadratic) + stats::plogis(-logit, log.p = TRUE),
                 rho = logit),
         parameter_scale = c(1, 1 / data$x_scale, 1, 1),
         objective = objective, upper_endpoint = logit >= 36)
  }
  grid <- seq(-36, 36, by = 1)
  values <- vapply(grid, evaluate, numeric(2))
  candidates <- grid[c(1L, length(grid), which.min(values[1, ]))]
  for (i in seq_len(length(grid) - 1L))
    if (values[2, i] < 0 && values[2, i + 1L] > 0)
      candidates <- c(candidates, stats::uniroot(
        function(x) evaluate(x)[2], grid[c(i, i + 1L)], tol = 1e-11)$root)
  best <- candidates[which.min(vapply(candidates,
                                      function(x) evaluate(x)[1], numeric(1)))]
  evaluate(best, TRUE)
}
