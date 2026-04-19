test_that("BoxConstrainedNewton keeps real-valued steps for slightly asymmetric Hessians", {
  target <- matrix(c(1, -2), ncol = 1)
  A <- diag(c(2, 1))
  skew <- matrix(c(0, 5, -5, 0), nrow = 2)

  fn <- function(par, gradient, hessian) {
    par <- as.matrix(par)
    diff <- par - target
    out <- list(objective = c(0.5 * t(diff) %*% A %*% diff))
    if (gradient)
      out$gradient <- A %*% diff
    if (hessian)
      out$hessian <- A + skew
    out
  }

  fit <- terradish:::BoxConstrainedNewton(
    c(0, 0),
    fn,
    control = NewtonRaphsonControl(maxit = 5, verbose = FALSE)
  )

  expect_true(all(is.finite(fit$par)))
  expect_false(is.complex(fit$par))
  expect_equal(c(fit$par), c(target), tolerance = 1e-8)
})

test_that("BoxConstrainedBFGS reuses the accepted line-search fit", {
  target <- 1
  call_counts <- new.env(parent = emptyenv())

  fn <- function(par, gradient, hessian) {
    par <- as.matrix(par)
    key <- formatC(par[1, 1], digits = 17, format = "fg")
    call_counts[[key]] <- if (is.null(call_counts[[key]])) 1L else call_counts[[key]] + 1L

    diff <- par - target
    out <- list(objective = c(0.5 * diff^2))
    if (gradient)
      out$gradient <- diff
    if (hessian)
      out$hessian <- matrix(1, 1, 1)
    out
  }

  fit <- terradish:::BoxConstrainedBFGS(
    0,
    fn,
    control = NewtonRaphsonControl(maxit = 3, verbose = FALSE)
  )

  expect_equal(c(fit$par), target, tolerance = 1e-8)
  final_key <- formatC(fit$par[1, 1], digits = 17, format = "fg")
  expect_equal(call_counts[[final_key]], 1L)
})
