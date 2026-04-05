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
