test_that("small objective changes cannot conceal a large projected gradient", {
  objective <- function(x, gradient, hessian) list(
    objective = sum(exp(x)), gradient = matrix(exp(x)), hessian = diag(c(exp(x)), nrow = length(x)))
  for (method in c("newton", "bfgs")) {
    optimizer <- if (method == "newton") BoxConstrainedNewton else BoxConstrainedBFGS
    control <- NewtonRaphsonControl(ftol = 10, del = 0.1,
      ls.control = if (method == "newton") HagerZhangControl(alphamax = 1) else ArmijoControl())
    expect_warning(fit <- optimizer(1, objective,
      control = control), "stalled")
    expect_equal(fit$convergence, 2)
    expect_equal(fit$criterion, "stalled")
  }
})

test_that("Newton steps use the free Hessian block at active bounds", {
  H <- matrix(c(100, 99, 99, 100), 2)
  objective <- function(x, gradient, hessian) list(
    objective = as.numeric(crossprod(x, H %*% x) / 2 + crossprod(c(100, -100), x)),
    gradient = H %*% x + c(100, -100), hessian = H)
  fit <- BoxConstrainedNewton(c(0, 0), objective, lower = c(0, -Inf))
  expect_equal(c(fit$par), c(0, 1), tolerance = 1e-10)
  expect_equal(fit$convergence, 0)
  expect_equal(fit$max_abs_projected_gradient, 0, tolerance = 1e-10)
  expect_lte(fit$iters, 3)
})
