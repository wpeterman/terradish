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

test_that("BoxConstrainedBFGS supports objective-only Armijo trial steps", {
  target <- 1
  calls <- data.frame(par = numeric(), gradient = logical())

  fn <- function(par, gradient, hessian) {
    par <- as.matrix(par)
    calls <<- rbind(
      calls,
      data.frame(par = par[1, 1], gradient = isTRUE(gradient))
    )

    diff <- par - target
    out <- list(objective = c(0.5 * diff^2))
    if (gradient)
      out$gradient <- diff
    if (hessian)
      out$hessian <- matrix(1, 1, 1)
    out
  }

  diagnostics <- terradish:::.terradish_new_diagnostics()
  control <- NewtonRaphsonControl(
    maxit = 3,
    verbose = FALSE,
    ls.control = ArmijoControl(initial = 4, contraction = 0.5)
  )
  control$diagnostics <- diagnostics

  fit <- terradish:::BoxConstrainedBFGS(
    0,
    fn,
    control = control
  )

  expect_equal(c(fit$par), target, tolerance = 1e-8)
  expect_true(any(!calls$gradient))
  expect_true(any(calls$par != 0 & !calls$gradient))
  expect_true(any(calls$par == target & calls$gradient))
  expect_gt(diagnostics$line_search_trials, 0L)
  expect_equal(diagnostics$line_search_trials,
               diagnostics$line_search_objective_only_trials)
  expect_equal(diagnostics$line_search_gradient_trials, 0L)
})

test_that("BoxConstrainedNewton rejects ArmijoControl explicitly", {
  fn <- function(par, gradient, hessian) {
    out <- list(objective = sum(par^2))
    if (gradient)
      out$gradient <- matrix(2 * par, ncol = 1)
    if (hessian)
      out$hessian <- diag(2, length(par))
    out
  }

  expect_error(
    terradish:::BoxConstrainedNewton(
      c(1, 1),
      fn,
      control = NewtonRaphsonControl(
        maxit = 1,
        verbose = FALSE,
        ls.control = ArmijoControl()
      )
    ),
    "supported for BFGS"
  )
})
