# The solver warning guard turns CHOLMOD "not positive definite" warnings into
# errors. Matrix raises those warnings from inside the compiled factorization,
# so the guard must let the solver call return before it signals the error.
# Raising the error from the handler jumped out of CHOLMOD before it reset its
# shared workspace, and a later supernodal factorization could then loop
# forever or crash.

test_that("solver warning guard raises only after the solver call returns", {
  finished <- FALSE
  expect_error(
    terradish:::.terradish_with_solver_warning_guard(
      {
        warning("CHOLMOD warning 'matrix not positive definite'")
        finished <- TRUE
      },
      context = "Test setup."
    ),
    "non-positive-definite reduced Laplacian"
  )
  # The expression ran to completion, so the error came after it returned.
  expect_true(finished)
})

test_that("solver warning guard keeps its message when the solver also errors", {
  # Matrix's Cholesky() signals its own error after the CHOLMOD warning once
  # the warning is muffled; the guard reports its own message instead.
  expect_error(
    terradish:::.terradish_with_solver_warning_guard(
      {
        warning("CHOLMOD warning 'matrix not positive definite'")
        stop("leading principal minor is not positive")
      },
      context = "Test setup."
    ),
    "^Test setup\\. CHOLMOD warning .*non-positive-definite reduced Laplacian\\.$"
  )
})

test_that("solver warning guard passes other errors through unchanged", {
  expect_error(
    terradish:::.terradish_with_solver_warning_guard(
      stop("unrelated failure"),
      context = "Test setup."
    ),
    "^unrelated failure$"
  )
})

test_that("solver warning guard passes non-fatal warnings and values through", {
  expect_warning(
    value <- terradish:::.terradish_with_solver_warning_guard(
      {
        warning("an unrelated warning")
        42
      },
      context = "Test setup."
    ),
    "an unrelated warning"
  )
  expect_identical(value, 42)
})

test_that("supernodal factorization still works after a guarded failure", {
  skip_on_cran()

  # A shifted 2-D grid Laplacian is sparse and positive definite and has many
  # supernodes. Making a late diagonal entry negative makes a supernode past
  # the first one fail, which is the case that left CHOLMOD's workspace dirty.
  m <- 40L
  tri <- Matrix::bandSparse(
    m,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, m - 1), rep(2, m), rep(-1, m - 1))
  )
  identity <- Matrix::Diagonal(m)
  spd <- Matrix::kronecker(identity, tri) + Matrix::kronecker(tri, identity) +
    Matrix::Diagonal(m * m, 1e-3)
  spd <- methods::as(Matrix::forceSymmetric(spd), "dsCMatrix")
  bad <- spd
  k <- floor(0.7 * nrow(bad))
  bad[k, k] <- -50
  bad <- methods::as(Matrix::forceSymmetric(bad), "dsCMatrix")

  expect_error(
    terradish:::.terradish_with_solver_warning_guard(
      Matrix::Cholesky(bad, LDL = FALSE, super = TRUE, perm = FALSE),
      context = "Test setup."
    )
  )

  # Before the fix this second factorization could hang or crash R.
  factor <- Matrix::Cholesky(spd, LDL = FALSE, super = TRUE, perm = FALSE)
  rhs <- seq_len(nrow(spd)) / nrow(spd)
  solution <- as.numeric(Matrix::solve(factor, rhs))
  expect_lt(max(abs(as.numeric(spd %*% solution) - rhs)), 1e-6)
})
