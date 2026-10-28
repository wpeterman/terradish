## PSOCK startup and cross-core fits are exercised in the extended test suite.
## Keep CRAN's routine check focused on single-core numerical contracts.
skip_on_cran()

test_that("source-loaded worker setup preserves each worker's directory", {
  if (!requireNamespace("pkgload", quietly = TRUE) ||
      !pkgload::is_dev_package("terradish"))
    skip("source-loaded workers require the development source tree")
  if (!"local_mocked_bindings" %in% getNamespaceExports("testthat"))
    skip("worker instrumentation requires testthat's local_mocked_bindings()")

  # The source fallback is supported when called from the package root.
  # Restore the test runner's directory on both success and error.
  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  setwd(getNamespaceInfo(asNamespace("terradish"), "path"))

  worker_checks <- 0L
  testthat::local_mocked_bindings(
    clusterEvalQ = function(cl, expr) {
      expression <- substitute(expr)
      result <- parallel::clusterCall(cl, function(code) {
        # Start somewhere other than the source root so an unintended setwd()
        # is observable, and restore the worker even if initialization fails.
        oldwd <- getwd()
        on.exit(setwd(oldwd), add = TRUE)
        scratch <- tempfile("worker-directory-")
        dir.create(scratch)
        on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
        setwd(scratch)
        before <- getwd()
        value <- eval(code, envir = globalenv())
        stopifnot(identical(getwd(), before))
        value
      }, expression)
      worker_checks <<- worker_checks + 1L
      result
    },
    .package = "terradish"
  )

  dat <- melip_fixture(1:6, aggregate_factor = 5L)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta <- rbind(c(altitude = -0.2, forestcover = 0.2),
                 c(altitude = 0.2, forestcover = -0.2))
  response <- dat$melip.Fst
  serial <- terradish_grid(theta, response ~ altitude + forestcover, surface,
                           loglinear_conductance, leastsquares, cores = 1)
  worker <- terradish_grid(theta, response ~ altitude + forestcover, surface,
                           loglinear_conductance, leastsquares, cores = 2)
  distance1 <- terradish_distance(theta, ~ altitude + forestcover, surface,
                                  loglinear_conductance, cores = 1)
  distance2 <- terradish_distance(theta, ~ altitude + forestcover, surface,
                                  loglinear_conductance, cores = 2)
  expect_equal(worker_checks, 2L)
  expect_true(all(is.finite(worker$loglik)))
  expect_true(all(is.finite(distance2$distance)))
  expect_equal(worker$loglik, serial$loglik, tolerance = 1e-8)
  expect_equal(distance2$distance, distance1$distance, tolerance = 1e-8)

  # The new nuisance-profile helpers must also be available to source-loaded
  # workers. Compare all supported likelihoods on the same small graph.
  S <- diag(nrow(response)) + matrix(0.1, nrow(response), nrow(response))
  for (model in list(mlpe, wishart_covariance, generalized_wishart)) {
    y <- if (identical(model, mlpe)) response else
      if (identical(model, wishart_covariance)) S else dist_from_cov(S)
    serial <- terradish_grid(theta, y ~ altitude + forestcover, surface,
      loglinear_conductance, model, nu = 30, cores = 1)
    worker <- terradish_grid(theta, y ~ altitude + forestcover, surface,
      loglinear_conductance, model, nu = 30, cores = 2)
    expect_true(all(is.finite(worker$loglik)))
    expect_equal(worker$loglik, serial$loglik, tolerance = 1e-8)
  }
  expect_equal(worker_checks, 5L)
})

test_that("terradish respects cores argument", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit1 <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
              cores = 1)
  )
  fit2 <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
              cores = 2)
  )

  expect_equal(unname(coef(fit1)), unname(coef(fit2)), tolerance = 1e-6)
  expect_equal(unname(fit1$loglik), unname(fit2$loglik), tolerance = 1e-6)
})

test_that("terradish_grid and terradish_distance match across cores", {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("terradish"))
  {
    skip("parallel worker consistency is validated against the installed package")
  }

  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta <- as.matrix(expand.grid(forestcover = c(-0.3, 0.3),
                                 altitude = c(-0.3, 0.3)))

  grid1 <- terradish_grid(theta, melip.Fst ~ altitude + forestcover,
                          data = surface,
                          conductance_model = loglinear_conductance,
                          measurement_model = leastsquares,
                          cores = 1)
  grid2 <- terradish_grid(theta, melip.Fst ~ altitude + forestcover,
                          data = surface,
                          conductance_model = loglinear_conductance,
                          measurement_model = leastsquares,
                          cores = 2)

  dist1 <- terradish_distance(theta, ~ altitude + forestcover,
                              data = surface,
                              conductance_model = loglinear_conductance,
                              cores = 1)
  dist2 <- terradish_distance(theta, ~ altitude + forestcover,
                              data = surface,
                              conductance_model = loglinear_conductance,
                              cores = 2)

  expect_equal(grid1$loglik, grid2$loglik, tolerance = 1e-8)
  expect_equal(dist1$distance, dist2$distance, tolerance = 1e-8)
})
