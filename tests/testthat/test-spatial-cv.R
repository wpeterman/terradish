test_that("terradish_folds is reproducible and preserves the random seed", {
  pts <- matrix(seq_len(24), ncol = 2)
  set.seed(42)
  before <- .Random.seed
  spatial_1 <- terradish_folds(pts, k = 3, seed = 8, projected = TRUE)
  spatial_2 <- terradish_folds(pts, k = 3, seed = 8, projected = TRUE)
  expect_equal(spatial_1, spatial_2)
  expect_equal(.Random.seed, before)
  expect_equal(length(unique(spatial_1)), 3L)

  random <- terradish_folds(pts, k = 3, method = "random", seed = 9)
  expect_lte(max(table(random)) - min(table(random)), 1)
  expect_error(terradish_folds(pts, k = 3), "projected")
})

test_that("fixed-domain cross-validation scores trained and baseline models", {
  dat <- melip_fixture(1:9)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords)
  folds <- terradish_folds(dat$coords, k = 3, method = "random", seed = 1)
  checkpoint <- tempfile(fileext = ".rds")

  cv <- suppressWarnings(terradish_cv_folds(
    surface,
    melip.Fst ~ altitude,
    folds = folds,
    model = "ls",
    keep_fits = "slim",
    checkpoint = checkpoint,
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
  ))

  expect_s3_class(cv, "terradish_cv_folds")
  expect_equal(nrow(cv$results), 3L)
  expect_true(all(cv$results$error == ""))
  expect_true(all(is.finite(cv$results$loglik)))
  expect_equal(cv$results$loglik_gain,
               cv$results$loglik - cv$results$baseline_loglik)
  expect_equal(length(cv$fits), 3L)
  expect_true(all(vapply(cv$fits, function(x) is.null(x$submodels), logical(1))))
  expect_gt(length(capture.output(print(cv))), 0L)

  resumed <- terradish_cv_folds(
    surface, melip.Fst ~ altitude, folds = folds, model = "ls",
    keep_fits = "slim", checkpoint = checkpoint,
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
  )
  expect_equal(resumed$results, cv$results)
  expect_equal(names(resumed$fits), names(cv$fits))
})

test_that("cv_model_selection rejects different held-out sites", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  item_1 <- list(train_mod = fx$fit, cv_loglik = fx$fit$loglik,
                 train_index = 1:5, test_index = 6:8)
  item_2 <- list(train_mod = fx$fit, cv_loglik = fx$fit$loglik,
                 train_index = 2:6, test_index = c(1, 7, 8))
  expect_error(cv_model_selection(list(item_1, item_2)), "identical")
})

test_that("measurement profiling reports convergence metadata", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  alg <- suppressWarnings(terradish_algorithm(
    f = fx$fit$submodels$f_internal,
    g = leastsquares,
    s = fx$surface,
    S = fx$data$melip.Fst,
    theta = fx$fit$mle$theta_internal,
    gradient = FALSE,
    hessian = FALSE,
    partial = FALSE,
    measurement_control = NewtonRaphsonControl(maxit = 2, verbose = FALSE)
  ))
  expect_named(alg$subproblem, c("convergence", "iters"))
  expect_true(alg$subproblem$iters >= 1L)
})
