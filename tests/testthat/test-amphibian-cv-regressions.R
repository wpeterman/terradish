amphibian_cv_fixture <- function() {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8,
                   ymin = 0, ymax = 8, crs = "EPSG:3857")
  terra::values(r) <- rep(0:1, each = 32)
  names(r) <- "cover"
  graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 58, 5)), directions = 4)
  graph$x$cover <- factor(graph$x$cover, levels = 0:1, labels = c("forest", "matrix"))
  contrasts(graph$x$cover) <- contr.treatment(levels(graph$x$cover), base = 1)
  response <- simulate_covariance_response(0.5, ~cover, graph,
    tau = 1, sigma = .2, nu = 200, seed = 71)$covariance
  list(graph = graph, S = response, folds = rep(1:3, 4))
}

test_that("factor checkpoint identity retains placement, reference and numeric compatibility", {
  x <- data.frame(cover = factor(c("forest", "matrix", "forest", "matrix")), z = 1:4)
  old <- .cv_signature(.cv_covariate_fingerprint(x))
  moved <- x; moved$cover <- moved$cover[c(2, 1, 3, 4)]
  expect_false(identical(old, .cv_signature(.cv_covariate_fingerprint(moved))))
  reference <- x; reference$cover <- relevel(reference$cover, "matrix")
  expect_false(identical(old, .cv_signature(.cv_covariate_fingerprint(reference))))
  encoded <- x; contrasts(encoded$cover) <- contr.sum(2)
  expect_false(identical(old, .cv_signature(.cv_covariate_fingerprint(encoded))))
  numeric <- data.frame(a = 1:4, b = 5:8)
  expect_identical(.cv_covariate_fingerprint(numeric), colSums(numeric))
})

test_that("categorical CV completes and refuses a changed spatial class layout", {
  fx <- amphibian_cv_fixture(); S <- fx$S
  checkpoint <- tempfile(fileext = ".rds"); on.exit(unlink(checkpoint))
  cv <- suppressWarnings(terradish_cv_folds(fx$graph, S ~ cover, fx$folds,
    model = wishart_covariance, nu = 200, baseline = FALSE, checkpoint = checkpoint,
    nuisance = "fixed", leverage = FALSE))
  expect_true(all(is.finite(cv$results$loglik)))
  expect_true(all(cv$results$nuisance_iterations == 0L))
  moved <- fx$graph
  moved$x$cover <- rev(moved$x$cover)
  expect_error(suppressWarnings(terradish_cv_folds(moved, S ~ cover, fx$folds,
    model = wishart_covariance, nu = 200, baseline = FALSE, checkpoint = checkpoint,
    nuisance = "fixed", leverage = FALSE)), "not compatible")
})

test_that("uniform predictive scores use unit conductance and fixed training nuisance", {
  fx <- amphibian_cv_fixture(); S <- fx$S
  train <- .terradish_subset_graph_focal(fx$graph, 1:8)
  fit <- terradish(S[1:8, 1:8] ~ 1, train, measurement_model = wishart_covariance,
    nu = 200, leverage = FALSE)
  expect_false(.no_structure_boundary(fit$fit))
  theta <- .terradish_cv_prediction_theta(fit, "fixed")
  expect_true(all(fit$submodels$f_internal(theta)$conductance == 1))
  scored <- .terradish_cv_score(theta, S ~ 1, train, S[1:8, 1:8],
    loglinear_conductance, wishart_covariance, 200, TRUE, NULL, phi = fit$fit$phi)
  expect_equal(scored$loglik, as.numeric(logLik(fit)), tolerance = 1e-8)
  expect_equal(scored$iterations, 0L)
  cv <- terradish_cv_folds(fx$graph, S ~ 1, fx$folds, model = wishart_covariance,
    nu = 200, baseline = FALSE, nuisance = "fixed", leverage = FALSE)
  expect_true(all(is.finite(cv$results$loglik)))
  expect_true(all(cv$results$nuisance_iterations == 0L))
})

test_that("public CV forwards an explicitly requested optimizer exactly once", {
  fx <- amphibian_cv_fixture(); S <- fx$S
  cv <- suppressWarnings(terradish_cv_folds(fx$graph, S ~ cover, fx$folds,
    model = wishart_covariance, nu = 200, baseline = FALSE, nuisance = "fixed",
    optimizer = "bfgs", keep_fits = "slim"))
  expect_true(all(cv$results$error == ""))
  expect_true(all(is.finite(cv$results$loglik)))
  expect_length(cv$fits, 3L)
  expect_true(all(vapply(cv$fits, function(fit)
    identical(fit$call$optimizer, "bfgs"), logical(1))))
})

test_that("explicit optimizer changes invalidate CV checkpoints", {
  fx <- amphibian_cv_fixture(); S <- fx$S
  path <- tempfile(fileext = '.rds'); on.exit(unlink(path))
  cv <- suppressWarnings(terradish_cv_folds(fx$graph, S ~ cover, fx$folds,
    model = wishart_covariance, nu = 200, baseline = FALSE, nuisance = 'fixed',
    optimizer = 'bfgs', checkpoint = path))
  expect_true(all(cv$results$error == ''))
  expect_error(suppressWarnings(terradish_cv_folds(fx$graph, S ~ cover, fx$folds,
    model = wishart_covariance, nu = 200, baseline = FALSE, nuisance = 'fixed',
    optimizer = 'newton', checkpoint = path)), 'not compatible')
})

test_that("Gaussian profiles refuse failed final measurement convergence", {
  r <- terra::rast(nrows = 7, ncols = 7, xmin = 0, xmax = 7,
                   ymin = 0, ymax = 7, crs = "EPSG:3857")
  set.seed(72); terra::values(r) <- rnorm(49); names(r) <- "x"
  graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 45, 6)), directions = 4)
  factory <- gaussian_smoothed_loglinear_conductance(graph, sigma_lower = .3, sigma_upper = 1.1)
  S <- simulate_covariance_response(c(x = .6, sigma.x = .8), ~x, graph,
    conductance_model = factory, tau = 1, sigma = .2, nu = 500, seed = 73)$covariance
  fit <- terradish(S ~ x, graph, conductance_model = factory,
    measurement_model = wishart_covariance, nu = 500, optimizer = "bfgs", leverage = FALSE)
  failed <- fit; failed$fit$subproblem$convergence <- 2L
  expect_error(gaussian_scale_profile(failed, "x", n = 5), "measurement parameters")
  local_mocked_bindings(terradish = function(...) {
    list(convergence = list(code = 0L), fit = list(subproblem = list(convergence = 2L)),
         loglik = fit$loglik)
  }, .package = "terradish")
  expect_error(gaussian_scale_profile(fit, "x", n = 5), "Profile optimization failed")
})
