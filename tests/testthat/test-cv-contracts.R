cv_contract_fixture <- function() {
  set.seed(19)
  r <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 9, ymin = 0, ymax = 9,
                    crs = "EPSG:3857")
  terra::values(r) <- rnorm(81); names(r) <- "x"
  z <- r; terra::values(z) <- rnorm(81); names(z) <- "z"
  r <- c(r, z)
  graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 80, 7)),
                                directions = 4, saveStack = TRUE)
  S <- simulate_covariance_response(.5, ~x, graph, tau = 1, sigma = .1,
                                     nu = 200, seed = 21)$covariance
  list(graph = graph, S = S, folds = rep(1:3, 4))
}

test_that("mixed conductance factories and Gaussian-only baselines work", {
  fx <- cv_contract_fixture(); S <- fx$S
  gaussian <- gaussian_smoothed_loglinear_conductance(fx$graph, sigma_lower = .3,
    sigma_upper = 1.3)
  cv <- terradish_cv_folds(fx$graph, list(linear = S ~ x, gaussian = S ~ x,
    spline = S ~ s(x, df = 3)), fx$folds, model = wishart_covariance, nu = 25,
    conductance_model = list(linear = loglinear_conductance, gaussian = gaussian,
                            spline = smooth_loglinear_conductance))
  expect_true(all(is.finite(cv$results$baseline_loglik)))
  expect_true(all(c("n_folds_ok", "common_folds_total", "mean_difference",
                    "se_difference", "within_one_se") %in% names(cv$summary)))
  single <- terradish_cv_folds(fx$graph, S ~ x, fx$folds, model = wishart_covariance,
    nu = 25, conductance_model = gaussian)
  expect_true(all(is.finite(single$results$baseline_loglik)))
})

test_that("checkpoint signatures reject changed responses and repeated folds pool", {
  fx <- cv_contract_fixture(); S <- fx$S
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  folds <- list(fx$folds, rev(fx$folds))
  cv <- terradish_cv_folds(fx$graph, list(one = S ~ x, two = S ~ x), folds,
    model = wishart_covariance, nu = 25, checkpoint = path)
  resumed <- terradish_cv_folds(fx$graph, list(one = S ~ x, two = S ~ x), folds,
    model = wishart_covariance, nu = 25, checkpoint = path)
  expect_equal(cv$results, resumed$results)
  expect_equal(cv$repeats, 2L)
  expect_equal(cv$summary$n_folds_ok, c(6L, 6L))
  expect_true(all(is.finite(cv$summary$between_repeat_sd)))
  expect_equal(cv$summary$mean_difference, c(0, 0))
  S <- S * 1.01
  expect_error(terradish_cv_folds(fx$graph, list(one = S ~ x, two = S ~ x), folds,
    model = wishart_covariance, nu = 25, checkpoint = path), "not compatible")
})

test_that("fixed predictive scores retain training nuisance parameters", {
  fx <- cv_contract_fixture(); S <- fx$S
  model <- wishart_covariates(data.frame(z = seq_len(12)))
  expect_error(terradish_cv_folds(fx$graph, list(base = S ~ x, env = S ~ x), fx$folds,
    model = list(base = wishart_covariance, env = model), nu = 25), "one measurement")
  cv <- terradish_cv_folds(fx$graph, list(base = S ~ x, env = S ~ x), fx$folds,
    model = list(base = wishart_covariance, env = model), nu = 25,
    nuisance = "fixed", keep_fits = "full")
  expect_true(all(cv$results$nuisance_iterations == 0L))
  fit <- cv$fits[["env::1::1"]]
  test <- which(fx$folds == 1)
  graph <- .terradish_subset_graph_focal(fx$graph, test)
  E <- terradish_distance(matrix(coef(fit), 1), ~x, graph,
                           covariance = TRUE)$covariance[, , 1]
  g <- .terradish_measurement_model(model, subset = test)
  expected <- -g(E, S[test, test], phi = c(fit$fit$phi), nu = 25)$objective
  expect_equal(cv$results$loglik[cv$results$model == "env" & cv$results$fold == "1"], expected)
  common <- cv$results[cv$results$common_fold, ]
  best <- cv$summary$model[1]
  best_score <- setNames(common$loglik[common$model == best],
                         common$fold[common$model == best])
  for (name in cv$summary$model) {
    rows <- common[common$model == name, ]
    delta <- rows$loglik - best_score[rows$fold]
    summary_row <- cv$summary[cv$summary$model == name, ]
    expect_equal(summary_row$mean_difference, mean(delta))
    expect_equal(summary_row$se_difference, sd(delta) / sqrt(length(delta)))
  }
})

test_that("fixed no-structure predictions do not require identified conductance", {
  fx <- cv_contract_fixture(); S <- fx$S
  f <- loglinear_conductance(~x, fx$graph$x)
  fit <- list(mle = list(theta_internal = NULL),
              fit = list(boundary = TRUE, no_structure_boundary = TRUE),
              submodels = list(f_internal = f))
  theta <- .terradish_cv_prediction_theta(fit, "fixed")
  expect_equal(theta, attr(f, "default"))
  expect_error(.terradish_cv_prediction_theta(fit, "reprofile"), "no conductance")
  a <- .terradish_cv_score(theta, S ~ x, fx$graph, S, loglinear_conductance,
    wishart_covariance, 25, FALSE, NULL, phi = c(0, log(.2)))
  b <- .terradish_cv_score(theta + 1, S ~ x, fx$graph, S, loglinear_conductance,
    wishart_covariance, 25, FALSE, NULL, phi = c(0, log(.2)))
  expect_true(is.finite(a$loglik))
  expect_equal(a$loglik, b$loglik)
})

test_that("Wishart CV totals scale with nu without changing rank", {
  fx <- cv_contract_fixture(); S <- fx$S
  a <- terradish_cv_folds(fx$graph, list(linear = S ~ x, both = S ~ x + z),
    fx$folds, model = wishart_covariance, nu = 25)
  b <- terradish_cv_folds(fx$graph, list(linear = S ~ x, both = S ~ x + z),
    fx$folds, model = wishart_covariance, nu = 500)
  expect_equal(b$summary$model, a$summary$model)
  expect_true(all(is.finite(a$summary$common_folds_total)))
  expect_gt(min(a$summary$n_common_folds), 0)
  expect_equal(b$summary$n_common_folds, a$summary$n_common_folds)
  expect_equal(b$summary$common_folds_total, 20 * a$summary$common_folds_total,
               tolerance = 1e-6)
  expect_equal(b$summary$total_loglik, 20 * a$summary$total_loglik, tolerance = 1e-6)
})

test_that("one failed fold makes the full total NA and preserves common-fold rankings", {
  fx <- cv_contract_fixture(); S <- fx$S
  selected <- which(fx$folds == 1)
  failure_value <- S[selected[1], selected[2]]
  failing <- function(E, S, phi, ...) {
    if (nrow(S) == 4 && identical(unname(S[1, 2]), unname(failure_value)))
      stop("deliberate held-out failure")
    if (missing(phi)) wishart_covariance(E, S, ...)
    else wishart_covariance(E, S, phi = phi, ...)
  }
  class(failing) <- class(wishart_covariance)
  attr(failing, "base_model") <- "wishart_covariance"
  expect_warning(cv <- terradish_cv_folds(fx$graph,
    list(good = S ~ x, bad = S ~ x), fx$folds,
    model = list(good = wishart_covariance, bad = failing), nu = 25,
    nuisance = "fixed", baseline = FALSE), "Failed cross-validation folds")
  bad <- cv$summary[cv$summary$model == "bad", ]
  expect_true(is.na(bad$total_loglik))
  expect_equal(bad$n_folds_ok, 2L)
  expect_equal(cv$summary$n_common_folds, c(2L, 2L))
  expect_true(all(is.finite(cv$summary$common_folds_total)))
})
