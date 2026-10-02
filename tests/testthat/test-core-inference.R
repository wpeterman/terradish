inference_fixture <- function() {
  set.seed(251)
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
  terra::values(r) <- rnorm(64)
  names(r) <- "x"
  coords <- terra::xyFromCell(r, c(1, 8, 13, 20, 27, 35, 42, 49, 57, 64))
  graph <- conductance_surface(r, coords, directions = 4, saveStack = TRUE)
  E <- terradish_distance(matrix(.4, 1), ~x, graph, covariance = TRUE)$covariance[, , 1]
  z <- data.frame(a = rnorm(10), b = rnorm(10))
  g <- wishart_covariates(z["a"])
  K <- attr(g, "kernel_covariates")[, , 1]
  set.seed(51)
  S <- rWishart(1, 1000, .8 * E + .6 * K + .2 * diag(10))[, , 1] / 1000
  list(graph = graph, S = S, g = g, z = z)
}

test_that("inference methods agree across full and slim fits", {
  fx <- inference_fixture(); S <- fx$S
  fit <- terradish(S ~ x, fx$graph, measurement_model = fx$g, nu = 1000)
  expect_equal(unname(diag(vcov(fit))), unname(summary(fit)$ztable[, "Std. Error"]^2))
  slim <- slim_terradish(fit)
  expect_equal(vcov(slim), vcov(fit))
  expect_equal(confint(slim), confint(fit))
  # sample size is the 10 focal sites, not the 45 site pairs
  expect_equal(nobs(fit), 10)
  expect_equal(as.numeric(BIC(fit)), -2 * fit$loglik + fit$df * log(10))
  expect_equal(rep(as.numeric(BIC(fit)), 2), suppressWarnings(aic_table(list(fit, fit), BIC = TRUE))$BIC)
  expect_error(confint(fit, parm = "missing"), "Unknown")
  expect_error(confint(fit, level = 1), "between")
  expect_warning(short <- terradish(S ~ x, fx$graph, measurement_model = fx$g,
    nu = 1000, control = NewtonRaphsonControl(maxit = 2)), "maxit")
  expect_false(short$convergence$code == 0)
  expect_warning(table <- aic_table(list(short, short)), "converg")
  expect_false(any(table$converged))
})

test_that("rescaling nu agrees with refitting and scales all covariance", {
  fx <- inference_fixture(); S <- fx$S
  fit <- terradish(S ~ x, fx$graph, measurement_model = fx$g, nu = 1000)
  rescaled <- terradish_rescale_nu(fit, 500)
  refit <- terradish(S ~ x, fx$graph, measurement_model = fx$g, nu = 500)
  expect_equal(coef(rescaled), coef(refit), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(rescaled)), as.numeric(logLik(refit)), tolerance = 1e-7)
  expect_equal(vcov(rescaled), vcov(refit), tolerance = 1e-5)
  expect_equal(rescaled$fit$phi_vcov_joint, refit$fit$phi_vcov_joint, tolerance = 1e-5)
  expect_equal(rescaled$nu_original, 1000)
  expect_equal(terradish_rescale_nu(rescaled, 200)$nu_original, 1000)
  expect_error(aic_table(list(fit, rescaled)), "same effective")
})

test_that("joint nuisance covariance and ratio uncertainty match numerical curvature", {
  skip_if_not_installed("numDeriv")
  fx <- inference_fixture(); S <- fx$S
  fit <- terradish(S ~ x, fx$graph, measurement_model = fx$g, nu = 1000)
  parameters <- c(coef(fit), fit$fit$phi[, 1])
  objective <- function(p) {
    E <- terradish_distance(matrix(p[1], 1), ~x, fx$graph,
                           covariance = TRUE)$covariance[, , 1]
    fx$g(E, S, phi = p[-1], nu = 1000, gradient = FALSE,
         hessian = FALSE, partial = FALSE)$objective
  }
  V <- solve(numDeriv::hessian(objective, parameters))
  expect_equal(unname(fit$fit$phi_vcov_joint), unname(V[-1, -1]), tolerance = 1e-4)
  ratio <- terradish_ibe_ratio(fit)
  derivative <- numDeriv::grad(function(p) p[3] / p[2], parameters)
  expected <- sqrt(as.numeric(crossprod(derivative, V %*% derivative)))
  expect_equal(ratio$SE, expected, tolerance = 1e-4)
  D <- dist_from_cov(S)
  mlpe_fit <- terradish(D ~ x, fx$graph,
    measurement_model = mlpe_covariates(pairwise_endpoint_covariates(fx$z["a"])))
  expect_gt(terradish_ibe_ratio(mlpe_fit)$ratio, 0)
  expect_gt(ratio$ratio, 0)
  expect_true(is.finite(terradish_ibe_ratio(mlpe_fit)$SE))
})

test_that("kernel boundary flags preserve conductance derivatives and use chi-bar tests", {
  expect_equal(.chibar_probability(3, 1), .5 * pchisq(3, 1, lower.tail = FALSE))
  expect_equal(.chibar_probability(3, 2), .5 * pchisq(3, 1, lower.tail = FALSE) +
                 .25 * pchisq(3, 2, lower.tail = FALSE))
  expect_equal(.chibar_probability(0, 2), 1)
  fx <- inference_fixture(); S <- fx$S
  E <- diag(10)
  boundary <- fx$g(E, S, phi = c(1, 0, log(.2)), nu = 1000)
  expect_true(boundary$boundary)
  expect_false(boundary$no_structure_boundary)
  full <- terradish(S ~ x, fx$graph, measurement_model = fx$g, nu = 1000)
  null <- terradish(S ~ x, fx$graph, measurement_model = wishart_covariance, nu = 1000)
  expect_warning(table <- anova(null, full), "zero bounds")
  expect_match(attr(table, "heading")[1], "chi-bar")
  expect_equal(table[2, "Pr(>Chi)"], .5 * pchisq(table[2, "ChiSq"], 1, lower.tail = FALSE))
  different <- full; different$comparison$measurement_columns <- matrix(1, 45, 1,
    dimnames = list(NULL, "other")); different$df <- different$df + 1
  expect_error(anova(full, different), "not nested")
  different <- full; different$comparison$graph <- "different"
  expect_error(anova(null, different), "same graph")
  missing_factory <- full
  missing_factory$submodels$f_factory <- NULL
  expect_error(anova(null, missing_factory), "retained model factories")
  null$submodels$f_factory <- NULL
  expect_error(anova(null, missing_factory), "retained model factories")
  bounded <- full; bounded$fit$phi[2, 1] <- 0; bounded$fit$boundary <- TRUE
  s <- summary(bounded)
  expect_true(s$phi_at_bound[2])
  expect_true(is.na(s$phi_table[2, 3]))
  expect_true(is.finite(s$phi_table[2, 4]))
  expect_false(is.null(s$ztable))
})

test_that("MLPE near-zero correlation has no Wald output", {
  E <- diag(5); S <- as.matrix(dist(1:5))
  raw <- mlpe(E, S, phi = c(0, 1, 0, -9))
  expect_true(raw$rho_boundary)
  fx <- inference_fixture(); S <- dist_from_cov(fx$S)
  fit <- terradish(S ~ x, fx$graph, measurement_model = mlpe)
  fit$fit$rho_boundary <- TRUE
  original_df <- fit$df
  s <- summary(fit)
  expect_equal(unname(s$phi["rho"]), 0)
  expect_true(all(is.na(s$phi_table["rho", -1])))
  expect_equal(s$df, original_df)
  expect_match(s$rho_note, "zero")
})
