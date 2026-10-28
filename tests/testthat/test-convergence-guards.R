test_that("AIC requires both convergence records, including for slim fits", {
  make_fit <- function(ll) list(formula = ~x, dim = c(focal = 8),
    loglik = ll, aic = -2 * ll + 4, df = 2,
    convergence = list(code = 0L),
    fit = list(subproblem = list(convergence = 0L)),
    submodels = list(g = mlpe))
  good <- make_fit(-10)
  for (part in c("outer", "nuisance", "missing_outer", "missing_nuisance")) {
    bad <- make_fit(-9)
    if (part == "outer") bad$convergence$code <- 1L
    if (part == "nuisance") bad$fit$subproblem$convergence <- 2L
    if (part == "missing_outer") bad$convergence <- NULL
    if (part == "missing_nuisance") bad$fit$subproblem <- NULL
    expect_false(.terradish_fit_converged(bad))
    for (criterion in c("AIC", "AICc", "BIC")) {
      expect_warning(tab <- aic_table(list(good, bad), mod_names = c("good", "bad"),
        AICc = criterion == "AICc", BIC = criterion == "BIC"), "Some models have not converged")
      expect_true(tab$converged[tab$model == "good"])
      expect_false(tab$converged[tab$model == "bad"])
    }
  }
  expect_true(.terradish_fit_converged(good))
  expect_no_warning(tab <- aic_table(list(good, make_fit(-9))))
  expect_true(all(tab$converged))
})

test_that("CV rejects failed outer and nuisance training optimization", {
  set.seed(19)
  r <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 9, ymin = 0, ymax = 9,
                   crs = "EPSG:3857")
  terra::values(r) <- rnorm(81); names(r) <- "x"
  graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 80, 7)),
                                directions = 4, saveStack = TRUE)
  S <- simulate_covariance_response(.5, ~x, graph, tau = 1, sigma = .1,
                                    nu = 200, seed = 21)$covariance
  folds <- rep(1:3, 4)
  # The environmental covariance model has an extra nuisance component and
  # does not use the ordinary two-component profile start. One Newton step
  # genuinely leaves this solve incomplete instead of starting at its optimum.
  nuisance_model <- wishart_covariates(data.frame(environment = seq_len(nrow(S))))
  # Uniform fits have no outer parameters: code 0 there must not mask a failed
  # nuisance solve. Fixed prediction does not optimize those parameters again.
  for (part in c("outer", "nuisance")) {
    args <- list(data = graph, formulas = if (part == "outer") S ~ x else S ~ 1,
      folds = folds, model = if (part == "outer") wishart_covariance else nuisance_model, nu = 25,
      nuisance = "fixed", baseline = FALSE)
    args[[if (part == "outer") "control" else "measurement_control"]] <-
      NewtonRaphsonControl(maxit = 1)
    cv <- expect_fixture_warnings(do.call(terradish_cv_folds, args),
      c("`maxit` reached for Newton steps", "Failed cross-validation folds"))
    expect_equal(cv$summary$n_folds_ok, 0L)
    expect_true(is.na(cv$summary$total_loglik))
    expect_true(all(grepl("training fit did not converge", cv$results$error)))
    expect_true(all(is.na(cv$results$loglik)))
  }
  good <- terradish_cv_folds(graph, S ~ 1, folds, model = wishart_covariance,
    nu = 25, nuisance = "fixed", baseline = FALSE, keep_fits = "full")
  expect_equal(good$summary$n_folds_ok, 3L)
  expect_true(is.finite(good$summary$total_loglik))
  expect_true(all(vapply(good$fits, .terradish_fit_converged, logical(1))))
  for (fit in good$fits)
    expect_true(.terradish_fit_converged(slim_terradish(fit)))
})
