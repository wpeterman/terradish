test_that("aic_table ranks fitted terradish models across supported criteria", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst

  fit1 <- suppressWarnings(
    terradish(melip.Fst ~ altitude,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )
  fit2 <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  aic_tab <- aic_table(list(fit1, fit2))
  expect_true(is.data.frame(aic_tab))
  expect_true(all(c("model", "K", "AIC", "Delta_AIC", "AIC_wt", "Cum.wt", "loglik") %in% names(aic_tab)))

  aicc_tab <- aic_table(list(fit1, fit2), AICc = TRUE)
  expect_true(is.data.frame(aicc_tab))
  expect_true(all(c("AICc", "Delta_AICc", "AICc_wt") %in% names(aicc_tab)))

  bic_tab <- aic_table(list(fit1, fit2), BIC = TRUE)
  expect_true(is.data.frame(bic_tab))
  expect_true(all(c("BIC", "Delta_BIC", "BIC_wt") %in% names(bic_tab)))
  # BIC penalizes with the number of focal sites, as AICc does, not site pairs
  n_sites <- fit1$dim[["focal"]]
  expect_equal(sort(bic_tab$BIC),
               sort(round(-2 * c(fit1$loglik, fit2$loglik) + c(fit1$df, fit2$df) * log(n_sites), 4)))

  expect_error(aic_table(list(fit1, fit2), AICc = TRUE, BIC = TRUE),
               "Set only one")
})


test_that("default model labels append [mlpe:n] when pairwise covariates are present", {
  set.seed(1)
  z <- pairwise_endpoint_covariates(matrix(rnorm(12), nrow = 4, ncol = 3))
  g_joint <- mlpe_covariates(z)

  fit_base <- list(
    formula = stats::as.formula(~ altitude),
    dim = c(vertices = 20, focal = 6, edge = 30),
    loglik = -10,
    aic = 30,
    df = 2,
    submodels = list(g = mlpe)
  )
  fit_joint <- list(
    formula = stats::as.formula(~ altitude + forestcover),
    dim = c(vertices = 20, focal = 6, edge = 30),
    loglik = -9,
    aic = 28,
    df = 3,
    submodels = list(g = g_joint)
  )

  tab <- aic_table(list(fit_base, fit_joint))
  expect_true(any(grepl("\\[mlpe:3\\]$", tab$model)))
  expect_true(any(tab$model == "altitude"))
  expect_true(any(tab$model == "altitude + forestcover [mlpe:3]"))

  custom <- aic_table(list(fit_base, fit_joint), mod_names = c("base", "joint"))
  expect_setequal(custom$model, c("base", "joint"))

})

test_that("aic_table rejects incomparable likelihoods, responses, and nu", {
  make_fit <- function(model, response, nu = NULL) {
    list(
      formula = stats::as.formula("response ~ altitude"),
      dim = c(vertices = 20, focal = nrow(response), edge = 30),
      loglik = -10,
      aic = 26,
      df = 3,
      fit = list(response = response),
      submodels = list(g = model),
      comparison = terradish:::.terradish_comparison_contract(model, nu, response)
    )
  }

  S <- diag(4)
  expect_error(
    aic_table(list(make_fit(mlpe, S),
                   make_fit(generalized_wishart, S, nu = 20))),
    "different likelihood families"
  )
  expect_error(
    aic_table(list(make_fit(mlpe, S), make_fit(mlpe, S + 1))),
    "same response matrix"
  )
  expect_error(
    aic_table(list(make_fit(generalized_wishart, S, nu = 20),
                   make_fit(generalized_wishart, S, nu = 40))),
    "same effective degrees"
  )

  unknown_model <- function(...) NULL
  class(unknown_model) <- c("terradish_measurement_model",
                            "radish_measurement_model")
  expect_error(
    aic_table(list(make_fit(unknown_model, S), make_fit(unknown_model, S))),
    "Could not identify every model's likelihood family"
  )
})


test_that("aic_table keeps log-likelihoods aligned when model labels repeat", {
  response <- diag(4)
  make_fit <- function(loglik, aic) {
    list(
      formula = stats::as.formula("response ~ altitude"),
      dim = c(vertices = 20, focal = 4, edge = 30),
      loglik = loglik,
      aic = aic,
      df = 3,
      fit = list(response = response),
      submodels = list(g = mlpe),
      comparison = terradish:::.terradish_comparison_contract(mlpe)
    )
  }

  tab <- aic_table(list(make_fit(-20, 46), make_fit(-10, 26)),
                   mod_names = c("duplicate", "duplicate"))
  expect_equal(tab$loglik, c(-10, -20))
})
