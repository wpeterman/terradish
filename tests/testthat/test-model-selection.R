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

  expect_error(aic_table(list(fit1, fit2), AICc = TRUE, BIC = TRUE),
               "Set only one")
})

test_that("cv_model_selection forwards AICc and BIC to aic_table", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  fit1 <- fx$fit
  melip.Fst <- fx$data$melip.Fst

  fit2 <- suppressWarnings(
    terradish(melip.Fst ~ altitude,
              data = fx$surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  cv_list <- list(
    list(train_mod = fit1, cv_loglik = fit1$loglik, full_mod = fit1),
    list(train_mod = fit2, cv_loglik = fit2$loglik, full_mod = fit2)
  )

  aicc_sel <- cv_model_selection(cv_list, aic = TRUE, AICc = TRUE)
  expect_true(is.list(aicc_sel))
  expect_true("AICc" %in% names(aicc_sel$AIC_tab))

  bic_sel <- cv_model_selection(cv_list, aic = TRUE, BIC = TRUE)
  expect_true(is.list(bic_sel))
  expect_true("BIC" %in% names(bic_sel$AIC_tab))
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

  cv_list <- list(
    list(train_mod = fit_base, cv_loglik = fit_base$loglik, full_mod = fit_base),
    list(train_mod = fit_joint, cv_loglik = fit_joint$loglik, full_mod = fit_joint)
  )
  cv_out <- cv_model_selection(cv_list, aic = TRUE)
  expect_true(any(cv_out$loglik_tab$model == "altitude + forestcover [mlpe:3]"))
  expect_true(any(cv_out$AIC_tab$model == "altitude + forestcover [mlpe:3]"))
})

test_that("aic_table rejects incomparable likelihoods, responses, pairs, and nu", {
  make_fit <- function(model, response, nu = NULL, pairs = NULL) {
    if (!is.null(pairs))
      model <- pair_subset_measurement_model(model, pairs)

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
  expect_error(
    aic_table(list(make_fit(mlpe, S, pairs = rbind(c(1, 2), c(2, 3))),
                   make_fit(mlpe, S, pairs = rbind(c(1, 2), c(3, 4))))),
    "same selected pairs"
  )

  unknown_model <- function(...) NULL
  class(unknown_model) <- c("terradish_measurement_model",
                            "radish_measurement_model")
  expect_error(
    aic_table(list(make_fit(unknown_model, S), make_fit(unknown_model, S))),
    "Could not identify every model's likelihood family"
  )
})

test_that("aic_table uses selected pair rows as the pair-subset BIC convention", {
  response <- diag(4)
  pairs <- rbind(c(1, 2), c(2, 3))
  pair_mlpe <- pair_subset_measurement_model(mlpe, pairs)
  make_fit <- function(loglik, df) {
    list(
      formula = stats::as.formula("response ~ altitude"),
      dim = c(vertices = 20, focal = 4, edge = 30),
      loglik = loglik,
      aic = -2 * loglik + 2 * df,
      df = df,
      fit = list(response = response),
      submodels = list(g = pair_mlpe),
      comparison = terradish:::.terradish_comparison_contract(
        pair_mlpe, response = response
      )
    )
  }

  tab <- aic_table(list(make_fit(-10, 3), make_fit(-12, 4)),
                   BIC = TRUE, mod_names = c("fit 1", "fit 2"))
  expect_equal(tab$BIC[match("fit 1", tab$model)],
               round(20 + 3 * log(nrow(pairs)), 4))
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
