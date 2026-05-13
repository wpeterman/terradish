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
