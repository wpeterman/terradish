test_that("radish_cv returns a cross-validation summary", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst

  out <- suppressWarnings(
    radish_cv(dat$coords, dat$covariates,
              melip.Fst ~ altitude + forestcover,
              model = "ls",
              prop_train = 0.8,
              seed = 1,
              fit_full = FALSE,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_true(is.list(out))
  expect_named(out, c("train_mod", "cv_loglik", "seed", "train_index", "test_index"))
  expect_length(out$cv_loglik, 1L)
})

test_that("cv_model_selection and radish_parameters summarize outputs", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst

  fit1 <- suppressWarnings(
    radish(melip.Fst ~ altitude,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )
  fit2 <- suppressWarnings(
    radish(melip.Fst ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  sel <- cv_model_selection(list(
    list(train_mod = fit1, cv_loglik = fit1$loglik, full_mod = fit1),
    list(train_mod = fit2, cv_loglik = fit2$loglik, full_mod = fit2)
  ), aic = TRUE)

  expect_true(is.data.frame(sel$loglik_tab))
  expect_true(is.data.frame(sel$AIC_tab))

  tmp <- tempfile()
  dir.create(tmp)
  saveRDS(fit2, file.path(tmp, "fit--ls.rds"))
  saveRDS(list(effect_size = coef(fit2)), file.path(tmp, "AllResults_list.rds"))
  params <- radish_parameters(tmp, radish_model = "ls", save_table = FALSE)
  expect_true(is.data.frame(params))
  expect_equal(nrow(params), length(coef(fit2)))
})
