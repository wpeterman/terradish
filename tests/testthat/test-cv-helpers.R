test_that("terradish_cv returns a cross-validation summary", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst

  out <- suppressWarnings(
    terradish_cv(dat$coords, dat$covariates,
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

test_that("terradish_cv forwards supported conductance models", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst

  out <- suppressWarnings(
    terradish_cv(dat$coords, dat$covariates,
                 melip.Fst ~ s(altitude, df = 3),
                 model = "ls",
                 conductance_model = smooth_loglinear_conductance,
                 prop_train = 0.8,
                 seed = 1,
                 fit_full = FALSE,
                 control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_length(out$cv_loglik, 1L)
  expect_true(is.finite(out$cv_loglik))
  expect_true(grepl("^s\\(altitude\\)\\.", names(coef(out$train_mod))[1]))
})

test_that("terradish_cv rebuilds fixed-graph conductance factories on splits", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst

  fixed_graph_factory <- function(formula, x)
    stop("factory should be rebuilt for each CV surface", call. = FALSE)
  class(fixed_graph_factory) <- c("terradish_conductance_model_factory",
                                  "radish_conductance_model_factory")
  attr(fixed_graph_factory, "requires_fixed_graph") <- TRUE
  attr(fixed_graph_factory, "rebuild_for_surface") <- function(formula, surface,
                                                               reference_model = NULL)
    loglinear_conductance(formula, surface$x)

  out <- suppressWarnings(
    terradish_cv(dat$coords, dat$covariates,
                 melip.Fst ~ altitude + forestcover,
                 model = "ls",
                 conductance_model = fixed_graph_factory,
                 prop_train = 0.8,
                 seed = 1,
                 fit_full = FALSE,
                 control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_length(out$cv_loglik, 1L)
  expect_true(is.finite(out$cv_loglik))
})

test_that("terradish_cv_replicates summarizes repeated held-out loglikelihood", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst

  out <- suppressWarnings(
    terradish_cv_replicates(dat$coords, dat$covariates,
                            melip.Fst ~ altitude + forestcover,
                            model = "ls",
                            seeds = c(1L, 2L),
                            fit_full = FALSE,
                            keep_fits = TRUE,
                            control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_true(is.list(out))
  expect_s3_class(out, "terradish_cv_replicates")
  expect_named(out, c("summary", "mean_loglik", "sd_loglik", "seeds", "fits"))
  expect_true(is.data.frame(out$summary))
  expect_equal(nrow(out$summary), 2L)
  expect_true(all(c("replicate", "seed", "cv_loglik") %in% names(out$summary)))
  expect_length(out$fits, 2L)
  expect_gt(length(capture.output(print(out))), 0)

  sm <- summary(out)
  expect_s3_class(sm, "summary.terradish_cv_replicates")
  expect_equal(sm$replicates, 2L)
  expect_gt(length(capture.output(print(sm))), 0)
})

test_that("cv_model_selection and terradish_parameters summarize outputs", {
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

  sel <- cv_model_selection(list(
    list(train_mod = fit1, cv_loglik = fit1$loglik, full_mod = fit1),
    list(train_mod = fit2, cv_loglik = fit2$loglik, full_mod = fit2)
  ), aic = TRUE)

  expect_true(is.data.frame(sel$loglik_tab))
  expect_true(is.data.frame(sel$AIC_tab))

  tmp <- tempfile()
  dir.create(tmp)
  suppressWarnings(saveRDS(fit2, file.path(tmp, "fit--ls.rds")))
  saveRDS(list(effect_size = coef(fit2)), file.path(tmp, "AllResults_list.rds"))
  results_meta <- terradish_results(tmp)
  expect_true(is.list(results_meta))
  expect_true("all_results" %in% names(results_meta))
  expect_true("all_dirs" %in% names(results_meta))

  params <- terradish_parameters(tmp, model = "ls", save_table = FALSE)
  expect_true(is.data.frame(params))
  expect_equal(nrow(params), length(coef(fit2)))

  expect_warning(
    legacy_params <- terradish_parameters(tmp, radish_model = "ls", save_table = FALSE),
    "deprecated"
  )
  expect_equal(legacy_params$parameter, params$parameter)
})
