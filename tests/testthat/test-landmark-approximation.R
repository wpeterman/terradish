test_that("landmark approximation subsets focal populations and rhs columns", {
  dat <- melip_fixture(1:10)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  subset <- terradish:::.terradish_landmark_subset(
    surface,
    dat$melip.Fst,
    approximation_control = list(
      n_landmarks = 4L,
      method = "spacefill",
      seed = 1
    )
  )

  expect_true(subset$used)
  expect_length(subset$index, 4L)
  expect_equal(nrow(subset$S), 4L)
  expect_equal(ncol(subset$S), 4L)
  expect_equal(length(subset$data$demes), 4L)
  expect_equal(ncol(subset$data$rhs), 4L)
})

test_that("terradish landmark approximation refines back to the full-data fit for mlpe", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  control <- NewtonRaphsonControl(maxit = 4, verbose = FALSE)

  fit_exact <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct"
    )
  )

  fit_landmark <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct",
      approximation = "landmark",
      approximation_control = list(
        n_landmarks = 4L,
        method = "spacefill",
        seed = 1,
        exact_refine = TRUE
      )
    )
  )

  expect_true(isTRUE(fit_landmark$approximation$used))
  expect_equal(fit_landmark$approximation$type, "landmark")
  expect_equal(fit_landmark$approximation$n_landmarks, 4L)
  expect_equal(dim(fit_landmark$fit$covariance), dim(fit_exact$fit$covariance))
  expect_equal(fit_landmark$loglik, fit_exact$loglik, tolerance = 1e-4)
  expect_lt(max(abs(fit_landmark$mle$theta - fit_exact$mle$theta)), 2e-2)
})

test_that("landmark exact refinement is guarded for leastsquares", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  control <- NewtonRaphsonControl(maxit = 4, verbose = FALSE)

  fit_landmark <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      control = control,
      solver = "direct",
      approximation = "landmark",
      approximation_control = list(
        n_landmarks = 4L,
        method = "spacefill",
        seed = 1,
        exact_refine = TRUE
      )
    )
  )

  expect_equal(fit_landmark$approximation$type, "landmark")
  expect_false(isTRUE(fit_landmark$approximation$used))
  expect_equal(fit_landmark$approximation$refine_guard, "disabled_for_leastsquares")
})

test_that("coarse-raster approximation refines back to the full-data fit", {
  dat <- melip_fixture(1:10)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  control <- NewtonRaphsonControl(maxit = 4, verbose = FALSE)

  fit_exact <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct"
    )
  )

  fit_coarse <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct",
      approximation = "coarse_raster",
      approximation_control = list(
        factor = 2L,
        exact_refine = TRUE
      )
    )
  )

  expect_equal(fit_coarse$approximation$type, "coarse_raster")
  expect_true(isTRUE(fit_coarse$approximation$used))
  expect_equal(fit_coarse$approximation$factor, 2L)
  expect_lt(fit_coarse$approximation$coarse_vertices,
            fit_coarse$approximation$full_vertices)
  expect_equal(fit_coarse$approximation$refine_control, "shared")
  expect_gte(fit_coarse$approximation$coarse_steps, 1L)
  expect_gte(fit_coarse$approximation$refine_steps, 1L)
  expect_equal(length(fit_coarse$approximation$stages), 1L)
  expect_equal(dim(fit_coarse$fit$covariance), dim(fit_exact$fit$covariance))
  expect_equal(fit_coarse$loglik, fit_exact$loglik, tolerance = 1e-4)
  expect_lt(max(abs(fit_coarse$mle$theta - fit_exact$mle$theta)), 2e-2)
})

test_that("coarse-raster approximation supports multilevel factor schedules", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit_coarse <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
      solver = "direct",
      approximation = "coarse_raster",
      approximation_control = list(
        factor = c(4L, 2L),
        exact_refine = TRUE,
        refine_control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
      )
    )
  )

  expect_equal(fit_coarse$approximation$type, "coarse_raster")
  expect_true(isTRUE(fit_coarse$approximation$used))
  expect_equal(fit_coarse$approximation$factor, c(4L, 2L))
  expect_equal(fit_coarse$approximation$n_levels, 2L)
  expect_equal(length(fit_coarse$approximation$coarse_vertices), 2L)
  expect_true(all(fit_coarse$approximation$coarse_vertices <
                    fit_coarse$approximation$full_vertices))
  expect_equal(fit_coarse$approximation$refine_control, "custom")
  expect_equal(fit_coarse$approximation$refine_steps, 1L)
  expect_equal(length(fit_coarse$approximation$stages), 2L)
  expect_true(all(vapply(fit_coarse$approximation$stages,
                         function(stage) stage$optimizer_steps >= 1L,
                         logical(1))))
})

test_that("auto optimizer resolves to BFGS above three parameters and BFGS reports steps", {
  expect_equal(terradish:::.resolve_terradish_optimizer("auto", 2L), "newton")
  expect_equal(terradish:::.resolve_terradish_optimizer("auto", 4L), "bfgs")

  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit_bfgs <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      optimizer = "bfgs",
      control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
      solver = "direct"
    )
  )

  expect_true("newton_steps" %in% names(fit_bfgs$cost))
  expect_gte(as.integer(fit_bfgs$cost[["newton_steps"]]), 1L)
  expect_false(is.null(fit_bfgs$fit$hessian))
  expect_false(is.null(fit_bfgs$mle$hessian))
  expect_equal(dim(fit_bfgs$mle$hessian),
               rep(length(fit_bfgs$mle$theta), 2L))
  expect_true(all(c("objective_evaluations", "gradient_evaluations",
                    "hessian_evaluations", "solver_setups",
                    "solver_solves", "line_search_trials") %in%
                    names(fit_bfgs$diagnostics)))
  expect_gte(fit_bfgs$diagnostics$objective_evaluations,
             fit_bfgs$diagnostics$gradient_evaluations)
  expect_gte(fit_bfgs$diagnostics$gradient_evaluations,
             fit_bfgs$diagnostics$hessian_evaluations)
  expect_gte(fit_bfgs$diagnostics$solver_setups,
             fit_bfgs$diagnostics$algorithm_calls)
  expect_gte(fit_bfgs$diagnostics$solver_solves,
             fit_bfgs$diagnostics$solver_setups)
  expect_gt(fit_bfgs$diagnostics$line_search_trials, 0L)
  expect_gt(fit_bfgs$diagnostics$line_search_gradient_trials, 0L)
})
