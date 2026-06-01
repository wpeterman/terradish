test_that("terradish_assess_settings profiles a graph and benchmarks direct solver settings", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords,
                                 directions = 8, saveStack = TRUE)

  assessment <- terradish_assess_settings(
    melip.Fst ~ altitude + forestcover,
    data = surface,
    measurement_model = leastsquares,
    optimizer_probe = FALSE,
    solver_probe = TRUE,
    factorization = "simplicial_ll",
    solve_backend = "matrix",
    verbose = FALSE
  )

  expect_s3_class(assessment, "terradish_setting_assessment")
  expect_equal(assessment$profile$n_vertices, nrow(surface$x))
  expect_equal(assessment$profile$n_focal, length(surface$demes))
  expect_equal(assessment$profile$n_parameters, 2L)
  expect_true(isTRUE(assessment$profile$can_coarse_raster))
  expect_equal(assessment$recommended$solver, "direct")
  expect_equal(assessment$recommended$solver_control$factorization,
               "simplicial_ll")
  expect_equal(assessment$recommended$solver_control$solve_backend,
               "matrix")
  expect_null(assessment$benchmarks$optimizer)
  expect_false(inherits(assessment$recommended$control$ls.control,
                        "terradish_armijo_control"))
  expect_output(print(assessment), "Recommended settings")
})

test_that("terradish_assess_settings can run an Armijo BFGS optimizer probe", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  assessment <- terradish_assess_settings(
    melip.Fst ~ altitude + forestcover,
    data = surface,
    measurement_model = leastsquares,
    optimizer_probe = TRUE,
    solver_probe = FALSE,
    optimizer_candidates = "bfgs_armijo",
    probe_maxit = 1L,
    verbose = FALSE
  )

  expect_equal(assessment$recommended$optimizer, "bfgs")
  expect_s3_class(assessment$recommended$control$ls.control,
                  "terradish_armijo_control")
  expect_equal(nrow(assessment$benchmarks$optimizer), 1L)
  expect_equal(assessment$benchmarks$optimizer$line_search, "armijo")
  expect_equal(assessment$benchmarks$optimizer$status, "OK")
  expect_gt(assessment$benchmarks$optimizer$line_search_objective_only_trials,
            0)
})

test_that("terradish_assess_settings reports smooth-model settings context", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  assessment <- terradish_assess_settings(
    melip.Fst ~ forestcover + s(altitude, df = 3),
    data = surface,
    conductance_model = smooth_loglinear_conductance,
    measurement_model = leastsquares,
    optimizer_probe = TRUE,
    solver_probe = FALSE,
    probe_maxit = 1L,
    verbose = FALSE
  )

  expect_true(isTRUE(assessment$profile$smooth_conductance))
  expect_equal(assessment$profile$n_smooth_terms, 1L)
  expect_equal(assessment$profile$n_smooth_basis_columns, 3L)
  expect_equal(assessment$defaults$optimizer, "bfgs")
  expect_true(grepl("defaults", assessment$comparison$summary, fixed = TRUE))
  expect_true(any(grepl("Smooth conductance formula expanded to",
                        assessment$notes,
                        fixed = TRUE)))
  expect_true(all(assessment$benchmarks$optimizer$status == "OK"))
})

test_that("settings comparison reports tuned diagnostic controls", {
  defaults <- list(
    optimizer = "newton",
    control = terradish:::.terradish_assessment_control("hager_zhang"),
    solver = "direct",
    solver_control = NULL,
    approximation = "none",
    approximation_control = NULL
  )
  recommended <- defaults
  recommended$solver_control <- list(
    factorization = "simplicial_ll",
    solve_backend = "matrix"
  )
  recommended$approximation_control <- list(factor = 2L)

  comparison <- terradish:::.terradish_assessment_compare_settings(
    defaults,
    recommended
  )

  expect_true(comparison$differs_from_defaults)
  expect_true("solver_control tuned from terradish defaults" %in%
                comparison$changes)
  expect_true("approximation_control tuned from terradish defaults" %in%
                comparison$changes)
})
