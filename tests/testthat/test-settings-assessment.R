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
