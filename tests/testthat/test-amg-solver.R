test_that("AMG solver matches direct solver on melip subproblem", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  fst <- ifelse(melip.Fst < 0, 0, melip.Fst)

  fit_direct <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    theta = c(-0.3, 0.3),
    partial = FALSE,
    solver = "direct"
  )

  fit_amg <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    theta = c(-0.3, 0.3),
    partial = FALSE,
    solver = "amg",
    solver_control = list(
      tol = 1e-8,
      maxit = 400L,
      coarse_enough = 1000L,
      estimate_spectral_radius = TRUE,
      power_iters = 4L
    )
  )

  expect_equal(fit_amg$solver_info$type, "amg")
  expect_true(all(fit_amg$solver_info$converged))
  expect_equal(fit_direct$objective, fit_amg$objective, tolerance = 1e-6)
  expect_equal(as.matrix(fit_direct$covariance), as.matrix(fit_amg$covariance), tolerance = 1e-6)
  expect_equal(fit_direct$gradient, fit_amg$gradient, tolerance = 1e-5)
  expect_equal(fit_direct$hessian, fit_amg$hessian, tolerance = 1e-4)
})

test_that("compiled Laplacian derivative products match sparse Matrix products", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  C <- model(c(altitude = -0.2, forestcover = 0.3))
  dconductance <- C$df__dtheta(1L)

  solver_state <- terradish:::.terradish_solver_setup(
    surface,
    C$conductance,
    solver = "direct"
  )
  G <- as.matrix(terradish:::.terradish_solver_solve(
    solver_state,
    surface$rhs
  )$solution)

  dQn_sparse <- Matrix::forceSymmetric(
    backpropagate_conductance_to_laplacian(dconductance, surface$adj)
  )
  dQnG_sparse <- as.matrix(dQn_sparse %*% G)
  dQnG_compiled <- laplacian_derivative_matrix_product(
    dconductance,
    surface$adj,
    G
  )

  expect_equal(dQnG_compiled, dQnG_sparse, tolerance = 1e-12)
})

test_that("compiled reduced-RHS products match explicit matrix multiplication", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  n_vertices <- nrow(surface$x)
  Zn <- terradish:::.graph_rhs(surface, n_vertices)

  left_input <- matrix(rnorm(ncol(Zn) * 5L), nrow = ncol(Zn), ncol = 5L)
  expect_equal(
    graph_rhs_matrix_product(surface$demes, n_vertices, left_input),
    as.matrix(Zn %*% left_input),
    tolerance = 1e-12
  )

  right_input <- matrix(rnorm(nrow(Zn) * 4L), nrow = nrow(Zn), ncol = 4L)
  expect_equal(
    graph_rhs_crossprod(surface$demes, n_vertices, right_input),
    as.matrix(crossprod(Zn, right_input)),
    tolerance = 1e-12
  )
})

test_that("adaptive AMG schedule stages early and final tolerances", {
  early <- terradish:::.terradish_solver_control_for_phase(
    solver = "amg",
    solver_control = list(
      adaptive = TRUE,
      tol_early = 1e-3,
      tol_mid = 1e-5,
      tol_final = 1e-8,
      maxit_early = 25L,
      maxit_mid = 50L,
      maxit_final = 100L,
      warmup_evals = 2L
    ),
    eval_count = 1L,
    final = FALSE
  )

  late <- terradish:::.terradish_solver_control_for_phase(
    solver = "amg",
    solver_control = list(
      adaptive = TRUE,
      tol_early = 1e-3,
      tol_mid = 1e-5,
      tol_final = 1e-8,
      maxit_early = 25L,
      maxit_mid = 50L,
      maxit_final = 100L,
      warmup_evals = 2L
    ),
    eval_count = 5L,
    final = FALSE
  )

  final <- terradish:::.terradish_solver_control_for_phase(
    solver = "amg",
    solver_control = list(
      adaptive = TRUE,
      tol_early = 1e-3,
      tol_mid = 1e-5,
      tol_final = 1e-8,
      maxit_early = 25L,
      maxit_mid = 50L,
      maxit_final = 100L,
      warmup_evals = 2L
    ),
    eval_count = 5L,
    final = TRUE
  )

  expect_equal(early$adaptive_phase, "early")
  expect_equal(early$tol, 1e-3)
  expect_equal(early$maxit, 25L)

  expect_equal(late$adaptive_phase, "mid")
  expect_equal(late$tol, 1e-5)
  expect_equal(late$maxit, 50L)

  expect_equal(final$adaptive_phase, "final")
  expect_equal(final$tol, 1e-8)
  expect_equal(final$maxit, 100L)
})

test_that("auto solver resolves conservatively on the bundled melip graph", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  resolution <- terradish:::.terradish_resolve_solver(
    s = surface,
    solver = "auto",
    n_vertices = nrow(surface$x)
  )

  expect_equal(resolution$type, "direct")
  expect_equal(resolution$requested_type, "auto")
  expect_equal(resolution$reason, "graph_not_large_enough_for_amg")
})

test_that("auto solver keeps moderately large graphs on the direct backend by default", {
  resolution <- terradish:::.terradish_resolve_solver(
    s = list(rhs = matrix(0, nrow = 25, ncol = 25)),
    solver = "auto",
    n_vertices = 1000000L
  )

  expect_equal(resolution$type, "direct")
  expect_equal(resolution$requested_type, "auto")
  expect_equal(resolution$reason, "prefer_direct_until_larger_graphs")
})

test_that("direct factorization chooser switches to supernodal only for large graphs", {
  default_control <- terradish:::.terradish_direct_control_defaults()

  expect_equal(
    terradish:::.terradish_resolve_direct_factorization(default_control, n_vertices = 22443L, n_rhs = 37L),
    "simplicial_ldl"
  )

  expect_equal(
    terradish:::.terradish_resolve_direct_factorization(default_control, n_vertices = 62500L, n_rhs = 25L),
    "supernodal_ll"
  )
})

test_that("terradish final AMG fit uses final tolerance", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE),
      solver = "amg",
      solver_control = list(
        adaptive = TRUE,
        tol_early = 1e-4,
        tol_mid = 1e-5,
        tol_final = 1e-7,
        maxit_early = 25L,
        maxit_mid = 50L,
        maxit_final = 75L,
        warmup_evals = 1L
      )
    )
  )

  expect_equal(fit$fit$solver_info$type, "amg")
  expect_equal(fit$fit$solver_info$adaptive_phase, "final")
  expect_equal(fit$fit$solver_info$target_tol, 1e-7)
  expect_equal(fit$fit$solver_info$target_maxit, 75L)
})

test_that("auto solver reports the resolved backend in terradish fits", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit_direct <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE),
      solver = "auto"
    )
  )

  fit_amg <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE),
      solver = "auto",
      solver_control = list(
        auto_direct_max_vertices = 1L,
        auto_amg_min_vertices = 2L,
        adaptive = FALSE,
        tol = 1e-8,
        maxit = 400L,
        coarse_enough = 1000L,
        estimate_spectral_radius = TRUE,
        power_iters = 4L
      )
    )
  )

  expect_equal(fit_direct$fit$solver_info$requested_type, "auto")
  expect_equal(fit_direct$fit$solver_info$type, "direct")
  expect_equal(fit_direct$fit$solver_info$auto_reason, "graph_not_large_enough_for_amg")

  expect_equal(fit_amg$fit$solver_info$requested_type, "auto")
  expect_equal(fit_amg$fit$solver_info$type, "amg")
  expect_true(grepl("amg", fit_amg$fit$solver_info$auto_reason))
  expect_equal(fit_direct$loglik, fit_amg$loglik, tolerance = 1e-5)
})

test_that("direct solver reports factorization mode and can reuse templates", {
  dat <- melip_fixture(1:10)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  fst <- ifelse(melip.Fst < 0, 0, melip.Fst)

  fit1 <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    theta = c(-0.3, 0.3),
    partial = FALSE,
    solver = "direct"
  )

  fit2 <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    phi = fit1$phi,
    theta = c(-0.28, 0.32),
    partial = FALSE,
    solver = "direct",
    solver_reuse_state = fit1$solver_reuse_state
  )

  fit_super <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    theta = c(-0.3, 0.3),
    partial = FALSE,
    solver = "direct",
    solver_control = list(
      factorization = "auto",
      supernodal_min_vertices = 1L,
      supernodal_max_rhs = 100L
    )
  )

  expect_equal(fit1$solver_info$type, "direct")
  expect_equal(fit1$solver_info$factorization, "simplicial_ldl")
  expect_false(isTRUE(fit1$solver_info$reused_factor_template))
  expect_false(is.null(fit1$solver_reuse_state))

  expect_equal(fit2$solver_info$type, "direct")
  expect_true(isTRUE(fit2$solver_info$reused_factor_template))

  expect_equal(fit_super$solver_info$type, "direct")
  expect_equal(fit_super$solver_info$factorization, "supernodal_ll")
})

test_that("AMG hierarchy reuse is threaded across nearby evaluations", {
  dat <- melip_fixture(1:10)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  fst <- ifelse(melip.Fst < 0, 0, melip.Fst)

  control <- list(
    adaptive = FALSE,
    tol = 1e-8,
    maxit = 400L,
    coarse_enough = 1000L,
    estimate_spectral_radius = TRUE,
    power_iters = 4L,
    reuse_preconditioner = TRUE
  )

  fit1 <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    theta = c(-0.3, 0.3),
    partial = FALSE,
    solver = "amg",
    solver_control = control
  )

  fit2 <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    fst,
    nu = 1000,
    phi = fit1$phi,
    theta = c(-0.28, 0.32),
    partial = FALSE,
    solver = "amg",
    solver_control = control,
    solver_warm_start = fit1$solver_warm_start,
    solver_reuse_state = fit1$solver_reuse_state
  )

  expect_equal(fit1$solver_info$type, "amg")
  expect_false(isTRUE(fit1$solver_info$reused_preconditioner))
  expect_false(is.null(fit1$solver_reuse_state))

  expect_equal(fit2$solver_info$type, "amg")
  expect_true(isTRUE(fit2$solver_info$reused_preconditioner))
  expect_true(fit2$solver_info$reuse_age >= 1L)
  expect_false(is.null(fit2$solver_reuse_state))
})

test_that("terradish_solver_benchmark reports direct solver timings", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  bench <- terradish_solver_benchmark(
    fx$surface,
    factorization = c("simplicial_ldl", "simplicial_ll"),
    n_replicates = 1
  )

  expect_s3_class(bench, "data.frame")
  expect_equal(nrow(bench), 2L)
  expect_true(all(c("solver", "factorization", "setup_time",
                    "solve_time", "total_time") %in% names(bench)))
  expect_true(all(bench$solver == "direct"))
  expect_true(all(bench$total_time >= 0))
})
