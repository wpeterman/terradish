# DRAGON (structured-coalescent directed engine) regression tests.
#
# The reference values below were produced by the independently validated Python
# engine (DRAGON/derisk_coalescent_vs_commute/dragon_engine.py + dragon_couple.py,
# adjoint checked vs finite differences to ~1e-7) and reproduced by a column-major
# re-implementation to ~1e-11. They pin the R port's forward map and analytic
# adjoint gradient for the uniform, drift, and FRAME-coupled coalescence models.

dragon_test_setup <- function() {
  DIM <- 5L; d <- 25L
  elev <- rep(0:4, times = DIM)                 # node directional potential (col index)
  gy   <- rep(0:4, each = DIM)
  z    <- (gy - mean(gy)) / sqrt(mean((gy - mean(gy))^2))   # population-sd standardized
  ep <- do.call(rbind, lapply(seq_len(d), function(k) {
    i <- k - 1L; cx <- i %% DIM; cy <- i %/% DIM
    rbind(if (cx < DIM - 1L) c(k, k + 1L),
          if (cy < DIM - 1L) c(k, k + DIM))
  }))
  Shat <- as.matrix(utils::read.csv(
    test_path("dragon_fixture_Shat_node.csv"), header = FALSE))
  dimnames(Shat) <- NULL
  list(ep = ep, elev = elev, z = z, Shat = Shat, nu = 600)
}

test_that("DRAGON forward + adjoint match the validated reference fixture", {
  s <- dragon_test_setup()

  m_u <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, mode = "uniform")
  g_u <- .dragon_grad(c(0.5, -0.3, -1.2, -3.5), m_u)
  expect_equal(g_u$objective, 25919.5022, tolerance = 1e-4)
  expect_equal(g_u$gradient,
               c(37939.978, -2188.060, -22451.529, -3617.734), tolerance = 1e-3)

  m_d <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, z = s$z, mode = "drift")
  g_d <- .dragon_grad(c(0.5, -0.3, 0.4, -1.2, -3.5), m_d)
  expect_equal(g_d$objective, 26932.4521, tolerance = 1e-4)
  expect_equal(g_d$gradient,
               c(38309.051, -3033.566, 6070.929, -23110.573, -4106.747),
               tolerance = 1e-3)

  m_c <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, mode = "coupled")
  g_c <- .dragon_grad(c(0.5, -0.3, 0.6, -1.2, -3.5), m_c)
  expect_equal(g_c$objective, 16829.4295, tolerance = 1e-4)
  expect_equal(g_c$gradient,
               c(7355.642, 853.424, -1843.282, -12131.232, -1315.064),
               tolerance = 1e-3)
})

test_that("DRAGON analytic gradient matches numDeriv (all coalescence modes)", {
  skip_if_not_installed("numDeriv")
  s <- dragon_test_setup()
  cases <- list(
    list(mode = "uniform", par = c(0.4, -0.2, -1.3, -3.6), z = NULL),
    list(mode = "drift",   par = c(0.4, -0.2, 0.3, -1.3, -3.6), z = s$z),
    list(mode = "coupled", par = c(0.4, -0.2, 0.5, -1.3, -3.6), z = NULL))
  for (cs in cases) {
    m <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, z = cs$z, mode = cs$mode)
    ga <- .dragon_grad(cs$par, m)$gradient
    gn <- numDeriv::grad(function(p) .dragon_forward(p, m)$objective, cs$par)
    expect_equal(ga, gn, tolerance = 1e-4,
                 info = paste("mode =", cs$mode))
  }
})

test_that("coupled reduces to uniform at alpha = 0", {
  s <- dragon_test_setup()
  m_u <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, mode = "uniform")
  m_c <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, mode = "coupled")
  o_u <- .dragon_forward(c(0.5, -0.3, -1.2, -3.5), m_u)$objective
  o_c <- .dragon_forward(c(0.5, -0.3, 0.0, -1.2, -3.5), m_c)$objective
  expect_equal(o_u, o_c, tolerance = 1e-6)
})

test_that("dragon() fits, prints, and exposes coef/logLik/AIC", {
  s <- dragon_test_setup()
  data <- list(edge_pairs = s$ep)
  class(data) <- "terradish_graph"
  fit <- dragon(directional = s$elev, data = data, S = s$Shat, nu = s$nu,
                coalescence = "uniform", n_start = 4L)
  expect_s3_class(fit, "dragon")
  expect_true(is.finite(fit$loglik))
  expect_true("gamma_dir" %in% names(coef(fit)))
  expect_true(is.finite(AIC(fit)))
})

test_that("dragon() handles a circulation (antisymmetric edge) covariate", {
  skip_if_not_installed("numDeriv")
  s <- dragon_test_setup()
  set.seed(1)
  circ <- runif(nrow(s$ep), -1, 1)                 # one value per undirected edge
  dedge <- matrix(c(circ, -circ), ncol = 1L)        # directed: forward +, reverse -

  # the analytic adjoint must be correct for the extra directional column (q = 2)
  m <- .dragon_model(s$ep, s$elev, s$Shat, nu = s$nu, dedge = dedge, mode = "coupled")
  par <- c(0.4, 0.3, -0.2, 0.5, -1.3, -3.6)         # gamma_dir, gamma_circ, c0, alpha, ltau, lnug
  ga <- .dragon_grad(par, m)$gradient
  gn <- numDeriv::grad(function(p) .dragon_forward(p, m)$objective, par)
  expect_equal(ga, gn, tolerance = 1e-4)

  # dragon() accepts `circulation`, names coefficients, and works potential-free
  data <- list(edge_pairs = s$ep); class(data) <- "terradish_graph"
  fit <- dragon(directional = s$elev, data = data, S = s$Shat, nu = s$nu,
                coalescence = "uniform", circulation = circ, n_start = 3L)
  expect_true(all(c("gamma_dir", "gamma_circ") %in% names(coef(fit))))

  fit0 <- dragon(directional = NULL, data = data, S = s$Shat, nu = s$nu,
                 coalescence = "uniform", circulation = circ, n_start = 3L)
  expect_true("gamma_circ" %in% names(coef(fit0)))
  expect_false("gamma_dir" %in% names(coef(fit0)))

  expect_error(dragon(directional = NULL, data = data, S = s$Shat, nu = s$nu),
               "directional.*circulation")
})

test_that("edge_flow builds an antisymmetric circulation covariate for dragon()", {
  s <- dragon_test_setup()
  gx <- s$elev; gy <- rep(0:4, each = 5L)
  vc <- cbind(gx, gy)
  data <- list(edge_pairs = s$ep, vertex_coordinates = vc)
  class(data) <- "terradish_graph"

  cen <- colMeans(vc)
  rot <- function(xy) cbind(-(xy[, 2] - cen[2]), xy[, 1] - cen[1])  # counter-clockwise curl
  circ <- edge_flow(rot, data)
  expect_length(circ, nrow(s$ep))
  expect_true(any(abs(circ) > 1e-8))                                # non-trivial flow

  # reversing each edge's orientation negates the covariate (antisymmetry)
  data_rev <- data; data_rev$edge_pairs <- s$ep[, 2:1, drop = FALSE]
  expect_equal(edge_flow(rot, data_rev), -circ)

  # the matrix and function forms agree, and the result drives dragon()
  expect_equal(edge_flow(rot(vc), data), circ)
  fit <- dragon(directional = s$elev, data = data, S = s$Shat, nu = s$nu,
                coalescence = "uniform", circulation = circ, n_start = 3L)
  expect_true("gamma_circ" %in% names(coef(fit)))
})
