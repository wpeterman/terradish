# Tier 3 directed (non-reversible) conductance engine. Uses small, well-conditioned
# synthetic lattices (the R per-absorber solver is exact but ill-conditions under
# strongly heterogeneous rates; full-resolution fitting is the C++ backend's job).

mk_directed_fixture <- function(DIM = 7L, directions = 8L) {
  r <- terra::rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM, ymin = 0, ymax = DIM)
  gx <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  gy <- terra::yFromCell(r, seq_len(terra::ncell(r)))
  v1   <- scale(gx + 0.5 * gy)[, 1]
  elev <- scale(gx)[, 1]
  covs <- c(terra::setValues(r, v1), terra::setValues(r, elev))
  names(covs) <- c("v1", "elev")
  fc <- c(1L, DIM, DIM * DIM, DIM * (DIM - 1L) + 1L, (DIM * DIM) %/% 2L, DIM * 3L + 2L)
  coords <- terra::xyFromCell(r, fc)
  surface <- conductance_surface(covs, coords, directions = directions)
  list(surface = surface, dir_cov = edge_gradient(covs[["elev"]], surface))
}

test_that("edge_gradient is antisymmetric over directed edges", {
  fx <- mk_directed_fixture(6L)
  ed <- fx$dir_cov$edges
  m <- nrow(ed) / 2
  # .directed_edges stacks forward then reverse, so edge k and k+m are reverses
  expect_equal(ed[seq_len(m), 1], ed[m + seq_len(m), 2])
  expect_equal(fx$dir_cov$d[seq_len(m)], -fx$dir_cov$d[m + seq_len(m)], tolerance = 1e-10)
})

test_that("gamma = 0 gives a reversible (symmetric-rate) generator", {
  fx <- mk_directed_fixture(6L)
  gen <- terradish:::.directed_generator(~ v1, fx$surface, fx$dir_cov)
  m <- nrow(gen$edges) / 2
  r0 <- gen$rates(c(0.4, 0))
  expect_equal(r0[seq_len(m)], r0[m + seq_len(m)], tolerance = 1e-10)
})

test_that("directed engine gradient matches numDeriv (generalized Wishart)", {
  skip_if_not_installed("numDeriv")
  fx <- mk_directed_fixture(7L)
  surface <- fx$surface; nf <- length(surface$demes)
  gen <- terradish:::.directed_generator(~ v1, surface, fx$dir_cov)
  nu <- 800
  E <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = c(0.5, 0.5))$covariance
  set.seed(7)
  Scov <- rWishart(1, df = nu, Sigma = (E + 0.2 * diag(nf)) / nu)[, , 1]
  S <- outer(diag(Scov), rep(1, nf)) + outer(rep(1, nf), diag(Scov)) - 2 * Scov
  diag(S) <- 0

  par <- c(0.3, 0.4)
  an <- terradish_directed_algorithm(gen, generalized_wishart, surface, S, par,
                                     nu = nu, gradient = TRUE)$gradient
  nd <- numDeriv::grad(function(pp)
    terradish_directed_algorithm(gen, generalized_wishart, surface, S, pp,
                                 nu = nu, gradient = FALSE)$objective, par)
  expect_equal(unname(an), nd, tolerance = 1e-4)
})

test_that("SparseLU directed backend matches Matrix reference", {
  fx <- mk_directed_fixture(7L)
  surface <- fx$surface; nf <- length(surface$demes)
  gen <- terradish:::.directed_generator(~ v1, surface, fx$dir_cov)
  par <- c(0.3, 0.45)

  ref <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = par,
                                      solver = "matrix")
  lu <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = par,
                                     solver = "sparse_lu_cpp")
  expect_equal(lu$covariance, ref$covariance, tolerance = 1e-8)
  expect_equal(lu$Hf, ref$Hf, tolerance = 1e-8)

  nu <- 800
  set.seed(17)
  Scov <- rWishart(1, df = nu,
                   Sigma = (ref$covariance + 0.2 * diag(nf)) / nu)[, , 1]
  S <- outer(diag(Scov), rep(1, nf)) + outer(rep(1, nf), diag(Scov)) - 2 * Scov
  diag(S) <- 0

  ref_fit <- terradish_directed_algorithm(gen, generalized_wishart, surface, S,
                                          par, nu = nu, gradient = TRUE,
                                          solver = "matrix")
  lu_fit <- terradish_directed_algorithm(gen, generalized_wishart, surface, S,
                                         par, nu = nu, gradient = TRUE,
                                         solver = "sparse_lu_cpp")
  expect_equal(lu_fit$objective, ref_fit$objective, tolerance = 1e-8)
  expect_equal(lu_fit$covariance, ref_fit$covariance, tolerance = 1e-8)
  expect_equal(unname(lu_fit$gradient), unname(ref_fit$gradient), tolerance = 1e-6)
})

test_that("directed covariance normalization avoids large-graph phi singularity", {
  fx <- mk_directed_fixture(60L)
  surface <- fx$surface; nf <- length(surface$demes)
  gen <- terradish:::.directed_generator(~ v1, surface, fx$dir_cov)
  par <- c(v1 = 4, gamma_elev = 4)

  E <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = par,
                                    solver = "sparse_lu_cpp")$covariance
  Escale <- max(abs(E))
  S <- E / Escale + 0.2 * diag(nf)
  ctrl <- NewtonRaphsonControl(verbose = FALSE, ftol = 1e-10, ctol = 1e-10)

  expect_gt(Escale, 1e12)
  expect_error(
    terradish:::radish_subproblem(wishart_covariance, E, S, nu = 1000,
                                  control = ctrl),
    "system is computationally singular",
    fixed = TRUE
  )

  scaled_subproblem <- terradish:::radish_subproblem(
    wishart_covariance, E / Escale, S, nu = 1000, control = ctrl)
  normalized <- terradish_directed_algorithm(
    gen, wishart_covariance, surface, S, par = par, nu = 1000,
    gradient = TRUE, solver = "sparse_lu_cpp")

  expect_true(is.finite(normalized$objective))
  expect_lt(abs(normalized$objective), 1e12)
  expect_true(all(is.finite(normalized$gradient)))
  expect_equal(normalized$objective, scaled_subproblem$loglikelihood,
               tolerance = 1e-8)
})

test_that("directed_rates and directed plots expose fitted edge bias", {
  fx <- mk_directed_fixture(6L)
  surface <- fx$surface
  theta <- c(v1 = 0.5)
  gamma <- c(gamma_elev = 0.6)
  vc <- diag(c(0.01, 0.02), 2)
  dimnames(vc) <- list(c(names(theta), names(gamma)),
                       c(names(theta), names(gamma)))
  fit <- list(
    formula = S ~ v1,
    theta = theta,
    gamma = gamma,
    vcov = vc,
    logconductance = surface$x$v1 * theta[["v1"]]
  )
  class(fit) <- c("terradish_directed", "terradish")

  rates <- directed_rates(fit, data = surface, directional = fx$dir_cov)
  m <- nrow(surface$edge_pairs)
  expect_s3_class(rates, "data.frame")
  expect_equal(nrow(rates), m)
  expect_true(all(rates$rate_ab > 0))
  expect_true(all(rates$rate_ba > 0))
  expect_equal(rates$log_rate_ratio, log(rates$rate_ab / rates$rate_ba),
               tolerance = 1e-10)
  expect_equal(rates$symmetric_rate, sqrt(rates$rate_ab * rates$rate_ba),
               tolerance = 1e-10)
  expect_true(all(rates$favored_from == rates$a | rates$favored_from == rates$b))
  expect_true(all(rates$favored_to == rates$a | rates$favored_to == rates$b))
  expect_true(all(is.finite(rates$abs_log_rate_ratio)))
  expect_true(all(is.finite(rates$log_rate_ratio_se)))

  p_directional <- plot(fit, type = "directional", data = surface,
                        directional = fx$dir_cov)
  p_combined <- plot(fit, type = "combined", data = surface,
                     directional = fx$dir_cov)
  expect_s3_class(p_directional, "ggplot")
  expect_s3_class(p_combined, "ggplot")
})

test_that("terradish_directed recovers symmetric and directional effects", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fx <- mk_directed_fixture(6L)
  surface <- fx$surface; nf <- length(surface$demes)
  gen <- terradish:::.directed_generator(~ v1, surface, fx$dir_cov)
  theta_true <- 0.5; gamma_true <- 0.6; nu <- 1000
  E <- terradish_directed_algorithm(gen, NULL, surface, NULL,
                                    par = c(theta_true, gamma_true))$covariance
  set.seed(11)
  Scov <- rWishart(1, df = nu, Sigma = (E + 0.2 * diag(nf)) / nu)[, , 1]
  S <- outer(diag(Scov), rep(1, nf)) + outer(rep(1, nf), diag(Scov)) - 2 * Scov
  diag(S) <- 0

  fit <- terradish_directed(S ~ v1, data = surface, directional = fx$dir_cov,
                            measurement_model = generalized_wishart, nu = nu,
                            estimate_vcov = FALSE,
                            solver = "sparse_lu_cpp")
  expect_s3_class(fit, "terradish_directed")
  expect_equal(unname(fit$theta), theta_true, tolerance = 0.35)
  expect_equal(unname(fit$gamma), gamma_true, tolerance = 0.35)
  expect_output(print(fit), "Directional")
  expect_length(coef(fit), 2L)
})

test_that("terradish_directed exposes batch-fit optimizer controls", {
  directed_args <- names(formals(terradish_directed))
  expect_true("estimate_vcov" %in% directed_args)
  expect_true("control" %in% directed_args)
  expect_true("solver" %in% directed_args)
})
