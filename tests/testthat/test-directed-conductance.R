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
                            measurement_model = generalized_wishart, nu = nu)
  expect_s3_class(fit, "terradish_directed")
  expect_equal(unname(fit$theta), theta_true, tolerance = 0.35)
  expect_equal(unname(fit$gamma), gamma_true, tolerance = 0.35)
  expect_output(print(fit), "Directional")
  expect_length(coef(fit), 2L)
})
