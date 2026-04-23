test_that("pair-subset least squares reproduces full-pair least squares", {
  set.seed(11)
  n <- 5L
  X <- matrix(rnorm(n * n), n, n)
  E <- crossprod(X) + diag(n)
  S <- matrix(0, n, n)
  S[lower.tri(S)] <- seq_len(sum(lower.tri(S))) / 20
  S <- S + t(S)
  pairs <- which(lower.tri(S), arr.ind = TRUE)
  phi <- c(alpha = 0.1, beta = 0.7, tau = -0.2)

  pair_model <- pair_subset_measurement_model(leastsquares, pairs)
  base <- leastsquares(E, S, phi = phi, partial = TRUE)
  pair_fit <- pair_model(E, S, phi = phi, partial = TRUE)
  probe <- matrix(rnorm(n * n), n, n)
  probe <- (probe + t(probe)) / 2

  expect_s3_class(pair_model, "terradish_measurement_model")
  expect_equal(pair_fit$objective, base$objective)
  expect_equal(pair_fit$gradient, base$gradient)
  expect_equal(pair_fit$hessian, base$hessian)
  expect_equal(pair_fit$gradient_E, base$gradient_E)
  expect_equal(pair_fit$partial_E, base$partial_E, ignore_attr = TRUE)
  expect_equal(pair_fit$jacobian_E(probe), base$jacobian_E(probe))
  expect_equal(pair_fit$jacobian_S(probe), base$jacobian_S(probe))
})

test_that("pair-subset MLPE reproduces full-pair MLPE", {
  set.seed(12)
  n <- 5L
  X <- matrix(rnorm(n * n), n, n)
  E <- crossprod(X) + diag(n)
  S <- matrix(0, n, n)
  S[lower.tri(S)] <- seq_len(sum(lower.tri(S))) / 30
  S <- S + t(S)
  pairs <- which(lower.tri(S), arr.ind = TRUE)
  phi <- c(alpha = 0.1, beta = 0.7, tau = -0.2, rho = qlogis(0.2))

  pair_model <- pair_subset_measurement_model(mlpe, pairs)
  base <- mlpe(E, S, phi = phi, partial = TRUE)
  pair_fit <- pair_model(E, S, phi = phi, partial = TRUE)
  probe <- matrix(rnorm(n * n), n, n)
  probe <- (probe + t(probe)) / 2

  expect_equal(pair_fit$objective, base$objective)
  expect_equal(pair_fit$gradient, base$gradient)
  expect_equal(pair_fit$hessian, base$hessian)
  expect_equal(pair_fit$gradient_E, base$gradient_E)
  expect_equal(pair_fit$partial_E, base$partial_E, ignore_attr = TRUE)
  expect_equal(pair_fit$jacobian_E(probe), base$jacobian_E(probe))
  expect_equal(pair_fit$jacobian_S(probe), base$jacobian_S(probe))
})

test_that("selected pairs retain all sites but use only requested observations", {
  set.seed(13)
  n <- 5L
  X <- matrix(rnorm(n * n), n, n)
  E <- crossprod(X) + diag(n)
  S <- matrix(0, n, n)
  S[lower.tri(S)] <- seq_len(sum(lower.tri(S))) / 25
  S <- S + t(S)
  rownames(S) <- colnames(S) <- paste0("site", seq_len(n))
  selected <- rbind(c("site1", "site2"),
                    c("site1", "site4"),
                    c("site3", "site5"))

  pair_model <- pair_subset_measurement_model(leastsquares, selected)
  pair_fit <- pair_model(E, S, phi = c(0.1, 0.8, -0.1), partial = TRUE)

  expect_equal(dim(pair_fit$fitted), c(n, n))
  expect_equal(dim(pair_fit$gradient_E), c(n, n))
  expect_equal(dim(pair_fit$partial_S), c(3L, n * n))
  expect_equal(pair_fit$jacobian_S(matrix(1, n, n))[2, 3], 0)
  expect_false(isTRUE(all.equal(pair_fit$objective,
                                leastsquares(E, S, phi = c(0.1, 0.8, -0.1))$objective)))
})

test_that("pair-subset measurement model can be used in terradish fits", {
  dat <- melip_fixture(1:5)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  pairs <- rbind(c(1, 2), c(1, 3), c(2, 4), c(3, 5))
  pair_model <- pair_subset_measurement_model(leastsquares, pairs)

  fit <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              measurement_model = pair_model,
              optimizer = "bfgs",
              control = NewtonRaphsonControl(maxit = 1, verbose = FALSE),
              leverage = FALSE,
              solver = "direct")
  )

  expect_s3_class(fit, "terradish")
  expect_equal(dim(fitted(fit)), dim(melip.Fst))
})
