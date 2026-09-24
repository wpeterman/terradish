test_that("eight-neighbor graphs and graph notices follow the public contract", {
  g <- robustness_surface()
  explicit <- conductance_surface(g$stack, g$vertex_coordinates[g$demes, ], directions = 8)
  expect_equal(g$adj, explicit$adj)
  coords <- g$vertex_coordinates[g$demes, ]
  expect_warning(conductance_surface(g$stack, rbind(coords, coords[1, ])), "share")
  islands <- g$stack
  islands[c(2, 9, 10)] <- NA
  expect_error(suppressWarnings(conductance_surface(islands, coords)), "disconnected component")
  expect_message(crop_to_focal_buffer(g$stack, coords, buffer = 1), "buffer")
  controls <- list(auto_direct_max_vertices = 20L, auto_amg_min_vertices = 30L,
                   auto_direct_max_rhs = 1L)
  expect_identical(.terradish_resolve_solver(g, "auto", 64L, controls)$type, "amg")
  expect_identical(.terradish_resolve_solver(g, "direct", 64L, controls)$type, "direct")
  state <- new.env(parent = emptyenv())
  expect_message(.terradish_parallel_notice(g, 2L, "direct", NULL, state), "PSOCK")
  expect_message(.terradish_parallel_notice(g, 2L, "direct", NULL, state), NA)
})

test_that("within-diagonal covariance gives an explicit diagnostic", {
  g <- robustness_surface()
  S <- simulate_covariance_response(c(x = 0.3), ~x, g, nu = 100, seed = 93)$covariance
  attr(S, "diagonal") <- "within"
  expect_warning(terradish(S ~ 1, g, measurement_model = wishart_covariance, nu = 100),
    "diagonal = 'within'.*gower")
  expect_error(suppressWarnings(terradish(S ~ 1, g,
    measurement_model = wishart_covariance)), "nu.*positive")
})

test_that("both families fit combined geographic and climate covariates", {
  g <- robustness_surface()
  g <- conductance_surface(g$stack, terra::xyFromCell(g$stack,
    c(1, 5, 8, 12, 20, 27, 33, 37, 45, 53, 57, 64)))
  coords <- g$vertex_coordinates[g$demes, ]
  z <- pairwise_covariates(pairwise_endpoint_covariates(data.frame(climate = sin(1:12))),
                          geographic = dist(coords))
  S <- simulate_covariance_response(c(x = 0.8), ~x, g, nu = 1000, seed = 301)$covariance
  kernels <- attr(wishart_covariates(z, model = "wishart_covariance"), "kernel_covariates")
  S <- S + 0.1 * kernels[, , 1] + 0.05 * kernels[, , 2]
  D <- dist_from_cov(S)
  for (family in c("mlpe", "wishart")) {
    model <- if (family == "mlpe") mlpe_covariates(z) else
      wishart_covariates(z, model = "wishart_covariance")
    response <- if (family == "mlpe") D else S
    fit <- suppressWarnings(terradish(response ~ x, g, measurement_model = model,
      nu = 100, control = NewtonRaphsonControl(maxit = 100)))
    expect_true(is.finite(fit$loglik))
    expect_equal(fit$convergence$code, 0L)
    expect_true(all(colnames(z) %in%
      sub("^(gamma_|lambda_)", "", rownames(fit$fit$phi))))
  }
})

test_that("stored scaling is applied by layer name to new rasters", {
  r <- robustness_surface()$stack
  r <- c(r, r^2)
  names(r) <- c("x", "y")
  reference <- scale_covariates(r)
  shifted <- r + 4
  applied <- scale_covariates(shifted, reference = reference)
  recipe <- attr(reference, "terradish_scale")
  expect_equal(terra::values(applied), sweep(sweep(terra::values(shifted), 2,
    vapply(recipe, `[[`, numeric(1), "center"), "-"), 2,
    vapply(recipe, `[[`, numeric(1), "scale"), "/"))
  reversed <- scale_covariates(shifted[[2:1]], reference = recipe)
  expect_equal(terra::values(reversed), terra::values(applied)[, 2:1])
  names(shifted) <- c("x", "unknown")
  expect_error(scale_covariates(shifted, reference = recipe), "reference|Reference")
})

test_that("combined pairwise inputs preserve ordering for both likelihoods", {
  sites <- data.frame(climate = c(2, 4, 1, 7, 5, 9))
  climate <- pairwise_endpoint_covariates(sites)
  geographic <- dist(cbind(1:6, c(1, 3, 2, 8, 5, 7)))
  z <- pairwise_covariates(climate, geographic = geographic)
  expect_equal(unname(z[, "geographic"]), as.numeric(geographic))
  expect_s3_class(z, "terradish_pairwise_covariates")
  expect_type(mlpe_covariates(z), "closure")
  expect_type(wishart_covariates(z), "closure")
  expect_error(mlpe_covariates(matrix(1, 5, 2)), "n\\(n-1\\)/2")
  expect_error(pairwise_covariates(bad = matrix(1, 3, 3)), "zero diagonal")
  expect_error(pairwise_covariates(climate, short = dist(1:3)), "same number")
  subset <- .subset_pairwise_endpoint_covariates(z, c(1, 3, 6))
  expect_equal(unname(subset[, "geographic"]), as.numeric(as.dist(as.matrix(geographic)[c(1, 3, 6), c(1, 3, 6)])))
})

test_that("spline monotonicity follows derivative signs over focal support", {
  x <- data.frame(x = seq(-2, 2, length.out = 51))
  f <- smooth_loglinear_conductance(~s(x, basis = "bs", df = 4), x)
  B <- .smooth_loglinear_model_matrix(~s(x, basis = "bs", df = 4), x, 4, "bs", 3, FALSE)
  linear <- qr.solve(B, x$x - mean(x$x))
  hump <- qr.solve(B, -(x$x^2 - mean(x$x^2)))
  names(linear) <- names(hump) <- names(attr(f, "default"))
  a <- .spline_monotonicity(f, linear, x)
  b <- .spline_monotonicity(f, hump, x)
  expect_true(a$monotone)
  expect_identical(a$direction, "increasing")
  expect_false(b$monotone)
  expect_equal(b$derivative_sign_changes, 1)
  expect_true(.spline_monotonicity(f, hump, x[x$x < 0, , drop = FALSE])$monotone)
})
