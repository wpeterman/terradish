test_that("wishart_covariance returns the full measurement-model interface", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Sigma <- 0.8 * E + 0.2 * diag(3)
  S <- rWishart(1, df = 20, Sigma = Sigma)[,,1] / 20

  start <- wishart_covariance(E, S, nu = 20)
  fit <- wishart_covariance(E, S, phi = start$phi, nu = 20)

  expect_equal(length(start$phi), 2)
  expect_true(is.finite(fit$objective))
  expect_equal(dim(fit$fitted), dim(S))
  expect_equal(dim(fit$gradient), c(2, 1))
  expect_equal(dim(fit$hessian), c(2, 2))
  expect_equal(dim(fit$gradient_E), dim(E))
  expect_equal(dim(fit$partial_E), c(length(E), 2))
  expect_equal(dim(fit$partial_S), c(2, sum(lower.tri(S, diag = TRUE))))
  expect_equal(dim(fit$jacobian_E(diag(nrow(E)))), dim(E))
  expect_equal(dim(fit$jacobian_S(diag(nrow(E)))), dim(S))
})

test_that("simulate_covariance_response is reproducible and returns Wishart inputs", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta_true <- c(altitude = 0.15, forestcover = -0.2)

  sim1 <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    nsim = 2,
    seed = 1
  )

  sim2 <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    nsim = 2,
    seed = 1
  )

  expect_equal(dim(sim1$covariance), c(6, 6, 2))
  expect_equal(sim1$covariance, sim2$covariance)
  expect_equal(dim(sim1$Sigma), c(6, 6))
  expect_equal(dim(sim1$E), c(6, 6))
  expect_true(isSymmetric(sim1$Sigma))
  expect_true(isSymmetric(sim1$E))
  expect_equal(unname(sim1$theta), matrix(unname(theta_true), nrow = 1))
  expect_equal(colnames(sim1$theta), names(theta_true))
})

test_that("terradish can optimize with wishart_covariance on covariance responses", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta_true <- c(altitude = 0.15, forestcover = -0.2)

  sim <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    seed = 1
  )
  genetic_cov <- sim$covariance

  fit <- suppressWarnings(
    terradish(genetic_cov ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = wishart_covariance,
              nu = 40,
              leverage = FALSE,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_s3_class(fit, "terradish")
  expect_s3_class(fit, "radish")
  expect_true(is.finite(fit$loglik))
  expect_true(is.matrix(fitted(fit, type = "response")))
  expect_equal(dim(fitted(fit, type = "response")), dim(genetic_cov))
})

test_that("plot methods handle covariance responses and smooth conductance", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta_true <- c(altitude = 0.15, forestcover = -0.2)

  sim <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    seed = 3
  )
  genetic_cov <- sim$covariance

  fit <- suppressWarnings(
    terradish(genetic_cov ~ s(altitude, df = 3),
              data = surface,
              conductance_model = smooth_loglinear_conductance,
              measurement_model = wishart_covariance,
              nu = 40,
              leverage = FALSE,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  fit_plot <- plot(fit, type = "fit")
  expect_s3_class(fit_plot, "ggplot")
  expect_equal(nrow(fit_plot$data), sum(lower.tri(genetic_cov, diag = TRUE)))
  expect_equal(fit_plot$labels$x, "Genetic covariance (fitted)")
  expect_equal(fit_plot$labels$y, "Genetic covariance (observed)")
  expect_equal(fit_plot$data$fitted,
               fitted(fit, type = "response")[lower.tri(genetic_cov, diag = TRUE)])

  surface_plot <- plot(fit, type = "surface", data = surface)
  expect_s3_class(surface_plot, "ggplot")

  conductance_raster <- conductance(surface, fit, quantile = 0.95)
  expect_s4_class(conductance_raster, "SpatRaster")
  expect_equal(names(conductance_raster), c("est", "lower95", "upper95"))
  conductance_values <- terra::values(conductance_raster, dataframe = TRUE)
  complete_cells <- stats::complete.cases(conductance_values)
  expect_true(any(complete_cells))
  expect_true(all(is.finite(as.matrix(conductance_values[complete_cells, ]))))
  expect_true(all(conductance_values$lower95[complete_cells] <=
                    conductance_values$est[complete_cells]))
  expect_true(all(conductance_values$est[complete_cells] <=
                    conductance_values$upper95[complete_cells]))

  marginal_plot <- plot(fit, type = "marginal", data = surface, n = 8)
  expect_s3_class(marginal_plot, "ggplot")
  expect_true("altitude (original scale)" %in%
                levels(marginal_plot$data$covariate))
  expect_true(all(is.finite(marginal_plot$data$est)))

  marginal_response_plot <- plot(fit, type = "marginal_response",
                                 data = surface, n = 4)
  expect_s3_class(marginal_response_plot, "ggplot")
  expect_equal(marginal_response_plot$scales$get_scales("y")$name,
               "Predicted genetic covariance")
  expect_true("altitude (original scale)" %in%
                levels(marginal_response_plot$data$covariate))
  expect_true(all(is.finite(marginal_response_plot$data$est)))
  expect_true(all(marginal_response_plot$data$lower <=
                    marginal_response_plot$data$est))
  expect_true(all(marginal_response_plot$data$est <=
                    marginal_response_plot$data$upper))
})

test_that("wishart_covariance supports leverage on full covariance responses", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta_true <- c(0.15, -0.2)

  sim <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    seed = 2
  )
  genetic_cov <- sim$covariance

  fit <- suppressWarnings(
    terradish(genetic_cov ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = wishart_covariance,
              nu = 40,
              leverage = TRUE,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_false(is.null(fit$leverage$S))
  expect_equal(dim(fit$leverage$S), c(nrow(genetic_cov), ncol(genetic_cov), 2))
  expect_true(isSymmetric(fit$leverage$S[, , 1]))
  expect_true(all(is.finite(fit$leverage$S)))
})
