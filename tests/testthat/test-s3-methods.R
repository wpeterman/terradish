test_that("S3 methods for terradish objects return consistent outputs", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  fit <- fx$fit

  expect_gt(length(capture.output(print(fit))), 0)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.terradish")
  expect_s3_class(sm, "summary.radish")
  expect_gt(length(capture.output(print(sm))), 0)
  expect_true(is.matrix(sm$phi_table))
  expect_equal(colnames(sm$phi_table), c("Estimate", "Std. Error"))
  expect_equal(rownames(sm$phi_table), names(sm$phi))
  expect_true(is.matrix(sm$phi_vcov))

  expect_equal(coef(fit), fit$mle$theta)
  expect_true(is.matrix(fitted(fit, type = "response")))
  expect_true(is.matrix(fitted(fit, type = "distance")))
  expect_true(is.matrix(fitted(fit, type = "covariance")))
  expect_true(is.matrix(resid(fit)))
  expect_equal(dim(simulate(fit, nsim = 2)), c(nrow(fit$fit$response), ncol(fit$fit$response), 2))
  expect_s3_class(logLik(fit), "logLik")
  expect_equal(AIC(fit), fit$aic)
})

test_that("anova compares fitted terradish models", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst

  fit1 <- suppressWarnings(
    terradish(melip.Fst ~ altitude,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )
  fit2 <- suppressWarnings(
    terradish(melip.Fst ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  tab <- anova(fit1, fit2)
  expect_true(is.matrix(tab) || is.data.frame(tab))
})

test_that("legacy radish wrapper warns and keeps compatibility classes", {
  fx <- fit_fixture(keep = 1:8,
                    control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  melip.Fst <- fx$data$melip.Fst
  surface <- fx$surface
  theta0 <- coef(fx$fit)

  legacy_fit <- function()
    radish(melip.Fst ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           theta = theta0,
           control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))

  expect_warning(
    withCallingHandlers(
      legacy_fit(),
      warning = function(w) {
        if (!grepl("deprecated", conditionMessage(w), fixed = TRUE))
          invokeRestart("muffleWarning")
      }
    ),
    "deprecated"
  )

  fit <- suppressWarnings(legacy_fit())
  expect_s3_class(fit, "terradish")
  expect_s3_class(fit, "radish")
})

test_that("plot methods return expected objects", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  fit_tab <- plot(fx$fit, type = "fit")
  expect_s3_class(fit_tab, "ggplot")
  expect_true(all(c("observed", "fitted") %in% names(fit_tab$data)))
  expect_gt(nrow(fit_tab$data), 0L)

  cond <- plot(fx$fit, type = "surface", data = fx$surface)
  expect_s3_class(cond, "ggplot")
  expect_true(all(c("x", "y", "conductance", "panel") %in% names(cond$data)))

  marg <- plot(fx$fit, type = "marginal", data = fx$surface, n = 20)
  expect_s3_class(marg, "ggplot")
  expect_true(all(c("covariate", "x", "est", "lower", "upper") %in%
                    names(marg$data)))
  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    levels(marg$data$covariate)))
  expect_true(all(vapply(split(marg$data$x, marg$data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
  expect_true(all(marg$data$lower <= marg$data$est))
  expect_true(all(marg$data$est <= marg$data$upper))
  expect_true(all(marg$data$lower <= marg$data$upper))
  expect_true(all((marg$data$upper - marg$data$lower) > 0))

  marg_response <- plot(fx$fit, type = "marginal_response",
                        data = fx$surface, n = 20)
  expect_s3_class(marg_response, "ggplot")
  expect_true(all(c("covariate", "x", "est", "lower", "upper") %in%
                    names(marg_response$data)))
  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    levels(marg_response$data$covariate)))
  expect_true(all(vapply(split(marg_response$data$x,
                               marg_response$data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
  expect_true(all(is.finite(marg_response$data$est)))
  expect_true(all(marg_response$data$lower <= marg_response$data$est))
  expect_true(all(marg_response$data$est <= marg_response$data$upper))
  expect_true(all((marg_response$data$upper - marg_response$data$lower) > 0))
  expect_gt(min(marg_response$data$upper - marg_response$data$lower), 1e-4)
  expect_equal(marg_response$scales$get_scales("y")$name,
               "Predicted genetic distance")
})

test_that("marginal plots can infer original covariate units from surface metadata", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  surface_scaled <- fx$surface
  attr(surface_scaled$stack, "terradish_scale") <- NULL

  marg_scaled <- plot(fx$fit, type = "marginal", data = surface_scaled, n = 20)
  marg_orig <- plot(fx$fit, type = "marginal", data = fx$surface, n = 20)

  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    levels(marg_orig$data$covariate)))
  expect_false(isTRUE(all.equal(range(marg_scaled$data$x),
                                range(marg_orig$data$x))))

  x_ranges <- tapply(marg_orig$data$x, marg_orig$data$covariate, range)
  expect_equal(unname(x_ranges[["altitude (original scale)"]]),
               c(-0.0001180684, 1.0632071528),
               tolerance = 1e-6)
  expect_equal(unname(x_ranges[["forestcover (original scale)"]]),
               c(0.004033356, 0.94775),
               tolerance = 1e-6)
  expect_true(all(vapply(split(marg_orig$data$x, marg_orig$data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
})

test_that("legacy-only radish classes dispatch S3 methods", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  fit <- fx$fit
  class(fit) <- "radish"

  surface <- fx$surface
  class(surface) <- "radish_graph"

  expect_gt(length(capture.output(print(fit))), 0)
  sm <- summary(fit)
  expect_s3_class(sm, "summary.radish")
  class(sm) <- "summary.radish"
  expect_gt(length(capture.output(print(sm))), 0)

  expect_equal(coef(fit), fx$fit$mle$theta)
  expect_true(is.matrix(fitted(fit, type = "response")))
  expect_true(is.matrix(residuals(fit)))
  expect_equal(dim(simulate(fit, nsim = 2)),
               c(nrow(fit$fit$response), ncol(fit$fit$response), 2))
  expect_s3_class(logLik(fit), "logLik")
  expect_equal(AIC(fit), fit$aic)
  expect_true(inherits(conductance(surface, fit), "SpatRaster"))

  expect_s3_class(plot(fit, type = "fit"), "ggplot")
})
