test_that("S3 methods for terradish objects return consistent outputs", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  fit <- fx$fit

  expect_gt(length(capture.output(print(fit))), 0)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.terradish")
  expect_s3_class(sm, "summary.radish")
  expect_gt(length(capture.output(print(sm))), 0)
  terradish_only <- fit
  class(terradish_only) <- "terradish"
  sm_terradish_only <- summary(terradish_only)
  expect_s3_class(sm_terradish_only, "summary.terradish")
  expect_match(capture.output(print(sm_terradish_only))[1],
               "Conductance surface with")
  expect_false(any(grepl("Length\\s+Class\\s+Mode",
                         capture.output(print(sm_terradish_only)))))
  expect_true(is.matrix(sm$phi_table))
  expect_equal(colnames(sm$phi_table),
               c("Estimate", "Std. Error", "Lower 95%", "Upper 95%"))
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
  collect_plot_data <- function(obj) {
    if (inherits(obj, "ggplot"))
      return(obj$data)
    do.call(rbind, lapply(obj, function(p) p$data))
  }
  first_plot <- function(obj) if (inherits(obj, "ggplot")) obj else obj[[1]]

  fit_tab <- plot(fx$fit, type = "fit")
  expect_s3_class(fit_tab, "ggplot")
  expect_true(all(c("observed", "fitted") %in% names(fit_tab$data)))
  expect_gt(nrow(fit_tab$data), 0L)
  expect_s3_class(fit_tab$theme$panel.border, "element_rect")
  expect_true(is.na(fit_tab$theme$panel.border$fill))
  expect_equal(fit_tab$theme$panel.border$colour, "black")
  expect_equal(fit_tab$theme$axis.title$face, "bold")
  expect_equal(fit_tab$theme$axis.title$size - fit_tab$theme$axis.text$size, 2)

  cond <- plot(fx$fit, type = "surface", data = fx$surface)
  expect_true(inherits(cond, "ggplot") || is.list(cond))
  cond_data <- collect_plot_data(cond)
  expect_true(all(c("x", "y", "conductance") %in% names(cond_data)))

  marg <- plot(fx$fit, type = "marginal", data = fx$surface, n = 20)
  expect_true(inherits(marg, "ggplot") || is.list(marg))
  marg_data <- collect_plot_data(marg)
  expect_true(all(c("covariate", "x", "est", "lower", "upper") %in%
                    names(marg_data)))
  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    unique(as.character(marg_data$covariate))))
  expect_true(all(vapply(split(marg_data$x, marg_data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
  expect_true(all(marg_data$lower <= marg_data$est))
  expect_true(all(marg_data$est <= marg_data$upper))
  expect_true(all(marg_data$lower <= marg_data$upper))
  expect_true(all((marg_data$upper - marg_data$lower) > 0))

  marg_response <- plot(fx$fit, type = "marginal_response",
                        data = fx$surface, n = 20)
  expect_true(inherits(marg_response, "ggplot") || is.list(marg_response))
  marg_response_data <- collect_plot_data(marg_response)
  marg_response_first <- first_plot(marg_response)
  expect_true(all(c("covariate", "x", "est", "lower", "upper") %in%
                    names(marg_response_data)))
  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    unique(as.character(marg_response_data$covariate))))
  expect_true(all(vapply(split(marg_response_data$x,
                               marg_response_data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
  expect_true(all(is.finite(marg_response_data$est)))
  expect_true(all(marg_response_data$lower <= marg_response_data$est))
  expect_true(all(marg_response_data$est <= marg_response_data$upper))
  expect_true(all((marg_response_data$upper - marg_response_data$lower) > 0))
  expect_gt(min(marg_response_data$upper - marg_response_data$lower), 1e-4)
  expect_equal(marg_response_first$scales$get_scales("y")$name,
               "Predicted genetic distance")
})

test_that("marginal plots support covariate panel selection", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  marg <- plot(fx$fit, type = "marginal", data = fx$surface, n = 8,
               marginal_covariates = "altitude")
  expect_s3_class(marg, "ggplot")
  expect_equal(unique(as.character(marg$data$covariate)),
               "altitude (original scale)")

  marg_response <- plot(fx$fit, type = "marginal_response", data = fx$surface,
                        n = 8,
                        marginal_covariates = "forestcover (original scale)")
  expect_s3_class(marg_response, "ggplot")
  expect_equal(unique(as.character(marg_response$data$covariate)),
               "forestcover (original scale)")

  expect_error(
    plot(fx$fit, type = "marginal", data = fx$surface, n = 8,
         marginal_covariates = "missing_covariate"),
    "Unknown `marginal_covariates`"
  )
})

test_that("marginal plots can infer original covariate units from surface metadata", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  surface_scaled <- fx$surface
  attr(surface_scaled$stack, "terradish_scale") <- NULL

  marg_scaled <- plot(fx$fit, type = "marginal", data = surface_scaled, n = 20)
  marg_orig <- plot(fx$fit, type = "marginal", data = fx$surface, n = 20)
  marg_scaled_data <- if (inherits(marg_scaled, "ggplot")) marg_scaled$data else do.call(rbind, lapply(marg_scaled, function(p) p$data))
  marg_orig_data <- if (inherits(marg_orig, "ggplot")) marg_orig$data else do.call(rbind, lapply(marg_orig, function(p) p$data))

  expect_true(all(c("altitude (original scale)",
                    "forestcover (original scale)") %in%
                    unique(as.character(marg_orig_data$covariate))))
  expect_false(isTRUE(all.equal(range(marg_scaled_data$x),
                                range(marg_orig_data$x))))

  x_ranges <- tapply(marg_orig_data$x, marg_orig_data$covariate, range)
  expect_equal(unname(x_ranges[["altitude (original scale)"]]),
               c(-0.0001180684, 1.0632071528),
               tolerance = 1e-6)
  expect_equal(unname(x_ranges[["forestcover (original scale)"]]),
               c(0.004033356, 0.94775),
               tolerance = 1e-6)
  expect_true(all(vapply(split(marg_orig_data$x, marg_orig_data$covariate),
                         function(z) !any(diff(z) < 0),
                         logical(1))))
})

test_that("marginal plots can clamp evaluation support to focal-site ranges", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface

  focal <- unique(surface$demes)
  non_focal <- setdiff(seq_len(nrow(surface$x)), focal)
  expect_gt(length(non_focal), 0L)

  surface_extreme <- surface
  surface_extreme$x$altitude[non_focal] <- min(surface$x$altitude[focal]) - 25

  collect_plot_data <- function(obj)
    if (inherits(obj, "ggplot")) obj$data else do.call(rbind, lapply(obj, function(p) p$data))

  marg_all <- plot(fx$fit, type = "marginal", data = surface_extreme, n = 20,
                   support = "none", clamp_covariates = "altitude")
  marg_focal <- plot(fx$fit, type = "marginal", data = surface_extreme, n = 20,
                     support = "focal", clamp_covariates = "altitude")
  marg_trim <- plot(fx$fit, type = "marginal", data = surface_extreme, n = 20,
                    support = "focal", support_probs = c(0.1, 0.9),
                    clamp_covariates = "altitude")

  data_all <- collect_plot_data(marg_all)
  data_focal <- collect_plot_data(marg_focal)
  data_trim <- collect_plot_data(marg_trim)

  alt_all <- data_all[data_all$covariate == "altitude (original scale)", "x"]
  alt_focal <- data_focal[data_focal$covariate == "altitude (original scale)", "x"]
  alt_trim <- data_trim[data_trim$covariate == "altitude (original scale)", "x"]

  expect_lt(min(alt_all), min(alt_focal))
  expect_gt(max(alt_trim) - min(alt_trim),
            0)
  expect_lt(max(alt_trim) - min(alt_trim),
            max(alt_focal) - min(alt_focal))

  expect_error(
    plot(fx$fit, type = "marginal", data = surface_extreme,
         support = "focal", support_probs = c(-0.1, 1)),
    "support_probs"
  )
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
