#' Plot methods for fitted terradish models
#'
#' Three plot types for objects returned by \code{\link{terradish}}:
#' \describe{
#'   \item{\code{"fit"}}{Observed vs. fitted pairwise genetic distances.}
#'   \item{\code{"surface"}}{Fitted conductance surface with asymptotic confidence
#'     interval bounds, displayed in three side-by-side panels.}
#'   \item{\code{"marginal"}}{Marginal effect of each raster covariate on
#'     conductance, varying one covariate at a time while holding all others at
#'     their mean, with a 95\% confidence band via the delta method.}
#' }
#'
#' @param x A fitted \code{terradish} object.
#' @param type One of \code{"fit"}, \code{"surface"}, or \code{"marginal"}.
#' @param data The \code{terradish_graph} used when fitting \code{x}. Required
#'   for \code{type = "surface"} and \code{type = "marginal"}.
#' @param covariates Optional \code{SpatRaster} of covariates on their
#'   \strong{original (unscaled)} scale. When supplied for
#'   \code{type = "marginal"}, the x-axis of each panel is back-transformed to
#'   the original units using the mean and standard deviation of each layer.
#'   Layer names must match the covariate names used in the model formula.
#' @param conductance_model The conductance model factory used when fitting
#'   \code{x} — e.g. \code{\link{loglinear_conductance}} (default) or
#'   \code{\link{linear_conductance}}. Required for \code{type = "marginal"}
#'   because a fresh model must be built at the evaluation points.
#' @param quantile Confidence level for interval bands. Default \code{0.95}.
#' @param n Number of evaluation points for each marginal effect curve.
#'   Default \code{200}.
#' @param ... Additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot}} (for \code{"fit"}) or ignored for the other
#'   types.
#'
#' @return Invisibly:
#' \itemize{
#'   \item \code{"fit"}: a \code{data.frame} with columns \code{observed} and
#'     \code{fitted} (lower-triangle values only).
#'   \item \code{"surface"}: the \code{SpatRaster} returned by
#'     \code{\link{conductance}}.
#'   \item \code{"marginal"}: a named list of \code{data.frame}s, one per
#'     continuous covariate, with columns \code{x}, \code{est}, \code{lower},
#'     and \code{upper}.
#' }
#'
#' @seealso \code{\link{terradish}}, \code{\link{conductance}},
#'   \code{\link{fitted.radish}}, \code{\link{loglinear_conductance}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data(melip)
#' melip.altitude    <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords      <- terra::unwrap(melip.coords)
#'
#' # Build surface with scaled covariates
#' covariates_scaled <- c(terra::scale(melip.altitude),
#'                        terra::scale(melip.forestcover))
#' names(covariates_scaled) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates_scaled, melip.coords, directions = 8)
#'
#' fit <- terradish(melip.Fst ~ altitude + forestcover, surface,
#'                  loglinear_conductance, mlpe)
#'
#' # Observed vs. fitted pairwise genetic distances
#' plot(fit, type = "fit")
#'
#' # Fitted conductance surface with 95% CI (three-panel figure)
#' plot(fit, type = "surface", data = surface)
#'
#' # Marginal effect plots on the original covariate scale
#' covariates_orig <- c(melip.altitude, melip.forestcover)
#' names(covariates_orig) <- c("altitude", "forestcover")
#' plot(fit, type = "marginal", data = surface, covariates = covariates_orig)
#' }
#'
#' @name plot.terradish
#' @method plot terradish
#' @importFrom graphics par plot.default abline lines rug
#' @importFrom stats lm
#' @importFrom terra global
#' @export
plot.terradish <- function(x,
                            type = c("fit", "surface", "marginal"),
                            data = NULL,
                            covariates = NULL,
                            conductance_model = loglinear_conductance,
                            quantile = 0.95,
                            n = 200L,
                            ...)
{
  type <- match.arg(type)
  switch(type,
    fit      = .plot_terradish_fit(x, ...),
    surface  = .plot_terradish_surface(x, data, quantile),
    marginal = .plot_terradish_marginal(x, data, covariates,
                                        conductance_model, quantile,
                                        as.integer(n))
  )
}

#' @rdname plot.terradish
#' @method plot radish
#' @export
plot.radish <- function(x, ...) plot.terradish(x, ...)


# ---- internal helpers -------------------------------------------------------

# Observed vs. fitted plot
.plot_terradish_fit <- function(fit,
                                 xlab = "Resistance distance (fitted)",
                                 ylab = "Genetic distance (observed)",
                                 main = "",
                                 pch  = 19,
                                 cex  = 0.6,
                                 col  = "#00000099",
                                 ...)
{
  obs  <- fit$response
  pred <- fitted(fit, "distance")

  # lower triangle only — avoid duplicate pairs and the zero diagonal
  obs_v  <- obs[lower.tri(obs)]
  pred_v <- pred[lower.tri(pred)]

  plot.default(pred_v, obs_v,
               xlab = xlab, ylab = ylab, main = main,
               pch = pch, cex = cex, col = col, ...)
  abline(lm(obs_v ~ pred_v), col = "steelblue", lwd = 2)

  invisible(data.frame(observed = obs_v, fitted = pred_v))
}


# Conductance surface + CI
.plot_terradish_surface <- function(fit, data, quantile)
{
  if (is.null(data) || !inherits(data, c("terradish_graph", "radish_graph")))
    stop('plot(type = "surface") requires `data`: ',
         'supply the terradish_graph used for fitting.',
         call. = FALSE)

  if (fit$fit$boundary || is.null(fit$mle$theta))
    stop("Cannot plot conductance surface: no conductance parameters ",
         "estimated (IBD or boundary model).",
         call. = FALSE)

  cond  <- conductance(data, fit, quantile = quantile)
  q_pct <- round(100 * quantile)

  lower_nm <- paste0("lower", q_pct)
  upper_nm <- paste0("upper", q_pct)

  # shared colour range across all three panels
  clim <- range(terra::values(cond[[lower_nm]]),
                terra::values(cond[[upper_nm]]),
                na.rm = TRUE)

  old_par <- par(mfrow = c(1, 3))
  on.exit(par(old_par), add = TRUE)

  plot(cond[["est"]],      range = clim, main = "Fitted conductance")
  plot(cond[[lower_nm]],   range = clim, main = paste0("Lower ", q_pct, "% CI"))
  plot(cond[[upper_nm]],   range = clim, main = paste0("Upper ", q_pct, "% CI"))

  invisible(cond)
}


# Marginal effect plots
.plot_terradish_marginal <- function(fit, data, covariates,
                                      conductance_model_factory, quantile, n)
{
  if (is.null(data) || !inherits(data, c("terradish_graph", "radish_graph")))
    stop('plot(type = "marginal") requires `data`: ',
         'supply the terradish_graph used for fitting.',
         call. = FALSE)

  if (fit$fit$boundary || is.null(fit$mle$theta))
    stop("Cannot plot marginal effects: no conductance parameters estimated ",
         "(IBD or boundary model).",
         call. = FALSE)

  if (!inherits(conductance_model_factory,
                c("terradish_conductance_model_factory",
                  "radish_conductance_model_factory")))
    stop("`conductance_model` must be a conductance model factory ",
         "(e.g. loglinear_conductance).",
         call. = FALSE)

  theta     <- coef(fit)
  formula   <- fit$formula

  # covariance matrix of theta (same calculation as summary.terradish)
  vcov_theta <- tryCatch(
    solve(fit$fit$hessian),
    error = function(e) MASS::ginv(fit$fit$hessian)
  )

  # raw variable names from the formula (not model-matrix column names)
  raw_vars <- all.vars(formula)

  # restrict to continuous covariates present in data$x
  available <- intersect(raw_vars, names(data$x))
  if (length(available) == 0)
    stop("No formula variables found in `data$x`.", call. = FALSE)

  is_cont <- vapply(data$x[, available, drop = FALSE], is.numeric, logical(1))
  if (any(!is_cont))
    warning("Skipping categorical covariates in marginal plots: ",
            paste(available[!is_cont], collapse = ", "),
            call. = FALSE)
  plot_vars <- available[is_cont]
  if (length(plot_vars) == 0)
    stop("No continuous covariates available to plot.", call. = FALSE)

  # mean values for "hold others constant"
  x_data  <- data$x[, available, drop = FALSE]
  x_means <- colMeans(x_data)

  # optional back-transformation to original scale
  back_transform <- NULL
  if (!is.null(covariates))
  {
    covariates <- .as_spatraster(covariates)
    matched    <- intersect(plot_vars, names(covariates))
    if (length(matched) == 0)
      warning("No names in `covariates` match model covariate names; ",
              "x-axis will show scaled values.",
              call. = FALSE)
    else
    {
      bt <- lapply(matched, function(nm) {
        list(
          mean = terra::global(covariates[[nm]], "mean", na.rm = TRUE)[[1]],
          sd   = terra::global(covariates[[nm]], "sd",   na.rm = TRUE)[[1]]
        )
      })
      names(bt)   <- matched
      back_transform <- bt
    }
  }

  # panel layout: up to 3 columns
  n_panels <- length(plot_vars)
  nc <- min(n_panels, 3L)
  nr <- ceiling(n_panels / nc)

  old_par <- par(mfrow = c(nr, nc))
  on.exit(par(old_par), add = TRUE)

  out <- vector("list", n_panels)
  names(out) <- plot_vars

  for (nm in plot_vars)
  {
    x_range <- range(x_data[[nm]])
    x_seq   <- seq(x_range[1], x_range[2], length.out = n)

    # data frame with nm varying, all other raw variables at their means
    new_data <- as.data.frame(
      matrix(x_means[available], nrow = n, ncol = length(available),
             byrow = TRUE, dimnames = list(NULL, available))
    )
    new_data[[nm]] <- x_seq

    # build fresh conductance model over the new evaluation points
    cond_fn   <- conductance_model_factory(formula, new_data)
    cond_vals <- cond_fn(theta)

    c_est <- cond_vals$conductance
    ci    <- cond_vals$confint(theta    = theta,
                               vcov     = vcov_theta,
                               quantile = quantile,
                               scale    = "conductance")
    c_lower <- ci[, "lower"]
    c_upper <- ci[, "upper"]

    # x-axis: original scale if back-transform available
    if (!is.null(back_transform) && nm %in% names(back_transform))
    {
      bt       <- back_transform[[nm]]
      x_plot   <- x_seq   * bt$sd + bt$mean
      rug_vals <- x_data[[nm]] * bt$sd + bt$mean
      xlab_k   <- paste0(nm, " (original scale)")
    } else {
      x_plot   <- x_seq
      rug_vals <- x_data[[nm]]
      xlab_k   <- paste0(nm, " (scaled)")
    }

    ylim <- range(c_lower, c_upper, finite = TRUE)

    plot.default(x_plot, c_est,
                 type = "l", lwd = 2,
                 xlab = xlab_k, ylab = "Conductance",
                 main = nm, ylim = ylim)
    lines(x_plot, c_lower, lty = 2, col = "grey40")
    lines(x_plot, c_upper, lty = 2, col = "grey40")
    rug(rug_vals, ticksize = 0.03)

    out[[nm]] <- data.frame(x     = x_plot,
                            est   = c_est,
                            lower = c_lower,
                            upper = c_upper)
  }

  invisible(out)
}
