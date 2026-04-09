#' Plot methods for fitted terradish models
#'
#' Three plot types for objects returned by \code{\link{terradish}}:
#' \describe{
#'   \item{\code{"fit"}}{Observed vs. fitted pairwise genetic distances.}
#'   \item{\code{"surface"}}{Fitted conductance surface with asymptotic confidence
#'     interval bounds, displayed in three side-by-side panels.}
#'   \item{\code{"marginal"}}{Marginal effect of each raster covariate on
#'     conductance, varying one covariate at a time while averaging over the
#'     observed values of all other model covariates, with a 95\% pointwise
#'     confidence band via the delta method.}
#'   \item{\code{"marginal_response"}}{Approximate response-scale marginal
#'     effects obtained by mapping the marginal conductance curve through the
#'     fitted measurement-model mean and adding predictive bands that combine
#'     \code{theta} uncertainty, conditional \code{phi} uncertainty, and the
#'     residual variance implied by \code{tau}.}
#' }
#'
#' @param x A fitted \code{terradish} object.
#' @param type One of \code{"marginal_response"} (default), \code{"fit"},
#'   \code{"surface"}, or \code{"marginal"}.
#' @param data The \code{terradish_graph} used when fitting \code{x}. Required
#'   for \code{type = "surface"}, \code{type = "marginal"}, and
#'   \code{type = "marginal_response"}.
#' @param covariates Optional \code{SpatRaster} of covariates on their
#'   \strong{original (unscaled)} scale. When supplied for
#'   \code{type = "marginal"}, the x-axis of each panel is back-transformed to
#'   the original units using the corresponding per-layer offset and divisor.
#'   Layer names must match the covariate names used in the model formula.
#' @param conductance_model The conductance model factory used when fitting
#'   \code{x} — e.g. \code{\link{loglinear_conductance}} (default) or
#'   \code{\link{linear_conductance}}. Required for
#'   \code{type = "marginal"} and \code{type = "marginal_response"} because a
#'   fresh model must be built at the evaluation points.
#' @param quantile Confidence level for interval bands. Default \code{0.95}.
#' @param n Number of evaluation points for each marginal effect curve.
#'   Defaults to \code{100} for \code{"marginal_response"} (each point
#'   requires a full Laplacian solve) and \code{200} for \code{"marginal"}.
#'   Supply an explicit integer to override.
#' @param ... Additional graphical parameters forwarded to
#'   \code{\link[graphics]{plot}} (for \code{"fit"}) or ignored for the other
#'   types.
#'
#' @return Invisibly:
#' \itemize{
#'   \item \code{"fit"}: a \code{ggplot} object with observed-vs-fitted
#'     pairwise distances.
#'   \item \code{"surface"}: a \code{ggplot} object with faceted conductance
#'     estimate and interval maps.
#'   \item \code{"marginal"}: a \code{ggplot} object with faceted marginal
#'     conductance effects and pointwise confidence bands.
#'   \item \code{"marginal_response"}: a \code{ggplot} object with faceted
#'     response-scale marginal effects and predictive bands.
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
#' covariates_scaled <- c(melip.altitude, melip.forestcover)
#' names(covariates_scaled) <- c("altitude", "forestcover")
#' covariates_scaled <- scale_covariates(covariates_scaled)
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
#' plot(fit, type = "marginal", data = surface)
#'
#' # Approximate response-scale marginal effects with predictive bands
#' plot(fit, type = "marginal_response", data = surface)
#' }
#'
#' @name plot.terradish
#' @method plot terradish
#' @importFrom ggplot2 aes coord_equal element_blank facet_wrap geom_abline
#'   geom_line geom_point geom_raster geom_ribbon geom_rug ggplot labs
#'   scale_fill_gradientn scale_y_continuous theme theme_bw
#' @importFrom stats lm
#' @importFrom terra global values
#' @export
plot.terradish <- function(x,
                            type = c("marginal_response", "fit", "surface",
                                     "marginal"),
                            data = NULL,
                            covariates = NULL,
                            conductance_model = loglinear_conductance,
                            quantile = 0.95,
                            n = NULL,
                            ...)
{
  type <- match.arg(type)
  n <- as.integer(if (is.null(n))
    if (identical(type, "marginal_response")) 100L else 200L
  else
    n)
  switch(type,
    fit      = .plot_terradish_fit(x, ...),
    surface  = .plot_terradish_surface(x, data, quantile),
    marginal = .plot_terradish_marginal(x, data, covariates,
                                        conductance_model, quantile, n),
    marginal_response = .plot_terradish_marginal(x, data, covariates,
                                                conductance_model, quantile, n,
                                                response_scale = TRUE)
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
                                 pch = 19,
                                 cex  = 0.6,
                                 col  = "#00000099",
                                 ...)
{
  obs  <- fit$fit$response
  pred <- fitted(fit, "distance")

  # lower triangle only — avoid duplicate pairs and the zero diagonal
  obs_v  <- obs[lower.tri(obs)]
  pred_v <- pred[lower.tri(pred)]
  plot_data <- data.frame(observed = obs_v, fitted = pred_v)
  lm_fit <- lm(observed ~ fitted, data = plot_data)

  ggplot(plot_data, aes(x = fitted, y = observed)) +
    geom_point(colour = col, shape = pch, size = 1.5 * cex, ...) +
    geom_abline(intercept = unname(coef(lm_fit)[1]),
                slope = unname(coef(lm_fit)[2]),
                colour = "steelblue", linewidth = 0.8) +
    labs(x = xlab, y = ylab, title = main) +
    theme_bw()
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

  panel_names <- c("Fitted conductance",
                   paste0("Lower ", q_pct, "% CI"),
                   paste0("Upper ", q_pct, "% CI"))

  plot_data <- do.call(rbind, Map(function(layer_name, panel_name) {
    layer_data <- as.data.frame(cond[[layer_name]], xy = TRUE, na.rm = TRUE)
    names(layer_data)[3] <- "conductance"
    layer_data$panel <- panel_name
    layer_data
  }, c("est", lower_nm, upper_nm), panel_names))
  plot_data$panel <- factor(plot_data$panel, levels = panel_names)

  ggplot(plot_data, aes(x = x, y = y, fill = conductance)) +
    geom_raster() +
    facet_wrap(~panel, nrow = 1) +
    coord_equal(expand = FALSE) +
    scale_fill_gradientn(colours = grDevices::terrain.colors(100), limits = clim,
                         name = "Conductance") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

.marginal_back_transform <- function(plot_vars, x_data, data_stack, covariates)
{
  if (is.null(covariates))
  {
    if (is.null(data_stack))
      return(NULL)

    scale_meta <- attr(data_stack, "terradish_scale")
    if (is.null(scale_meta))
      return(NULL)

    matched <- intersect(plot_vars, names(scale_meta))
    if (length(matched) == 0)
      return(NULL)

    out <- lapply(matched, function(nm) {
      list(intercept = unname(scale_meta[[nm]][["center"]]),
           slope = unname(scale_meta[[nm]][["scale"]]),
           label = paste0(nm, " (original scale)"),
           rug = unname(scale_meta[[nm]][["center"]] +
                        scale_meta[[nm]][["scale"]] * x_data[[nm]]))
    })
    names(out) <- matched
    return(out)
  }

  covariates <- .as_spatraster(covariates)
  matched <- intersect(plot_vars, names(covariates))
  if (length(matched) == 0)
  {
    warning("No names in `covariates` match model covariate names; ",
            "x-axis will show scaled values.",
            call. = FALSE)
    return(NULL)
  }

  out <- lapply(matched, function(nm) {
    if (!is.null(data_stack) &&
        inherits(data_stack, "SpatRaster") &&
        nm %in% names(data_stack) &&
        terra::ncell(data_stack[[nm]]) == terra::ncell(covariates[[nm]]))
    {
      scaled_values <- terra::values(data_stack[[nm]], dataframe = FALSE)[, 1]
      original_values <- terra::values(covariates[[nm]], dataframe = FALSE)[, 1]
      keep <- is.finite(scaled_values) & is.finite(original_values)
      if (sum(keep) >= 2L && stats::sd(scaled_values[keep]) > 0)
      {
        map_fit <- lm(original_values[keep] ~ scaled_values[keep])
        return(list(intercept = unname(coef(map_fit)[1]),
                    slope = unname(coef(map_fit)[2]),
                    label = paste0(nm, " (original scale)"),
                    rug = unname(coef(map_fit)[1] +
                                 coef(map_fit)[2] * x_data[[nm]])))
      }
    }

    bt_mean <- terra::global(covariates[[nm]], "mean", na.rm = TRUE)[[1]]
    bt_sd <- terra::global(covariates[[nm]], "sd", na.rm = TRUE)[[1]]
    list(intercept = bt_mean,
         slope = bt_sd,
         label = paste0(nm, " (original scale)"),
         rug = x_data[[nm]] * bt_sd + bt_mean)
  })
  names(out) <- matched
  out
}

.marginal_effect_summary <- function(formula, x_data, x_seq, covariate,
                                     conductance_model_factory, theta,
                                     vcov_theta, quantile)
{
  z <- qnorm((1 + quantile) / 2)
  estimates <- lower <- upper <- numeric(length(x_seq))

  for (i in seq_along(x_seq))
  {
    eval_data <- x_data
    eval_data[[covariate]] <- x_seq[i]

    # Append the original data so the temporary model matrix is full rank even
    # when the focal covariate is constant over the prediction block.
    eval_data <- rbind(eval_data, x_data)

    cond_fn <- conductance_model_factory(formula, eval_data)
    cond_vals <- cond_fn(theta)
    keep <- seq_len(nrow(x_data))

    if (identical(conductance_model_factory, loglinear_conductance))
    {
      eta <- log(cond_vals$conductance[keep])
      eta_hat <- mean(eta)
      grad_i <- vapply(seq_along(theta),
                       function(k) mean(cond_vals$df__dtheta(k)[keep] /
                                          cond_vals$conductance[keep]),
                       numeric(1))
      se_i <- sqrt(max(drop(t(grad_i) %*% vcov_theta %*% grad_i), 0))

      estimates[i] <- exp(eta_hat)
      lower[i] <- exp(eta_hat - z * se_i)
      upper[i] <- exp(eta_hat + z * se_i)
    }
    else
    {
      estimates[i] <- mean(cond_vals$conductance[keep])
      grad_i <- vapply(seq_along(theta),
                       function(k) mean(cond_vals$df__dtheta(k)[keep]),
                       numeric(1))
      se_i <- sqrt(max(drop(t(grad_i) %*% vcov_theta %*% grad_i), 0))

      lower[i] <- max(estimates[i] - z * se_i, 0)
      upper[i] <- estimates[i] + z * se_i
    }
  }

  data.frame(est = estimates, lower = pmin(lower, upper), upper = pmax(lower, upper))
}

.marginal_response_components <- function(fit)
{
  phi <- c(fit$fit$phi[, 1])
  if (length(phi) == 0 || !all(c("alpha", "beta", "tau") %in% names(phi)))
    stop('plot(type = "marginal_response") currently requires a ',
         'regression-style measurement model with nuisance parameters ',
         '`alpha`, `beta`, and `tau`.',
         call. = FALSE)

  vcov_phi <- .safe_hessian_inverse(fit$fit$phi_hessian)
  grad_template <- numeric(length(phi))
  names(grad_template) <- names(phi)
  grad_template["alpha"] <- 1

  extra_coef <- setdiff(names(phi), c("alpha", "beta", "tau", "rho"))
  offset <- unname(phi["alpha"])

  if (length(extra_coef) > 0)
  {
    pairwise_covariates <- attr(fit$submodels$g, "pairwise_covariates")
    if (is.null(pairwise_covariates))
      stop('plot(type = "marginal_response") cannot infer endpoint-covariate ',
           'offsets for this measurement model.',
           call. = FALSE)

    z_means <- colMeans(.as_pairwise_covariate_matrix(pairwise_covariates))
    if (!all(extra_coef %in% names(z_means)))
      stop('plot(type = "marginal_response") found nuisance coefficients ',
           'without matching endpoint-covariate columns.',
           call. = FALSE)

    grad_template[extra_coef] <- z_means[extra_coef]
    offset <- offset + sum(phi[extra_coef] * z_means[extra_coef])
  }

  tau_var <- if ("tau" %in% rownames(vcov_phi))
    max(vcov_phi["tau", "tau"], 0)
  else
    0

  list(offset = offset,
       beta = unname(phi["beta"]),
       grad_template = grad_template,
       vcov_phi = vcov_phi,
       residual_var = exp(-unname(phi["tau"]) + 0.5 * tau_var))
}

.marginal_response_summary <- function(formula, x_data, x_seq, covariate,
                                       conductance_model_factory, theta,
                                       vcov_theta, response_components,
                                       quantile, graph_data)
{
  z       <- qnorm((1 + quantile) / 2)
  N       <- nrow(x_data)
  n_focal <- length(graph_data$demes)
  estimates <- lower <- upper <- numeric(length(x_seq))

  # Analytical weight matrix for d(mean_R_ij)/d(E_kl).
  # mean_R = (2 / (n*(n-1))) * sum_{i<j} (E_ii + E_jj - 2*E_ij), so:
  #   d(mean_R)/d(E_ii)  =  2/n        (each diagonal enters n-1 pairs)
  #   d(mean_R)/d(E_ij)  = -4/(n*(n-1))  (each off-diagonal enters one pair)
  W      <- matrix(-4 / (n_focal * (n_focal - 1L)), n_focal, n_focal)
  diag(W) <- 2 / n_focal

  # Pre-compute the fixed RHS matrix (does not change with conductance).
  Zn <- .graph_rhs(graph_data, N)

  # Reuse the Cholesky symbolic factorization across evaluation points:
  # only the numeric values change between iterations, not the sparsity pattern.
  solver_reuse_state <- NULL

  for (i in seq_along(x_seq))
  {
    # Set the focal covariate to x_seq[i] across all raster cells; keep all
    # other covariates at their observed values.  Append the original rows so
    # that model.matrix stays full-rank for I()-wrapped or interacted terms.
    eval_data        <- x_data
    eval_data[[covariate]] <- x_seq[i]
    eval_data        <- rbind(eval_data, x_data)
    keep             <- seq_len(N)

    cond_fn   <- conductance_model_factory(formula, eval_data)
    cond_vals <- cond_fn(theta)

    conductance    <- cond_vals$conductance[keep]
    df__dtheta_mat <- .conductance_df_dtheta_matrix(cond_vals, theta)[keep, , drop = FALSE]

    # Solve the reduced Laplacian for the current conductance surface,
    # reusing the symbolic Cholesky factorization from the previous iteration.
    solver_state <- .terradish_solver_setup(
      graph_data, conductance,
      solver             = "direct",
      solver_reuse_state = solver_reuse_state
    )
    solve_result <- .terradish_solver_solve(solver_state, Zn)
    G  <- as.matrix(solve_result$solution)   # (N-1) x n_focal
    tG <- t(G)                               # n_focal x (N-1)

    # Carry the Cholesky factor forward: update() reuses the symbolic
    # factorization so only the numeric refactorization is repeated.
    solver_reuse_state <- list(
      type      = "direct",
      factor    = solver_state$factor,
      signature = solver_state$signature
    )

    # Genuine effective resistance distances among focal populations.
    E_eval <- graph_rhs_crossprod(graph_data$demes, N, G)  # n_focal x n_focal
    R_eval <- dist_from_cov(as.matrix(E_eval))             # n_focal x n_focal
    mean_R <- mean(R_eval[lower.tri(R_eval)])

    # Fitted response on the observed genetic-distance scale.
    estimates[i] <- response_components$offset +
      response_components$beta * mean_R

    # Gradient of mean_R w.r.t. theta via the adjoint (backpropagation) method.
    # This mirrors exactly what terradish_algorithm does for dl/d(theta), but
    # substituting W (the gradient of mean_R w.r.t. E) in place of the
    # measurement-model gradient dl/dE — no additional Laplacian solves needed.
    W_dQnG <- W %*% tG
    dl_dC  <- backpropagate_laplacian_to_conductance(W_dQnG, tG, graph_data$adj)
    grad_mean_R_theta <- c(crossprod(df__dtheta_mat, c(dl_dC)))

    # Delta-method variance: theta uncertainty propagated through beta * mean_R,
    # plus phi (intercept/slope/precision) uncertainty, plus residual variance.
    grad_phi         <- response_components$grad_template
    grad_phi["beta"] <- grad_phi["beta"] + mean_R

    var_theta <- response_components$beta^2 *
      max(drop(t(grad_mean_R_theta) %*% vcov_theta %*% grad_mean_R_theta), 0)
    var_phi   <- max(drop(t(grad_phi) %*% response_components$vcov_phi %*%
                            grad_phi), 0)
    pred_se   <- sqrt(var_theta + var_phi +
                      max(response_components$residual_var, 0))

    lower[i] <- estimates[i] - z * pred_se
    upper[i] <- estimates[i] + z * pred_se
  }

  data.frame(est = estimates, lower = pmin(lower, upper), upper = pmax(lower, upper))
}


# Marginal effect plots
.plot_terradish_marginal <- function(fit, data, covariates,
                                      conductance_model_factory, quantile, n,
                                      response_scale = FALSE)
{
  if (is.null(data) || !inherits(data, c("terradish_graph", "radish_graph")))
    stop('Marginal plots require `data`: ',
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
  if (isTRUE(attr(conductance_model_factory, "requires_fixed_graph", exact = TRUE)))
    stop("Marginal-effect plots are not yet implemented for conductance models ",
         "that require the original raster graph during evaluation.",
         call. = FALSE)

  theta     <- coef(fit)
  formula   <- fit$formula

  # covariance matrix of theta (same calculation as summary.terradish)
  vcov_theta <- .safe_hessian_inverse(fit$fit$hessian)

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

  # Average over the observed values of non-focal covariates so marginal
  # uncertainty reflects the empirical surface, not a single mean-covariate row.
  x_data  <- data$x[, available, drop = FALSE]

  response_components <- if (isTRUE(response_scale))
    .marginal_response_components(fit)
  else
    NULL

  back_transform <- .marginal_back_transform(plot_vars, x_data,
                                             data$stack, covariates)
  curves <- vector("list", length(plot_vars))
  rugs <- vector("list", length(plot_vars))
  names(curves) <- names(rugs) <- plot_vars

  for (nm in plot_vars)
  {
    x_range <- range(x_data[[nm]])
    x_seq   <- seq(x_range[1], x_range[2], length.out = n)

    curve_i <- if (isTRUE(response_scale))
      .marginal_response_summary(
        formula = formula,
        x_data = x_data,
        x_seq = x_seq,
        covariate = nm,
        conductance_model_factory = conductance_model_factory,
        theta = theta,
        vcov_theta = vcov_theta,
        response_components = response_components,
        quantile = quantile,
        graph_data = data
      )
    else
      .marginal_effect_summary(
        formula = formula,
        x_data = x_data,
        x_seq = x_seq,
        covariate = nm,
        conductance_model_factory = conductance_model_factory,
        theta = theta,
        vcov_theta = vcov_theta,
        quantile = quantile
      )

    # x-axis: original scale if back-transform available
    if (!is.null(back_transform) && nm %in% names(back_transform))
    {
      bt       <- back_transform[[nm]]
      x_plot   <- bt$intercept + bt$slope * x_seq
      rug_vals <- bt$rug
      panel_lab <- bt$label
    } else {
      x_plot   <- x_seq
      rug_vals <- x_data[[nm]]
      panel_lab <- paste0(nm, " (scaled)")
    }

    x_order <- order(x_plot)
    x_plot <- x_plot[x_order]
    curve_i <- curve_i[x_order, , drop = FALSE]

    curves[[nm]] <- data.frame(covariate = panel_lab,
                               x = x_plot,
                               est = curve_i$est,
                               lower = curve_i$lower,
                               upper = curve_i$upper)
    rugs[[nm]] <- data.frame(covariate = panel_lab, x = rug_vals)
  }

  curve_data <- do.call(rbind, curves)
  rug_data <- do.call(rbind, rugs)
  curve_data$covariate <- factor(curve_data$covariate,
                                 levels = vapply(curves, function(z) z$covariate[1],
                                                 character(1)))
  rug_data$covariate <- factor(rug_data$covariate,
                               levels = levels(curve_data$covariate))

  ggplot(curve_data, aes(x = x, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "grey80", alpha = 0.5) +
    geom_line(linewidth = 0.6, colour = "black") +
    geom_rug(data = rug_data, aes(x = x), inherit.aes = FALSE,
             sides = "b", linewidth = 0.3) +
    facet_wrap(~covariate, scales = "free_x", ncol = min(length(plot_vars), 3L)) +
    scale_y_continuous(name = if (isTRUE(response_scale))
      "Predicted genetic distance"
    else
      "Conductance") +
    labs(x = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())
}
