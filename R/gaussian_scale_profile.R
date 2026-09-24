.fixed_gaussian_factory <- function(factory, parameter, value) {
  force(factory); force(parameter); force(value)
  fixed <- function(formula, x) {
    model <- factory(formula, x)
    scale <- .conductance_model_parameter_scale(model)[parameter]
    for (name in c("default", "lower", "upper")) {
      values <- attr(model, name, exact = TRUE)
      values[parameter] <- value / scale
      attr(model, name) <- values
    }
    model
  }
  attributes(fixed) <- attributes(factory)
  fixed
}

#' Profile a Gaussian smoothing scale
#'
#' Fixes one scale at each grid value and reoptimizes all remaining conductance
#' and measurement parameters on the fitted graph.
#' @param fit An unslimmed Gaussian scale fit with its original model closures.
#' @param layer Name of the raster layer whose scale is profiled.
#' @param n Number of grid points across the fitted scale bounds, at least five.
#'   The fitted estimate is also included if it is not already on the grid.
#' @param level Confidence level between zero and one.
#' @return A list with layer, estimate, level, interval, grid, and fit. The grid
#'   contains scale values and profile log likelihoods. The returned fit stores
#'   the interval, so \code{confint(result$fit)} uses it at the matching level.
#' @details The interval is the connected likelihood-ratio support region
#'   containing the fitted estimate, using the one-degree-of-freedom chi-square
#'   cutoff. Grid crossings are refined by root finding. Increase \code{n} to
#'   check for additional support regions. Limits stop at the original bounds.
#'   At a boundary the usual chi-square calibration is an approximation.
#'   Failed optimizations stop the calculation instead of silently dropping
#'   profile points. This can be substantially slower than Wald intervals.
#' @seealso \code{\link{gaussian_scale_summary}},
#'   \code{\link{gaussian_smoothed_loglinear_conductance}},
#'   \code{\link{terradish_inference}}
#' @examples
#' \donttest{
#' set.seed(24)
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
#'                  ymin = 0, ymax = 10, crs = "EPSG:3857")
#' terra::values(r) <- rnorm(100)
#' names(r) <- "x"
#' graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 98, 5)),
#'                              directions = 4)
#' factory <- gaussian_smoothed_loglinear_conductance(graph,
#'   sigma_lower = 0.4, sigma_upper = 1.5)
#' response <- simulate_covariance_response(c(x = 0.7, sigma.x = 0.8),
#'   ~x, graph, conductance_model = factory, tau = 1, sigma = 0.1,
#'   nu = 10000, seed = 14)$covariance
#' fit <- terradish(response ~ x, graph, conductance_model = factory,
#'   measurement_model = wishart_covariance, nu = 10000, optimizer = "newton")
#' profile <- gaussian_scale_profile(fit, "x", n = 5)
#' profile$interval
#' confint(profile$fit)
#' # The stored profile interval replaces the matching scale Wald interval.
#' }
#' @export
gaussian_scale_profile <- function(fit, layer, n = 25, level = 0.95) {
  .terradish_require_submodels(fit, "Gaussian scale profiling")
  info <- fit$gaussian_scale_info
  if (length(layer) != 1L || is.null(info) || !layer %in% info$scale_vars)
    stop("`layer` must name an estimated Gaussian scale.", call. = FALSE)
  if (length(n) != 1L || !is.finite(n) || n < 5 || n != as.integer(n))
    stop("`n` must be an integer of at least five.", call. = FALSE)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be between zero and one.", call. = FALSE)
  factory <- fit$submodels$f_factory
  graph <- get0("surface", envir = environment(factory), inherits = FALSE)
  if (is.null(graph)) stop("The original Gaussian graph is unavailable; refit before profiling.", call. = FALSE)
  S <- fit$fit$response
  if (is.null(S)) stop("The fit must retain its response matrix.", call. = FALSE)
  parameter <- paste0("sigma.", layer)
  estimate <- unname(coef(fit)[parameter])
  if (!is.finite(estimate)) stop("The scale is not identified in this fit.", call. = FALSE)
  values <- sort(unique(c(exp(seq(log(info$lower[layer]), log(info$upper[layer]),
                                 length.out = n)), estimate)))
  values <- pmin(pmax(values, info$lower[layer]), info$upper[layer])
  formula <- reformulate(attr(terms(fit$formula), "term.labels"), response = "S")
  cache <- new.env(parent = emptyenv())
  profile_loglik <- function(value) {
    if (value == estimate) return(fit$loglik)
    key <- format(value, digits = 17)
    if (!is.null(cache[[key]])) return(cache[[key]])
    initial <- coef(fit)
    initial[parameter] <- value
    candidate <- terradish(formula, data = graph,
      conductance_model = .fixed_gaussian_factory(factory, parameter, value),
      measurement_model = fit$submodels$g, nu = fit$comparison$nu,
      theta = initial, optimizer = "newton", leverage = FALSE,
      control = NewtonRaphsonControl(maxit = 200, ctol = 1e-6, ftol = 1e-8))
    if (candidate$convergence$code != 0L)
      stop("Profile optimization failed at scale ", value, ".", call. = FALSE)
    cache[[key]] <- candidate$loglik
    candidate$loglik
  }
  ll <- vapply(values, profile_loglik, numeric(1))
  if (max(ll) > fit$loglik + 1e-5 * max(1, abs(fit$loglik)))
    stop("A profile point improves the fitted optimum; refit before interpreting its interval.", call. = FALSE)
  support <- 2 * (fit$loglik - ll) - stats::qchisq(level, 1)
  center <- match(estimate, values)
  left <- center
  right <- center
  while (left > 1L && support[left - 1L] <= 0) left <- left - 1L
  while (right < length(values) && support[right + 1L] <= 0) right <- right + 1L
  cutoff <- function(value) 2 * (fit$loglik - profile_loglik(value)) - stats::qchisq(level, 1)
  tolerance <- 1e-6 * diff(range(values))
  lower <- if (left == 1L) values[1L] else
    stats::uniroot(cutoff, values[c(left - 1L, left)], tol = tolerance)$root
  upper <- if (right == length(values)) utils::tail(values, 1L) else
    stats::uniroot(cutoff, values[c(right, right + 1L)], tol = tolerance)$root
  profile <- list(layer = layer, estimate = estimate, level = level,
    interval = c(lower = lower, upper = upper),
    grid = data.frame(sigma = values, logLik = ll))
  fit$scale_profiles[[layer]] <- profile
  profile$fit <- fit
  profile
}
