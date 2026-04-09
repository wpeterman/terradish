.gaussian_scale_strip_identity <- function(expr)
{
  if (is.call(expr))
  {
    if (length(expr) == 2L && identical(expr[[1]], as.name("I")))
      return(.gaussian_scale_strip_identity(expr[[2]]))

    for (i in 2:length(expr))
      expr[[i]] <- .gaussian_scale_strip_identity(expr[[i]])
  }

  expr
}

.gaussian_scale_formula_expression <- function(label)
{
  expr_text <- gsub(":", "*", label, fixed = TRUE)
  expr <- tryCatch(str2lang(expr_text),
                   error = function(e)
                     stop("Could not parse Gaussian scale-aware term `", label, "`",
                          call. = FALSE))
  .gaussian_scale_strip_identity(expr)
}

.gaussian_scale_term_spec <- function(label, raw_vars)
{
  expr <- .gaussian_scale_formula_expression(label)

  first <- stats::setNames(vector("list", length(raw_vars)), raw_vars)
  second <- stats::setNames(vector("list", length(raw_vars)), raw_vars)

  for (v in raw_vars)
  {
    first[[v]] <- tryCatch(
      .gaussian_scale_strip_identity(D(expr, v)),
      error = function(e)
        stop("Gaussian scale-aware conductance does not support term `", label,
             "` because its derivative with respect to `", v,
             "` could not be constructed analytically.",
             call. = FALSE)
    )
    second[[v]] <- stats::setNames(vector("list", length(raw_vars)), raw_vars)
  }

  for (v in raw_vars)
    for (w in raw_vars)
      second[[v]][[w]] <- tryCatch(
        .gaussian_scale_strip_identity(D(first[[v]], w)),
        error = function(e)
          stop("Gaussian scale-aware conductance does not support term `", label,
               "` because its second derivative with respect to `", v,
               "` and `", w, "` could not be constructed analytically.",
               call. = FALSE)
      )

  list(label = label, expr = expr, first = first, second = second)
}

.gaussian_scale_formula_spec <- function(formula, stack_names)
{
  stopifnot(inherits(formula, "formula"))

  tt <- terms(formula)
  raw_vars <- unique(all.vars(delete.response(tt)))
  term_labels <- attr(tt, "term.labels")

  if (!length(raw_vars))
    stop("Gaussian scale-aware conductance requires at least one raster covariate",
         call. = FALSE)
  if (!length(term_labels))
    stop("Gaussian scale-aware conductance requires at least one non-intercept raster term",
         call. = FALSE)
  if (!all(raw_vars %in% stack_names))
    stop("All formula variables must be present in the retained raster stack",
         call. = FALSE)

  list(
    raw_vars = raw_vars,
    term_labels = term_labels,
    term_specs = lapply(term_labels, .gaussian_scale_term_spec, raw_vars = raw_vars)
  )
}

.validate_gaussian_scale_vars <- function(scale_vars, vars)
{
  if (is.null(scale_vars))
    return(vars)

  scale_vars <- as.character(scale_vars)
  if (!all(scale_vars %in% vars))
    stop("`scale_vars` must be a subset of the raster variables used in the formula",
         call. = FALSE)

  unique(scale_vars)
}

.coerce_gaussian_scale_bound <- function(x, scale_vars, default, arg_name)
{
  if (!length(scale_vars))
    return(stats::setNames(numeric(), character()))

  if (is.null(x))
    return(stats::setNames(rep(default, length(scale_vars)), scale_vars))

  x <- c(x)
  if (is.null(names(x)))
  {
    if (length(x) == 1L)
      x <- rep(x, length(scale_vars))
    if (length(x) != length(scale_vars))
      stop("`", arg_name, "` must have length 1 or one entry per scaled raster",
           call. = FALSE)
    names(x) <- scale_vars
  }
  else
  {
    if (!all(scale_vars %in% names(x)))
      stop("Named `", arg_name, "` values must include all `scale_vars`",
           call. = FALSE)
    x <- x[scale_vars]
  }

  if (any(!is.finite(x)) || any(x <= 0))
    stop("`", arg_name, "` must contain finite positive values", call. = FALSE)

  x
}

.gaussian_scale_extent_diagonal <- function(stack)
{
  ext <- terra::ext(stack)
  width <- ext$xmax - ext$xmin
  height <- ext$ymax - ext$ymin
  sqrt(width^2 + height^2)
}

.gaussian_scale_sigma_bounds <- function(surface, scale_vars,
                                         sigma_lower = NULL,
                                         sigma_upper = NULL)
{
  if (!length(scale_vars))
    return(list(lower = stats::setNames(numeric(), character()),
                upper = stats::setNames(numeric(), character())))

  stack <- .as_spatraster(surface$stack)
  mean_res <- mean(abs(terra::res(stack[[1]])))
  extent_diag <- .gaussian_scale_extent_diagonal(stack)

  lower_default <- mean_res / 2
  upper_default <- max(extent_diag, lower_default * 4)

  lower <- .coerce_gaussian_scale_bound(sigma_lower, scale_vars,
                                        lower_default, "sigma_lower")
  upper <- .coerce_gaussian_scale_bound(sigma_upper, scale_vars,
                                        upper_default, "sigma_upper")

  if (any(lower >= upper))
    stop("Each `sigma_lower` bound must be strictly less than `sigma_upper`",
         call. = FALSE)

  list(lower = lower, upper = upper)
}

.fft_pad_matrix <- function(x, nrow_pad, ncol_pad)
{
  out <- matrix(0, nrow_pad, ncol_pad)
  out[seq_len(nrow(x)), seq_len(ncol(x))] <- x
  out
}

.fft_convolution_crop <- function(fft_x, fft_k, prep)
{
  conv_full <- Re(stats::fft(fft_x * fft_k, inverse = TRUE)) /
    (prep$nrow_pad * prep$ncol_pad)
  conv_full[prep$i1:prep$i2, prep$j1:prep$j2, drop = FALSE]
}

.gaussian_scale_prepare_layer <- function(layer, rowcol)
{
  mat <- terra::as.matrix(layer, wide = TRUE)
  mask <- is.finite(mat)
  values0 <- mat
  values0[!mask] <- 0

  nr <- nrow(mat)
  nc <- ncol(mat)
  nrow_pad <- nr + nr - 1L
  ncol_pad <- nc + nc - 1L
  center_row <- ceiling(nr / 2)
  center_col <- ceiling(nc / 2)
  y_res <- abs(terra::res(layer)[2])
  x_res <- abs(terra::res(layer)[1])

  row_offset <- (seq_len(nr) - center_row) * y_res
  col_offset <- (seq_len(nc) - center_col) * x_res
  dist2 <- outer(row_offset^2, col_offset^2, `+`)

  list(
    nrow = nr,
    ncol = nc,
    nrow_pad = nrow_pad,
    ncol_pad = ncol_pad,
    i1 = floor(nr / 2) + 1L,
    i2 = floor(nr / 2) + nr,
    j1 = floor(nc / 2) + 1L,
    j2 = floor(nc / 2) + nc,
    dist2 = dist2,
    rowcol = rowcol,
    na_mask = !mask,
    value_fft = stats::fft(.fft_pad_matrix(values0, nrow_pad, ncol_pad)),
    mask_fft = stats::fft(.fft_pad_matrix(mask * 1, nrow_pad, ncol_pad)),
    raw_active = mat[cbind(rowcol[, 1], rowcol[, 2])]
  )
}

.gaussian_scale_kernel_fft <- function(prep, sigma)
{
  sigma2 <- sigma^2
  ratio <- prep$dist2 / sigma2
  kernel <- exp(-prep$dist2 / (2 * sigma2))
  dkernel <- kernel * ratio / sigma
  d2kernel <- kernel * (ratio^2 - 3 * ratio) / sigma2

  list(
    kernel_fft = stats::fft(.fft_pad_matrix(kernel, prep$nrow_pad, prep$ncol_pad)),
    dkernel_fft = stats::fft(.fft_pad_matrix(dkernel, prep$nrow_pad, prep$ncol_pad)),
    d2kernel_fft = stats::fft(.fft_pad_matrix(d2kernel, prep$nrow_pad, prep$ncol_pad))
  )
}

.gaussian_scale_standardize <- function(x, dx, d2x)
{
  keep <- is.finite(x)
  x_keep <- x[keep]
  dx_keep <- dx[keep]
  d2x_keep <- d2x[keep]

  mu <- mean(x_keep)
  dmu <- mean(dx_keep)
  d2mu <- mean(d2x_keep)

  centered <- x_keep - mu
  dcentered <- dx_keep - dmu
  d2centered <- d2x_keep - d2mu

  if (length(x_keep) <= 1L)
  {
    scale <- 1
    dscale <- 0
    d2scale <- 0
  }
  else
  {
    variance <- sum(centered^2) / (length(x_keep) - 1L)
    if (!is.finite(variance) || variance <= 0)
    {
      scale <- 1
      dscale <- 0
      d2scale <- 0
    }
    else
    {
      dvariance <- 2 * sum(centered * dcentered) / (length(x_keep) - 1L)
      d2variance <- 2 * sum(dcentered^2 + centered * d2centered) /
        (length(x_keep) - 1L)

      scale <- sqrt(variance)
      dscale <- dvariance / (2 * scale)
      d2scale <- d2variance / (2 * scale) - dvariance^2 / (4 * scale^3)
    }
  }

  out <- x
  dout <- dx
  d2out <- d2x

  out[keep] <- centered / scale
  dout[keep] <- dcentered / scale - centered * dscale / scale^2
  d2out[keep] <- d2centered / scale -
    2 * dcentered * dscale / scale^2 -
    centered * d2scale / scale^2 +
    2 * centered * dscale^2 / scale^3

  list(value = out, deriv = dout, second = d2out)
}

.gaussian_scale_layer_values <- function(prep, sigma, standardize = TRUE)
{
  if (!is.finite(sigma) || sigma <= 0)
    stop("All sigma parameters must be finite and strictly positive", call. = FALSE)

  kernel <- .gaussian_scale_kernel_fft(prep, sigma)

  num <- .fft_convolution_crop(prep$value_fft, kernel$kernel_fft, prep)
  den <- .fft_convolution_crop(prep$mask_fft, kernel$kernel_fft, prep)
  dnum <- .fft_convolution_crop(prep$value_fft, kernel$dkernel_fft, prep)
  dden <- .fft_convolution_crop(prep$mask_fft, kernel$dkernel_fft, prep)
  d2num <- .fft_convolution_crop(prep$value_fft, kernel$d2kernel_fft, prep)
  d2den <- .fft_convolution_crop(prep$mask_fft, kernel$d2kernel_fft, prep)

  smooth <- num / den
  dsmooth <- (dnum * den - num * dden) / den^2
  d2smooth <- d2num / den -
    num * d2den / den^2 -
    2 * dnum * dden / den^2 +
    2 * num * dden^2 / den^3

  smooth[den <= 0] <- NA_real_
  dsmooth[den <= 0] <- NA_real_
  d2smooth[den <= 0] <- NA_real_
  smooth[prep$na_mask] <- NA_real_
  dsmooth[prep$na_mask] <- NA_real_
  d2smooth[prep$na_mask] <- NA_real_

  idx <- cbind(prep$rowcol[, 1], prep$rowcol[, 2])
  value <- smooth[idx]
  deriv <- dsmooth[idx]
  second <- d2smooth[idx]

  if (isTRUE(standardize))
    return(.gaussian_scale_standardize(value, deriv, second))

  list(value = value, deriv = deriv, second = second)
}

.gaussian_scale_conductance_default <- function(beta_names, scale_vars,
                                                sigma_bounds, surface)
{
  sigma_names <- paste0("sigma.", scale_vars)
  out <- rep(0, length(beta_names) + length(scale_vars))
  names(out) <- c(beta_names, sigma_names)

  if (length(scale_vars))
  {
    start_sigma <- pmin(
      pmax(mean(abs(terra::res(surface$stack[[1]]))), sigma_bounds$lower),
      sigma_bounds$upper
    )
    out[sigma_names] <- start_sigma
  }

  out
}

.gaussian_scale_eval_environment <- function(values)
{
  env <- list2env(values, parent = baseenv())
  env$I <- function(x) x
  env
}

.gaussian_scale_eval_expression <- function(expr, env)
  eval(expr, envir = env)

.gaussian_scale_probability_label <- function(probability)
  gsub("\\.", "_", sprintf("%g", 100 * probability))

.gaussian_scale_unit_suffix <- function(unit)
  gsub("[^[:alnum:]_]+", "_", tolower(unit))

#' Gaussian scale-aware log-linear conductance model
#'
#' Creates a conductance-model factory that smooths selected raster layers with a
#' native Gaussian FFT operator and treats the corresponding Gaussian scale
#' parameters as part of the optimized conductance parameter vector.
#'
#' @param surface A \code{\link{conductance_surface}} object created with
#'   \code{saveStack = TRUE}. The retained raster stack is used as the native
#'   smoothing domain during optimization.
#' @param scale_vars Optional character vector naming the raster layers whose
#'   Gaussian scales of effect should be estimated. Defaults to all numeric
#'   raster variables present in the model formula.
#' @param standardize Should each smoothed layer be centered and scaled across
#'   the active graph cells at every parameter evaluation? Default \code{TRUE}.
#' @param sigma_lower,sigma_upper Optional lower and upper bounds for the
#'   Gaussian scale parameters on their natural distance scale. Scalar values
#'   are recycled across all scaled rasters; named vectors may be used to set
#'   different bounds by raster name. When omitted, \code{sigma_lower} defaults
#'   to half of the mean raster cell width and \code{sigma_upper} defaults to
#'   the diagonal length of the retained raster extent.
#'
#' @details
#' This is a scale-aware extension of \code{\link{loglinear_conductance}}.
#' Given a model formula such as
#' \code{~ altitude * forestcover + I(altitude^2)}, the fitted parameter vector
#' contains both conductance coefficients and one natural-scale \code{sigma}
#' parameter for each scaled raster layer. The spatial kernel is Gaussian and is
#' evaluated with an FFT-based normalized convolution over the retained raster
#' stack from \code{surface}. Using a fixed full-grid kernel keeps the sigma
#' derivatives smooth during optimization.
#'
#' The current implementation supports numeric raster formulas built from raw
#' raster names, interactions, and polynomial/arithmetic terms that can be
#' differentiated analytically by \code{\link[stats:D]{D}} after smoothing the base
#' raster layers. Factor-valued raster layers are still rejected.
#'
#' For this conductance model, \code{\link{terradish}} prefers
#' \code{optimizer = "bfgs"} when \code{optimizer = "auto"}, and leverage /
#' \code{partial_X} calculations are currently disabled.
#'
#' @return A function of class \code{terradish_conductance_model_factory} that
#'   can be supplied as the \code{conductance_model=} argument to
#'   \code{\link{terradish}}.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' covariates <- c(melip.altitude, melip.forestcover)
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, melip.coords, directions = 8,
#'                                saveStack = TRUE)
#'
#' gaussian_model <- gaussian_smoothed_loglinear_conductance(
#'   surface,
#'   sigma_lower = 100,
#'   sigma_upper = 5000
#' )
#' fit <- terradish(melip.Fst ~ altitude * forestcover + I(altitude^2),
#'                  data = surface,
#'                  conductance_model = gaussian_model,
#'                  measurement_model = mlpe,
#'                  optimizer = "bfgs")
#' }
#'
#' @export
gaussian_smoothed_loglinear_conductance <- function(surface,
                                                    scale_vars = NULL,
                                                    standardize = TRUE,
                                                    sigma_lower = NULL,
                                                    sigma_upper = NULL)
{
  stopifnot(inherits(surface, c("terradish_graph", "radish_graph")))
  if (is.null(surface$stack))
    stop("`surface` must retain its raster stack; use `conductance_surface(..., saveStack = TRUE)`",
         call. = FALSE)

  stack <- .as_spatraster(surface$stack)
  active_cells <- terra::cellFromXY(stack[[1]], surface$vertex_coordinates)
  rowcol <- terra::rowColFromCell(stack[[1]], active_cells)

  factory <- function(formula, x)
  {
    spec <- .gaussian_scale_formula_spec(formula, names(stack))
    vars <- spec$raw_vars
    chosen_scale_vars <- .validate_gaussian_scale_vars(scale_vars, vars)

    if (nrow(x) != nrow(surface$x))
      stop("Gaussian scale-aware conductance currently requires the original graph design matrix; ",
           "marginal-effect prediction grids and coarse-raster approximations are not yet supported",
           call. = FALSE)

    if (any(vapply(vars, function(nm) is.factor(stack[[nm]]), logical(1))))
      stop("Gaussian scale-aware conductance currently supports numeric raster layers only",
           call. = FALSE)

    layer_preps <- lapply(vars, function(nm)
      .gaussian_scale_prepare_layer(stack[[nm]], rowcol))
    names(layer_preps) <- vars

    beta_names <- spec$term_labels
    sigma_names <- paste0("sigma.", chosen_scale_vars)
    sigma_bounds <- .gaussian_scale_sigma_bounds(
      surface = surface,
      scale_vars = chosen_scale_vars,
      sigma_lower = sigma_lower,
      sigma_upper = sigma_upper
    )
    default <- .gaussian_scale_conductance_default(
      beta_names = beta_names,
      scale_vars = chosen_scale_vars,
      sigma_bounds = sigma_bounds,
      surface = surface
    )

    lower <- rep(-Inf, length(default))
    upper <- rep(Inf, length(default))
    names(lower) <- names(upper) <- names(default)
    if (length(chosen_scale_vars))
    {
      lower[sigma_names] <- sigma_bounds$lower[chosen_scale_vars]
      upper[sigma_names] <- sigma_bounds$upper[chosen_scale_vars]
    }

    conductance_model <- function(theta)
    {
      theta <- c(theta)
      stopifnot(length(theta) == length(default))
      names(theta) <- names(default)

      beta <- theta[beta_names]
      sigma <- theta[sigma_names]

      base_value <- stats::setNames(vector("list", length(vars)), vars)
      base_deriv <- stats::setNames(vector("list", length(vars)), vars)
      base_second <- stats::setNames(vector("list", length(vars)), vars)

      for (nm in vars)
      {
        if (nm %in% chosen_scale_vars)
        {
          scaled <- .gaussian_scale_layer_values(
            prep = layer_preps[[nm]],
            sigma = sigma[[paste0("sigma.", nm)]],
            standardize = standardize
          )
          base_value[[nm]] <- scaled$value
          base_deriv[[nm]] <- scaled$deriv
          base_second[[nm]] <- scaled$second
        }
        else
        {
          values <- layer_preps[[nm]]$raw_active
          zeros <- rep(0, length(values))
          if (isTRUE(standardize))
          {
            centered <- .gaussian_scale_standardize(values, zeros, zeros)
            values <- centered$value
          }
          base_value[[nm]] <- values
          base_deriv[[nm]] <- zeros
          base_second[[nm]] <- zeros
        }
      }

      eval_env <- .gaussian_scale_eval_environment(base_value)
      n <- nrow(surface$x)
      p <- length(beta_names)
      q <- length(chosen_scale_vars)

      X <- matrix(0, nrow = n, ncol = p)
      colnames(X) <- beta_names

      dX <- stats::setNames(vector("list", q), chosen_scale_vars)
      d2X <- stats::setNames(vector("list", q), chosen_scale_vars)
      for (j in chosen_scale_vars)
      {
        dX[[j]] <- matrix(0, nrow = n, ncol = p)
        colnames(dX[[j]]) <- beta_names
        d2X[[j]] <- stats::setNames(vector("list", q), chosen_scale_vars)
        for (k in chosen_scale_vars)
        {
          d2X[[j]][[k]] <- matrix(0, nrow = n, ncol = p)
          colnames(d2X[[j]][[k]]) <- beta_names
        }
      }

      for (col in seq_along(spec$term_specs))
      {
        term <- spec$term_specs[[col]]
        X[, col] <- .gaussian_scale_eval_expression(term$expr, eval_env)

        for (j in chosen_scale_vars)
        {
          partial1 <- .gaussian_scale_eval_expression(term$first[[j]], eval_env)
          dX[[j]][, col] <- partial1 * base_deriv[[j]]
        }

        for (j in chosen_scale_vars)
        {
          partial1_j <- .gaussian_scale_eval_expression(term$first[[j]], eval_env)
          for (k in chosen_scale_vars)
          {
            partial2 <- .gaussian_scale_eval_expression(term$second[[j]][[k]], eval_env)
            d2X[[j]][[k]][, col] <- if (identical(j, k))
              partial2 * base_deriv[[j]]^2 + partial1_j * base_second[[j]]
            else
              partial2 * base_deriv[[j]] * base_deriv[[k]]
          }
        }
      }

      if (qr(X)$rank < ncol(X))
        stop("Gaussian scale-aware conductance design matrix is rank deficient at the current sigma",
             call. = FALSE)

      linpred <- as.vector(X %*% beta)
      conductance <- exp(linpred)

      linpred_grad <- matrix(0, nrow = n, ncol = length(default))
      colnames(linpred_grad) <- names(default)
      linpred_grad[, beta_names] <- X

      if (length(chosen_scale_vars))
      {
        for (j in seq_along(chosen_scale_vars))
        {
          sigma_var <- chosen_scale_vars[j]
          linpred_grad[, sigma_names[j]] <- as.vector(dX[[sigma_var]] %*% beta)
        }
      }

      df__dtheta_matrix <- conductance * linpred_grad

      d2_linpred <- function(k, l)
      {
        par_k <- names(default)[k]
        par_l <- names(default)[l]

        beta_k <- match(par_k, beta_names, nomatch = 0L)
        beta_l <- match(par_l, beta_names, nomatch = 0L)
        sigma_k <- match(par_k, sigma_names, nomatch = 0L)
        sigma_l <- match(par_l, sigma_names, nomatch = 0L)

        if (beta_k > 0L && sigma_l > 0L)
          return(dX[[chosen_scale_vars[sigma_l]]][, beta_k])
        if (sigma_k > 0L && beta_l > 0L)
          return(dX[[chosen_scale_vars[sigma_k]]][, beta_l])
        if (sigma_k > 0L && sigma_l > 0L)
          return(as.vector(d2X[[chosen_scale_vars[sigma_k]]][[chosen_scale_vars[sigma_l]]] %*% beta))

        rep(0, n)
      }

      confint <- function(theta, vcov, quantile = 0.95,
                          scale = c("conductance", "linpred"))
      {
        scale <- match.arg(scale)
        linpred_sd <- sqrt(pmax(rowSums((linpred_grad %*% vcov) * linpred_grad), 0))
        ci <- linpred + qnorm((1 - quantile) / 2) * linpred_sd %*% t(c(1, -1))
        colnames(ci) <- c("lower", "upper")
        attr(ci, "quantile") <- quantile
        if (identical(scale, "linpred"))
          return(ci)
        exp(ci)
      }

      list(
        conductance = conductance,
        confint = confint,
        df__dx = function(k)
          stop("partial_X is not available for Gaussian scale-aware conductance",
               call. = FALSE),
        df__dtheta = function(k) df__dtheta_matrix[, k],
        df__dtheta_matrix = df__dtheta_matrix,
        d2f__dtheta_dtheta = function(k, l)
          conductance * (linpred_grad[, k] * linpred_grad[, l] + d2_linpred(k, l)),
        d2f__dtheta_dx = function(k, l)
          stop("partial_X is not available for Gaussian scale-aware conductance",
               call. = FALSE)
      )
    }

    class(conductance_model) <- c("terradish_conductance_model",
                                  "radish_conductance_model")
    attr(conductance_model, "default") <- default
    attr(conductance_model, "lower") <- lower
    attr(conductance_model, "upper") <- upper
    attr(conductance_model, "supports_partial") <- FALSE
    attr(conductance_model, "link") <- "log"
    attr(conductance_model, "scale_vars") <- chosen_scale_vars
    attr(conductance_model, "gaussian_scale") <- TRUE
    attr(conductance_model, "gaussian_scale_info") <- list(
      scale_vars = chosen_scale_vars,
      lower = lower[sigma_names],
      upper = upper[sigma_names],
      resolution = abs(terra::res(stack[[1]])),
      extent_diagonal = .gaussian_scale_extent_diagonal(stack),
      is_lonlat = terra::is.lonlat(stack),
      standardize = standardize
    )
    conductance_model
  }

  class(factory) <- c("terradish_conductance_model_factory",
                      "radish_conductance_model_factory")
  attr(factory, "default") <- NULL
  attr(factory, "preferred_optimizer") <- "bfgs"
  attr(factory, "supports_partial") <- FALSE
  attr(factory, "requires_fixed_graph") <- TRUE
  attr(factory, "link") <- "log"
  attr(factory, "gaussian_scale") <- TRUE
  factory
}

#' Summarize fitted Gaussian scales of effect
#'
#' Converts fitted Gaussian \code{sigma} parameters into map-unit, cell-based,
#' and Gaussian-kernel interpretation summaries for a fitted model produced with
#' \code{\link{gaussian_smoothed_loglinear_conductance}}.
#'
#' @param object A fitted \code{\link{terradish}} object using the Gaussian
#'   scale-aware conductance model.
#' @param probabilities Probabilities used to summarize the one-dimensional
#'   half-width and two-dimensional radial extent of the Gaussian kernel.
#' @param distance_per_map_unit Optional scalar conversion factor used to express
#'   \code{sigma} and the derived radii in user-supplied distance units. For
#'   example, if the raster CRS is in metres, use \code{distance_per_map_unit =
#'   0.001} and \code{distance_unit = "km"} to report kilometres.
#' @param distance_unit Optional label used when
#'   \code{distance_per_map_unit} is supplied.
#'
#' @details
#' The returned table always reports \code{sigma} on the raster's native map
#' scale and in raster cell widths. For each probability \code{p},
#' \code{axis_*} gives the one-dimensional half-width
#' \eqn{qnorm((1 + p) / 2) * sigma}, while \code{radial_*} gives the isotropic
#' two-dimensional radius \eqn{sigma * sqrt(-2 * log(1 - p))} containing
#' proportion \code{p} of the Gaussian kernel mass.
#'
#' If the retained raster is in longitude/latitude, the native-unit results are
#' in degrees. In that case, use a projected raster for direct distance
#' interpretation, or provide an approximate \code{distance_per_map_unit}
#' conversion if you need a quick descriptive summary.
#'
#' @return A data frame with one row per fitted \code{sigma} parameter.
#'
#' @examples
#' \dontrun{
#' fit <- terradish(melip.Fst ~ forestcover,
#'                  data = surface,
#'                  conductance_model = gaussian_smoothed_loglinear_conductance(surface),
#'                  measurement_model = leastsquares,
#'                  optimizer = "bfgs")
#' gaussian_scale_summary(fit, distance_per_map_unit = 0.001, distance_unit = "km")
#' }
#'
#' @export
gaussian_scale_summary <- function(object,
                                   probabilities = c(0.5, 0.95),
                                   distance_per_map_unit = NULL,
                                   distance_unit = NULL)
{
  stopifnot(inherits(object, c("terradish", "radish")))

  info <- attr(object$submodels$f, "gaussian_scale_info", exact = TRUE)
  if (is.null(info))
    stop("`object` was not fitted with `gaussian_smoothed_loglinear_conductance()`",
         call. = FALSE)

  theta <- coef(object)
  sigma_names <- paste0("sigma.", info$scale_vars)
  if (!length(sigma_names))
    stop("No Gaussian sigma parameters were estimated for this model",
         call. = FALSE)
  if (!all(sigma_names %in% names(theta)))
    stop("Could not locate the fitted sigma coefficients in `object`",
         call. = FALSE)

  probabilities <- c(probabilities)
  if (!length(probabilities) ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1))
    stop("`probabilities` must contain values strictly between 0 and 1",
         call. = FALSE)

  if (!is.null(distance_per_map_unit))
  {
    if (!is.numeric(distance_per_map_unit) ||
        length(distance_per_map_unit) != 1L ||
        !is.finite(distance_per_map_unit) ||
        distance_per_map_unit <= 0)
      stop("`distance_per_map_unit` must be one finite positive number",
           call. = FALSE)
    if (is.null(distance_unit))
      distance_unit <- "distance"
  }

  sigma <- theta[sigma_names]
  out <- data.frame(
    covariate = info$scale_vars,
    sigma = unname(sigma),
    sigma_lower = unname(info$lower[sigma_names]),
    sigma_upper = unname(info$upper[sigma_names]),
    native_unit = rep(if (isTRUE(info$is_lonlat)) "degrees" else "map_units",
                      length(sigma_names)),
    sigma_cells_x = unname(sigma / info$resolution[1]),
    sigma_cells_y = unname(sigma / info$resolution[2]),
    stringsAsFactors = FALSE
  )

  for (prob in probabilities)
  {
    label <- .gaussian_scale_probability_label(prob)
    out[[paste0("axis_", label)]] <- stats::qnorm((1 + prob) / 2) * out$sigma
    out[[paste0("radial_", label)]] <- sqrt(-2 * log(1 - prob)) * out$sigma
  }

  if (!is.null(distance_per_map_unit))
  {
    suffix <- .gaussian_scale_unit_suffix(distance_unit)
    out[[paste0("sigma_", suffix)]] <- out$sigma * distance_per_map_unit
    out[[paste0("sigma_lower_", suffix)]] <- out$sigma_lower * distance_per_map_unit
    out[[paste0("sigma_upper_", suffix)]] <- out$sigma_upper * distance_per_map_unit
    for (prob in probabilities)
    {
      label <- .gaussian_scale_probability_label(prob)
      out[[paste0("axis_", label, "_", suffix)]] <-
        out[[paste0("axis_", label)]] * distance_per_map_unit
      out[[paste0("radial_", label, "_", suffix)]] <-
        out[[paste0("radial_", label)]] * distance_per_map_unit
    }
  }

  if (isTRUE(info$is_lonlat) && is.null(distance_per_map_unit))
    warning("The retained raster is in longitude/latitude, so sigma is reported in degrees. ",
            "Use a projected raster or provide `distance_per_map_unit` for a descriptive conversion.",
            call. = FALSE)

  out
}
