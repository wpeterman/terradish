assemble_model_matrix <- function(formula, spdat, check_rank = TRUE)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(spdat))

  # check if formula is consistent with data, remove response, add intercept
  formula_covariates <- attr(delete.response(terms(formula)), "factors")
  if (length(formula_covariates) > 0)
  {
    # Use all.vars() so that in-line transformations such as I(x^2) or
    # interactions x:z do not appear as required column names; only the
    # underlying raw variables need to be present in the data frame.
    stopifnot(all.vars(formula) %in% colnames(spdat))
    formula <- reformulate(colnames(formula_covariates), env = environment(formula))

    # if any layers are not in formula, remove them
    missing_covariates <- !(colnames(spdat) %in% all.vars(formula))
    if (any(missing_covariates))
    {
      unused_covariates <- colnames(spdat)[missing_covariates]
      warning("Removed unused spatial covariates: ", 
              paste(unused_covariates, collapse = " "))
      spdat <- spdat[,!missing_covariates,drop=FALSE]
    }
  }
  else
    formula <- formula(~1)

  # get model matrix and check for rank deficiency
  # NOTE: sparse via Matrix::sparse.model.matrix?
  frame <- stats::model.frame(formula, data = spdat, na.action = stats::na.fail)
  spdat <- model.matrix(formula, data = frame)
  prediction_spec <- list(terms = terms(frame),
    xlevels = lapply(frame[vapply(frame, is.factor, logical(1))], levels),
    contrasts = attr(spdat, "contrasts"))
  if (check_rank) stopifnot(qr(spdat)$rank == ncol(spdat))
  if (ncol(spdat) > 1) #unless IBD, remove intercept
    spdat <- spdat[,colnames(spdat) != "(Intercept)", drop=FALSE]

  attr(spdat, "prediction_spec") <- prediction_spec
  spdat
}

.validate_conductance_values <- function(conductance, context = "Conductance model")
{
  if (!all(is.finite(conductance)))
    stop(context, " produced non-finite conductance values at the current parameters.",
         call. = FALSE)
  if (!all(conductance > 0))
    stop(context, " requires strictly positive conductance values at the current parameters.",
         call. = FALSE)
  conductance
}

.smooth_loglinear_eval_arg <- function(arg, default, envir)
{
  if (is.null(arg))
    return(default)
  eval(arg, envir = envir)
}

.smooth_loglinear_basis <- function(label, x, default_df, default_basis,
                                    default_degree, default_intercept,
                                    envir, spec = NULL)
{
  if (is.null(spec))
  {
    call <- str2lang(label)
    if (!is.call(call) || !identical(call[[1]], as.name("s")) || length(call) < 2)
      stop("Smooth terms must be written as s(variable, ...).", call. = FALSE)

    variable <- all.vars(call[[2]])
    if (length(variable) != 1L || !identical(variable, as.character(call[[2]])))
      stop("Smooth terms currently support one raw column name, e.g. s(elevation).",
           call. = FALSE)

    args <- as.list(call[-1])
    df <- .smooth_loglinear_eval_arg(args$df, NULL, envir)
    if (is.null(df))
      df <- .smooth_loglinear_eval_arg(args$k, default_df, envir)
    basis <- .smooth_loglinear_eval_arg(args$basis, NULL, envir)
    if (is.null(basis))
      basis <- .smooth_loglinear_eval_arg(args$bs, default_basis, envir)
    degree <- .smooth_loglinear_eval_arg(args$degree, default_degree, envir)
    intercept <- .smooth_loglinear_eval_arg(args$intercept,
                                            default_intercept, envir)

    df <- as.integer(df)
    if (length(df) != 1L || is.na(df) || df < 1L)
      stop("Smooth term df/k must be a positive integer.", call. = FALSE)
    basis <- match.arg(as.character(basis), c("ns", "bs"))
    degree <- as.integer(degree)
    intercept <- isTRUE(intercept)
  }
  else
  {
    variable <- spec$variable
    df <- spec$df
    basis <- spec$basis
    degree <- spec$degree
    intercept <- spec$intercept
  }

  if (!variable %in% colnames(x))
    stop("Smooth variable not found in x: ", variable, call. = FALSE)

  z <- x[[variable]]
  if (!is.numeric(z))
    stop("Smooth variable must be numeric: ", variable, call. = FALSE)

  out <- if (is.null(spec))
  {
    switch(
      basis,
      ns = splines::ns(z, df = df, intercept = intercept),
      bs = splines::bs(z, df = df, degree = degree, intercept = intercept)
    )
  }
  else
  {
    switch(
      basis,
      ns = splines::ns(z, knots = spec$knots,
                       Boundary.knots = spec$Boundary.knots,
                       intercept = intercept),
      bs = splines::bs(z, knots = spec$knots,
                       Boundary.knots = spec$Boundary.knots,
                       degree = degree, intercept = intercept)
    )
  }
  out <- as.matrix(out)
  # Preserve the fitted basis and its graph-wide origin when rebuilding on
  # focal support, plotting grids, or a coarser raster.
  knots <- attr(out, "knots", exact = TRUE)
  boundary_knots <- attr(out, "Boundary.knots", exact = TRUE)
  centers <- if (is.null(spec$centers)) colMeans(out) else spec$centers
  out <- sweep(out, 2L, centers, "-")
  if (!is.null(spec$columns) && length(spec$columns) == ncol(out))
    colnames(out) <- spec$columns
  else
    colnames(out) <- paste0("s(", variable, ").", seq_len(ncol(out)))
  attr(out, "smooth_spec") <- list(
    label = label,
    variable = variable,
    df = df,
    basis = basis,
    degree = degree,
    intercept = intercept,
    knots = knots,
    Boundary.knots = boundary_knots,
    centers = centers,
    columns = colnames(out)
  )
  out
}

.smooth_loglinear_model_matrix <- function(formula, x, df, basis, degree,
                                           intercept, smooth_specs = NULL,
                                           param_spec = NULL)
{
  fitting_basis <- is.null(smooth_specs)
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(x))

  tt <- terms(formula, specials = "s", keep.order = TRUE)
  labels <- attr(tt, "term.labels")
  smooth_idx <- attr(tt, "specials")$s
  smooth_labels <- labels[smooth_idx]
  param_labels <- labels[setdiff(seq_along(labels), smooth_idx)]
  if (!is.null(smooth_specs) && length(smooth_specs) != length(smooth_labels))
    stop("Stored smooth specifications do not match the formula smooth terms.",
         call. = FALSE)

  if (length(param_labels))
  {
    param_formula <- reformulate(param_labels, env = environment(formula))
    param_vars <- all.vars(param_formula)
    if (is.null(param_spec)) {
      param_x <- assemble_model_matrix(param_formula, x[, param_vars, drop = FALSE],
                                       check_rank = fitting_basis)
      param_spec <- attr(param_x, "prediction_spec", exact = TRUE)
      param_spec$columns <- colnames(param_x)
    } else {
      frame <- stats::model.frame(param_spec$terms, x, xlev = param_spec$xlevels,
                                   na.action = stats::na.fail)
      param_x <- model.matrix(param_spec$terms, frame, contrasts.arg = param_spec$contrasts)
      param_x <- param_x[, param_spec$columns, drop = FALSE]
    }
  }
  else if (!length(smooth_labels))
    param_x <- assemble_model_matrix(~1, x)
  else
    param_x <- matrix(numeric(0), nrow = nrow(x), ncol = 0L)

  smooth_x <- lapply(
    seq_along(smooth_labels),
    function(i) {
      .smooth_loglinear_basis(
        smooth_labels[[i]],
        x = x,
        default_df = df,
        default_basis = basis,
        default_degree = degree,
        default_intercept = intercept,
        envir = environment(formula),
        spec = if (is.null(smooth_specs)) NULL else smooth_specs[[i]]
      )
    }
  )
  smooth_specs <- lapply(smooth_x, attr, "smooth_spec", exact = TRUE)
  smooth_x <- if (length(smooth_x)) do.call(cbind, smooth_x) else NULL

  out <- cbind(param_x, smooth_x)
  if (!ncol(out))
    stop("No conductance covariates were produced from formula.", call. = FALSE)
  if (fitting_basis && qr(out)$rank < ncol(out))
    stop("Smooth conductance model matrix is rank deficient.", call. = FALSE)
  rownames(out) <- NULL
  attr(out, "smooth_specs") <- smooth_specs
  attr(out, "param_spec") <- param_spec
  out
}

.smooth_loglinear_conductance_from_matrix <- function(x, df, basis, degree,
                                                      intercept, smooth_specs)
{
  default <- rep(0, ncol(x))
  names(default) <- colnames(x)

  conductance_model <- function(theta)
  {
    stopifnot(length(theta) == ncol(x))

    conductance <- as.vector(exp(x %*% theta))
    conductance <- .validate_conductance_values(
      conductance,
      context = "smooth_loglinear_conductance()"
    )
    df__dtheta_matrix <- conductance * x

    df__dx             <- function(k)    c(conductance * theta[k])
    df__dtheta         <- function(k)    c(df__dtheta_matrix[, k])
    d2f__dtheta_dtheta <- function(k, l) conductance * x[, k] * x[, l]
    d2f__dtheta_dx     <- function(k, l) conductance * ((k == l) + x[, k] * theta[l])

    confint <- function(theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
    {
      scale <- match.arg(scale)
      cond_sd <- sqrt(rowSums((x %*% vcov) * x))
      ci <- log(conductance) + qnorm((1 - quantile)/2) * cond_sd %*% t(c(1, -1))
      colnames(ci) <- c("lower", "upper")
      attr(ci, "quantile") <- quantile
      if (scale == "linpred")
        return (ci)
      else if (scale == "conductance")
        return (exp(ci))
    }

    list(conductance        = conductance,
         confint            = confint,
         df__dx             = df__dx,
         df__dtheta         = df__dtheta,
         df__dtheta_matrix  = df__dtheta_matrix,
         d2f__dtheta_dtheta = d2f__dtheta_dtheta,
         d2f__dtheta_dx     = d2f__dtheta_dx)
  }

  class(conductance_model) <- c("terradish_conductance_model",
                                "radish_conductance_model")
  attr(conductance_model, "default") <- default
  attr(conductance_model, "link") <- "log"
  attr(conductance_model, "smooth_loglinear") <- TRUE
  attr(conductance_model, "smooth_loglinear_info") <- list(
    df = df,
    basis = basis,
    degree = degree,
    intercept = intercept,
    columns = names(default),
    smooth_specs = smooth_specs
  )
  attr(conductance_model, "plot_factory") <- .smooth_loglinear_factory(
    df = df,
    basis = basis,
    degree = degree,
    intercept = intercept,
    smooth_specs = smooth_specs,
    param_spec = attr(x, "param_spec", exact = TRUE)
  )
  conductance_model
}

#' Conductance model factories
#'
#' Functions that generate objects of class \code{"terradish_conductance_model"}
#' that represent mappings from spatial data (e.g. rasters) to conductance.
#'
#' @name terradish_conductance_model_factory
#' @seealso \code{\link{loglinear_conductance}},
#'   \code{\link{smooth_loglinear_conductance}}
terradish_conductance_model_factory <- NULL

#' Legacy radish conductance-model class alias
#'
#' The legacy \code{"radish_conductance_model_factory"} class name is retained
#' for backward compatibility with older objects and workflows.
#'
#' @name radish_conductance_model_factory
#' @keywords internal
NULL

#' Log-link conductance model
#'
#' Returns a function of class \code{"terradish_conductance_model"} that
#' represents a log-linear mapping from spatial covariates to conductance.
#' This is the recommended conductance model for most applications.
#'
#' @param formula Model formula describing which spatial covariates drive
#'   conductance. The left-hand side is ignored; only the right-hand side terms
#'   are used (e.g. \code{~ altitude + forestcover}).
#' @param x Data frame of spatial covariates extracted from a
#'   \code{\link{conductance_surface}} object (typically \code{surface$x}).
#'   Every variable named on the right-hand side of \code{formula} must be a
#'   column of \code{x}.
#'
#' @details
#' The conductance at grid cell \code{i} is:
#'
#' \deqn{C_i = \exp(\theta_1 x_{i1} + \theta_2 x_{i2} + \ldots)}
#'
#' where \eqn{x_{ij}} is the value of covariate \eqn{j} at cell \eqn{i} and
#' \eqn{\theta_j} is the corresponding conductance parameter.
#'
#' The intercept is intentionally omitted: multiplying all conductances by a
#' common factor rescales resistance inversely, and the measurement model's
#' resistance coefficient absorbs that factor. The absolute conductance level
#' is therefore not identifiable.
#'
#' \strong{Interpreting \eqn{\theta}:}
#' \itemize{
#'   \item \eqn{\theta_j > 0}: higher values of covariate \eqn{j} increase
#'     fitted relative conductance and tend to lower graph resistance distance,
#'     conditional on the other terms and graph domain.
#'   \item \eqn{\theta_j < 0}: higher values lower fitted relative conductance.
#'   \item \eqn{\theta_j = 0}: the fitted conductance is locally unchanged by
#'     that term, conditional on the model.
#'   \item A one-standard-deviation increase in covariate \eqn{j} multiplies
#'     fitted conductance by \eqn{\exp(\theta_j)} when the covariate was
#'     standardized and the term is linear.
#' }
#' These are model-based associations, not direct estimates of movement rate,
#' habitat suitability, or a causal landscape effect.
#'
#' The exponential link guarantees strictly positive conductances for any real
#' \eqn{\theta}.
#'
#' Categorical covariates must be stored as \code{factor} columns in \code{x}
#' (see \code{\link{conductance_surface}} for how to encode them). They are
#' dummy-coded using the default R contrasts via
#' \code{\link[stats]{model.matrix}}, with one level dropped as a reference.
#' In-formula transformations such as \code{I(x^2)} and interaction terms
#' \code{x * z} are supported.
#'
#' @return A function of class \code{"terradish_conductance_model"} that
#'   accepts a numeric vector of conductance parameters \code{theta} and
#'   returns a list with elements \code{conductance} (a vector of per-cell
#'   conductance values), \code{confint} (a function for confidence intervals),
#'   and derivative functions used internally by the optimizer.
#'
#' @seealso \code{\link{gaussian_smoothed_loglinear_conductance}},
#'   \code{\link{conductance_surface}}, \code{\link{terradish}}
#'
#' @examples
#' x <- data.frame(altitude = c(-1, 0, 1), forestcover = c(0.2, 0.6, 0.5))
#' model <- loglinear_conductance(~ altitude + forestcover, x)
#'
#' # Evaluate at specific parameter values
#' fit <- model(c(altitude = 0.3, forestcover = -0.2))
#' fit$conductance  # per-cell conductance values
#'
#' # Interaction and polynomial terms work too
#' x2 <- data.frame(altitude = c(-1, -0.25, 0.5, 1),
#'                  fc = c(0.2, 0.6, 0.5, 0.8))
#' model2 <- loglinear_conductance(~ altitude + I(altitude^2) + fc, x2)
#'
#' @export

loglinear_conductance <- function(formula, x)
{
  .loglinear_conductance_from_matrix(assemble_model_matrix(formula, x))
}

.loglinear_conductance_from_matrix <- function(x)
{

  # default starting values
  default <- rep(0, ncol(x))
  names(default) <- colnames(x)

  conductance_model <- function(theta)
  {
    stopifnot(length(theta) == ncol(x))

    conductance        <- as.vector(exp(x %*% theta))
    conductance        <- .validate_conductance_values(
      conductance,
      context = "loglinear_conductance()"
    )
    df__dtheta_matrix  <- conductance * x

    # first- and second-order derivatives
    df__dx             <- function(k)    conductance * theta[k]
    df__dtheta         <- function(k)    df__dtheta_matrix[, k]
    d2f__dtheta_dtheta <- function(k, l) conductance * x[,k] * x[,l]
    d2f__dtheta_dx     <- function(k, l) conductance * ((k==l) + x[,k] * theta[l])

    # asymptotic confidence intervals
    confint <- function(theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
    {
      scale <- match.arg(scale)
      cond_sd <- sqrt(rowSums((x %*% vcov) * x))
      ci <- log(conductance) + qnorm((1 - quantile)/2) * cond_sd %*% t(c(1, -1))
      colnames(ci) <- c("lower", "upper")
      attr(ci, "quantile") <- quantile 
      if (scale == "linpred") 
        return (ci)
      else if (scale == "conductance")
        return (exp(ci))
    }

    # default starting values
    default <- function()
    {
      out <- rep(0, ncol(x))
      names(out) <- colnames(x)
      out
    }

    list(conductance        = conductance,
         confint            = confint,
         df__dx             = df__dx,
         df__dtheta         = df__dtheta,
         df__dtheta_matrix  = df__dtheta_matrix,
         d2f__dtheta_dtheta = d2f__dtheta_dtheta, 
         d2f__dtheta_dx     = d2f__dtheta_dx)
  }

  class(conductance_model) <- c("terradish_conductance_model",
                                "radish_conductance_model")
  attr(conductance_model, "default") <- default
  attr(conductance_model, "prediction_spec") <- attr(x, "prediction_spec", exact = TRUE)
  conductance_model
}
class(loglinear_conductance) <- c("terradish_conductance_model_factory",
                                  "radish_conductance_model_factory")

#' Smooth log-link conductance model
#'
#' Returns a function of class \code{"terradish_conductance_model"} that
#' represents a log-linear conductance surface after expanding selected
#' numeric covariates into spline basis columns. The result is an unpenalized
#' regression spline with a fixed number of coefficients, fitted through the
#' selected terradish measurement likelihood.
#'
#' @param formula Model formula describing which spatial covariates drive
#'   conductance. Smooth terms are written as \code{s(variable, df = 4)} or
#'   \code{s(variable, k = 4)}. Ordinary linear terms can be mixed with smooth
#'   terms, e.g. \code{~ forestcover + s(altitude, df = 4)}.
#' @param x Data frame of spatial covariates extracted from a
#'   \code{\link{conductance_surface}} object (typically \code{surface$x}).
#' @param df Default degrees of freedom for smooth terms that do not supply
#'   \code{df} or \code{k}.
#' @param basis Default spline basis. \code{"ns"} uses
#'   \code{\link[splines]{ns}} natural splines; \code{"bs"} uses
#'   \code{\link[splines]{bs}} B-splines.
#' @param degree B-spline polynomial degree used when \code{basis = "bs"}.
#' @param intercept Logical; include intercept columns within each spline
#'   basis? The default is \code{FALSE}, matching terradish's usual convention
#'   of omitting a global conductance intercept.
#'
#' @details
#' The model first expands all \code{s()} terms into basis columns and then
#' applies the same positive log-link used by \code{\link{loglinear_conductance}}:
#'
#' \deqn{C_i = \exp(B_i \theta)}
#'
#' Each spline column is centered over the active graph cells. The fitted
#' knots, boundary knots, degree, and column centers are retained and reused
#' for plotting, focal-support prediction, and coarse warm starts. Centering
#' changes the arbitrary conductance reference level but preserves the
#' likelihood and fitted conductance coefficients. Uncertainty bands are
#' relative to the landscape-mean spline contribution to log conductance.
#'
#' where \eqn{B_i} is the row of the expanded spline/parametric design matrix
#' for grid cell \eqn{i}. The same conductance factory works with MLPE and
#' Wishart measurement models.
#'
#' Smoothing parameters are not estimated in this first implementation. The
#' effective smoothness is controlled by fixed basis dimension through
#' \code{df} or \code{k}; use small values for stable first-pass model
#' comparison.
#' AIC and AICc at nominal Wishart information can over-select flexible
#' curves. Prefer spatial CV for selecting conductance terms, and examine
#' sensitivity to \code{nu} before interpreting uncertainty. Tails beyond
#' the sampled covariate range are poorly identified. A selected curve can
#' describe how a process maps onto the graph without identifying an
#' ecological response mechanism.
#'
#' \code{summary(fit)$spline_monotonicity} describes each fitted smooth over
#' the range at focal sites. It reports whether the curve is monotone, its
#' direction, and the number of derivative sign changes. This is a shape
#' diagnostic, not an uncertainty test.
#'
#' @return A function of class \code{"terradish_conductance_model"} that
#'   accepts a numeric vector of conductance parameters \code{theta} and
#'   returns conductance values plus derivative functions used internally by
#'   the optimizer. Derivatives are with respect to the expanded basis columns.
#'
#' @seealso \code{\link{loglinear_conductance}}, \code{\link{mlpe}},
#'   \code{\link{conductance_surface}}, \code{\link{terradish}}
#'
#' @examples
#' x <- data.frame(altitude = seq(-1, 1, length.out = 8),
#'                 forestcover = c(0.2, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2, 0.1))
#' model <- smooth_loglinear_conductance(~ forestcover + s(altitude, df = 3), x)
#' theta <- attr(model, "default")
#' fit <- model(theta)
#' fit$conductance
#'
#' @export

smooth_loglinear_conductance <- function(formula, x, df = 4L,
                                         basis = c("ns", "bs"),
                                         degree = 3L,
                                         intercept = FALSE)
{
  basis <- match.arg(basis)
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)
  x <- .smooth_loglinear_model_matrix(formula, x, df = df, basis = basis,
                                      degree = degree, intercept = intercept)
  smooth_specs <- attr(x, "smooth_specs", exact = TRUE)
  .smooth_loglinear_conductance_from_matrix(
    x = x, df = df, basis = basis, degree = degree,
    intercept = intercept, smooth_specs = smooth_specs
  )
}
class(smooth_loglinear_conductance) <- c("terradish_conductance_model_factory",
                                         "radish_conductance_model_factory")
attr(smooth_loglinear_conductance, "link") <- "log"

.smooth_loglinear_factory <- function(df = 4L, basis = c("ns", "bs"),
                                      degree = 3L, intercept = FALSE,
                                      smooth_specs = NULL, param_spec = NULL)
{
  basis <- match.arg(basis)
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)

  factory <- function(formula, x)
  {
    x <- .smooth_loglinear_model_matrix(
      formula, x, df = df, basis = basis, degree = degree,
      intercept = intercept, smooth_specs = smooth_specs, param_spec = param_spec
    )
    .smooth_loglinear_conductance_from_matrix(
      x = x, df = df, basis = basis, degree = degree,
      intercept = intercept,
      smooth_specs = attr(x, "smooth_specs", exact = TRUE)
    )
  }

  class(factory) <- c("terradish_conductance_model_factory",
                      "radish_conductance_model_factory")
  attr(factory, "link") <- "log"
  attr(factory, "smooth_loglinear") <- TRUE
  attr(factory, "smooth_loglinear_info") <- list(
    df = df,
    basis = basis,
    degree = degree,
    intercept = intercept
  )
  factory
}
