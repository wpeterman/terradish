assemble_model_matrix <- function(formula, spdat)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(spdat))

  # check if formula is consistent with data, remove response, add intercept
  formula_covariates <- attr(delete.response(terms(formula)), "factors")
  if (length(formula_covariates) > 0)
  {
    # Use all.vars() so that in-line transformations such as I(x^2) or
    # interactions x:z do not appear as required column names — only the
    # underlying raw variables need to be present in the data frame.
    stopifnot(all.vars(formula) %in% colnames(spdat))
    formula <- reformulate(colnames(formula_covariates))

    # if any layers are not in formula, remove them
    missing_covariates <- !(colnames(spdat) %in% rownames(formula_covariates))
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
  spdat <- model.matrix(formula, data = spdat)
  stopifnot(qr(spdat)$rank == ncol(spdat))
  if (ncol(spdat) > 1) #unless IBD, remove intercept
    spdat <- spdat[,colnames(spdat) != "(Intercept)", drop=FALSE]

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
    knots = attr(out, "knots", exact = TRUE),
    Boundary.knots = attr(out, "Boundary.knots", exact = TRUE),
    columns = colnames(out)
  )
  out
}

.smooth_loglinear_model_matrix <- function(formula, x, df, basis, degree,
                                           intercept, smooth_specs = NULL)
{
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
    param_formula <- reformulate(param_labels)
    param_vars <- all.vars(param_formula)
    param_x <- assemble_model_matrix(param_formula, x[, param_vars, drop = FALSE])
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
  if (qr(out)$rank < ncol(out))
    stop("Smooth conductance model matrix is rank deficient.", call. = FALSE)
  rownames(out) <- NULL
  attr(out, "smooth_specs") <- smooth_specs
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
    smooth_specs = smooth_specs
  )
  conductance_model
}

#' Conductance model factories
#'
#' Functions that generate objects of class \code{"terradish_conductance_model"}
#' that represent mappings from spatial data (e.g. rasters) to conductance.
#'
#' @name terradish_conductance_model_factory
#' @seealso \code{\link{linear_conductance}}, \code{\link{loglinear_conductance}},
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
#' constant does not change the effective resistance distances, so an intercept
#' is non-identifiable.
#'
#' \strong{Interpreting \eqn{\theta}:}
#' \itemize{
#'   \item \eqn{\theta_j > 0}: higher values of covariate \eqn{j} increase
#'     conductance (easier movement, lower resistance distance).
#'   \item \eqn{\theta_j < 0}: higher values act as a barrier.
#'   \item \eqn{\theta_j = 0}: the covariate has no effect on conductance.
#'   \item A one-standard-deviation increase in covariate \eqn{j} multiplies
#'     conductance by \eqn{\exp(\theta_j)}.
#' }
#'
#' The exponential link guarantees strictly positive conductances for any real
#' \eqn{\theta}, making \code{loglinear_conductance} more numerically stable
#' than \code{\link{linear_conductance}} when parameters stray far from zero.
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
#' @seealso \code{\link{linear_conductance}},
#'   \code{\link{gaussian_smoothed_loglinear_conductance}},
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
  x <- assemble_model_matrix(formula, x)

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
  conductance_model
}
class(loglinear_conductance) <- c("terradish_conductance_model_factory",
                                  "radish_conductance_model_factory")

#' Smooth log-link conductance model
#'
#' Returns a function of class \code{"terradish_conductance_model"} that
#' represents a log-linear conductance surface after expanding selected
#' numeric covariates into spline basis columns. This provides a GAM-like
#' conductance model while retaining terradish's existing likelihood and MLPE
#' measurement machinery.
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
#' where \eqn{B_i} is the row of the expanded spline/parametric design matrix
#' for grid cell \eqn{i}. This is GAM-like in the conductance surface, not a
#' replacement for the MLPE response model. The \code{\link{mlpe}} measurement
#' model can be used with this conductance factory in the same way it is used
#' with \code{\link{loglinear_conductance}}.
#'
#' Smoothing parameters are not estimated in this first implementation. The
#' effective smoothness is controlled by fixed basis dimension through
#' \code{df} or \code{k}; use small values for stable first-pass model
#' comparison.
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
                                      smooth_specs = NULL)
{
  basis <- match.arg(basis)
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)

  factory <- function(formula, x)
  {
    x <- .smooth_loglinear_model_matrix(
      formula, x, df = df, basis = basis, degree = degree,
      intercept = intercept, smooth_specs = smooth_specs
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

#' Identity-link conductance model
#'
#' Returns a function of class \code{"terradish_conductance_model"} that
#' represents a linear mapping from spatial covariates to conductance.
#'
#' @param formula Model formula describing which spatial covariates drive
#'   conductance. The left-hand side is ignored; only the right-hand side terms
#'   are used.
#' @param x Data frame of spatial covariates extracted from a
#'   \code{\link{conductance_surface}} object (typically \code{surface$x}).
#'
#' @details
#' The conductance at grid cell \code{i} is:
#'
#' \deqn{C_i = \theta_1 x_{i1} + \theta_2 x_{i2} + \ldots}
#'
#' The intercept is omitted because it is non-identifiable (multiplying all
#' conductances by a constant leaves resistance distances unchanged).
#'
#' \strong{When to prefer \code{linear_conductance} over
#' \code{loglinear_conductance}:}
#' \itemize{
#'   \item When you want conductance to be a direct, additive mixture of
#'     raster layers (e.g. habitat suitability scores that already live on a
#'     natural additive scale).
#'   \item When theory predicts a linear relationship.
#' }
#'
#' \strong{Caution:} conductance must be strictly positive. The optimizer does
#' not automatically enforce this; choose starting values and parameter bounds
#' so that \eqn{C_i > 0} throughout the optimization. Fitting is safer when
#' covariates are non-negative and parameters are constrained to be positive.
#' For unrestricted parameters, \code{\link{loglinear_conductance}} is more
#' numerically robust.
#'
#' Default starting values are all 1 (rather than 0 as in
#' \code{loglinear_conductance}) to ensure positive conductances at the start.
#'
#' Categorical covariates and in-formula transformations are supported via
#' \code{\link[stats]{model.matrix}}, the same as in
#' \code{\link{loglinear_conductance}}.
#'
#' @return A function of class \code{"terradish_conductance_model"} that
#'   accepts a numeric vector of conductance parameters \code{theta} and
#'   returns a list with elements \code{conductance}, \code{confint}, and
#'   internal derivative functions.
#'
#' @seealso \code{\link{loglinear_conductance}}, \code{\link{terradish}}
#'
#' @examples
#' x <- data.frame(altitude = c(1, 2, 3), forestcover = c(2, 5, 4))
#' model <- linear_conductance(~ altitude + forestcover, x)
#' fit <- model(c(altitude = 0.5, forestcover = 1))
#' fit$conductance
#'
#' @export

linear_conductance <- function(formula, x)
{
  x <- assemble_model_matrix(formula, x)

  # default starting values
  default <- rep(1, ncol(x))
  names(default) <- colnames(x)

  conductance_model <- function(theta)
  {
    stopifnot(length(theta) == ncol(x))

    conductance        <- as.vector(x %*% theta)
    conductance        <- .validate_conductance_values(
      conductance,
      context = "linear_conductance()"
    )

    ones <- matrix(1, nrow(x), 1)
    df__dtheta_matrix <- x

    # asymptotic confidence intervals
    confint <- function(theta, vcov, quantile = 0.95, scale = c("conductance", "linpred"))
    {
      scale <- match.arg(scale)
      cond_sd <- sqrt(rowSums((x %*% vcov) * x))
      ci <- conductance + qnorm((1 - quantile)/2) * cond_sd %*% t(c(1, -1))
      colnames(ci) <- c("lower", "upper")
      attr(ci, "quantile") <- quantile 
      if (scale == "linpred") 
        return (ci)
      else if (scale == "conductance")
        return (ci)
    }

    # first- and second-order derivatives
    df__dx             <- function(k)    ones * theta[k]
    df__dtheta         <- function(k)    df__dtheta_matrix[, k]
    d2f__dtheta_dtheta <- function(k, l) 0. * ones
    d2f__dtheta_dx     <- function(k, l) (k==l) * ones

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
  conductance_model
}
class(linear_conductance) <- c("terradish_conductance_model_factory",
                               "radish_conductance_model_factory")
