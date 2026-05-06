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

#' Conductance model factories
#'
#' Functions that generate objects of class \code{"terradish_conductance_model"}
#' that represent mappings from spatial data (e.g. rasters) to conductance.
#'
#' Legacy \code{"radish_*"} class names are retained for backward
#' compatibility.
#'
#' @name terradish_conductance_model_factory
#' @aliases radish_conductance_model_factory
#' @seealso \code{\link{linear_conductance}}, \code{\link{loglinear_conductance}}
terradish_conductance_model_factory <- NULL

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
#' x2 <- data.frame(altitude = c(-1, 0, 1), fc = c(0.2, 0.6, 0.5))
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

    stopifnot(all(conductance > 0))
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

    stopifnot(all(conductance > 0))

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
