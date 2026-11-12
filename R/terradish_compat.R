#' Deprecated radish names
#'
#' \pkg{terradish} grew out of the \pkg{radish} package. The original function
#' and class names still work so that older scripts and saved fits keep
#' running, but they are deprecated. New code should use the
#' \code{terradish} names in the table below.
#'
#' @details
#' \tabular{ll}{
#'   \strong{Legacy name} \tab \strong{Use instead} \cr
#'   \code{radish()} \tab \code{\link{terradish}()} \cr
#'   \code{radish_grid()} \tab \code{\link{terradish_grid}()} \cr
#'   \code{radish_distance()} \tab \code{\link{terradish_distance}()} \cr
#'   \code{radish_algorithm()} \tab \code{\link{terradish_algorithm}()} \cr
#'   \code{terra_radish_scale_optim()} \tab \code{\link{terradish_scale_optim}()} \cr
#'   methods for class \code{"radish"} (\code{print()}, \code{summary()},
#'     \code{coef()}, \code{fitted()}, \code{residuals()}, \code{simulate()},
#'     \code{anova()}, \code{logLik()}, \code{AIC()}, \code{plot()}) \tab
#'     the same methods for class \code{"terradish"}; see
#'     \code{\link{terradish_methods}} and \code{\link{plot.terradish}} \cr
#'   \code{conductance()} for class \code{"radish_graph"} \tab
#'     \code{\link{conductance}()} for class \code{"terradish_graph"} \cr
#'   class \code{"radish_conductance_model_factory"} \tab
#'     \code{\link{terradish_conductance_model_factory}} \cr
#'   class \code{"radish_measurement_model"} \tab
#'     \code{\link{terradish_measurement_model}} \cr
#' }
#'
#' The wrappers \code{radish()}, \code{radish_grid()},
#' \code{radish_distance()}, \code{radish_algorithm()}, and
#' \code{terra_radish_scale_optim()} pass their arguments unchanged to the
#' replacement and signal a deprecation warning. Objects created by \pkg{terradish}
#' carry both the new and the legacy class names, so either set of methods
#' dispatches on them.
#'
#' Arguments, defaults, and returned values are documented on the
#' replacement pages. Full compatibility with every historical \pkg{radish}
#' entry point is not promised.
#'
#' @name radish_compatibility
#' @seealso \code{\link{terradish}}, \code{\link{terradish-package}}
NULL

.terradish_deprecate <- function(old, new)
{
  .Deprecated(new = new, package = "terradish", old = old)
}

.terradish_forward_call <- function(call, new)
{
  call[[1L]] <- as.name(new)
  eval(call, parent.frame(2L))
}

.terradish_set_class <- function(x, ...)
{
  class(x) <- unique(c(..., class(x)))
  x
}
