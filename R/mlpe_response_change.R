#' Response-scale change for MLPE pairwise covariates
#'
#' Computes the expected change in the genetic-distance response associated
#' with changing one or more pairwise covariates in an
#' \code{\link{mlpe_covariates}} measurement model.  This is useful for
#' interpreting IBE coefficients, which are profiled nuisance parameters rather
#' than conductance-surface coefficients.
#'
#' @param object A fitted \code{terradish} model whose measurement model was
#'   created with \code{\link{mlpe_covariates}}.
#' @param covariate Optional character vector naming the pairwise covariates to
#'   summarize.  Defaults to all pairwise covariates stored in the measurement
#'   model.
#' @param probs Numeric vector of length two giving the lower and upper
#'   quantiles used to define the contrast when \code{values} is \code{NULL}.
#'   Defaults to the 10th and 90th percentiles.
#' @param values Optional numeric vector of length two giving explicit lower
#'   and upper covariate values for a single \code{covariate}.  Alternatively,
#'   provide a named list with one length-two numeric vector per covariate.
#' @param conf.level Confidence level for the Wald confidence interval.
#'
#' @details
#' For a model fitted with \code{mlpe_covariates()}, the pairwise mean structure
#' includes terms of the form \eqn{\gamma Z_{ij}}.  This function reports
#' \eqn{\gamma (z_{high} - z_{low})}, the expected response-scale change in
#' genetic distance across the requested covariate contrast, conditional on the
#' fitted conductance surface.
#'
#' The interval uses the conditional standard error for the MLPE covariate
#' coefficient from \code{summary(object)} and treats the selected covariate
#' contrast as fixed.
#'
#' @return A data frame with one row per covariate and columns describing the
#'   low and high covariate values, the contrast, estimated response-scale
#'   change, confidence limits, and confidence level.
#'
#' @seealso \code{\link{mlpe_covariates}},
#'   \code{\link{pairwise_endpoint_covariates}}, \code{\link{terradish}}
#'
#' @examples
#' \dontrun{
#' change <- mlpe_response_change(fit_joint)
#' change
#'
#' mlpe_response_change(fit_joint, covariate = "absdiff_altitude",
#'                      probs = c(0.25, 0.75))
#' }
#'
#' @export
mlpe_response_change <- function(object,
                                 covariate = NULL,
                                 probs = c(0.1, 0.9),
                                 values = NULL,
                                 conf.level = 0.95)
{
  if (!inherits(object, c("terradish", "radish")))
    stop("`object` must be a fitted terradish model")

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be a single number between 0 and 1")

  pairwise_covariates <- attr(object$submodels$g, "pairwise_covariates",
                              exact = TRUE)
  if (is.null(pairwise_covariates))
    stop("`object` does not contain pairwise covariates from `mlpe_covariates()`")

  zmat <- .as_pairwise_covariate_matrix(pairwise_covariates)
  if (is.null(colnames(zmat)))
    stop("pairwise covariates must have column names")

  sm <- summary(object, conf.level = conf.level)
  if (is.null(sm$phi_table))
    stop("`object` does not contain a nuisance-parameter summary table")

  available <- intersect(colnames(zmat), rownames(sm$phi_table))
  if (!length(available))
    stop("no pairwise covariate coefficients were found in `summary(object)$phi_table`")

  if (is.null(covariate))
    covariate <- available
  else
  {
    missing <- setdiff(covariate, available)
    if (length(missing))
      stop("unknown pairwise covariate(s): ", paste(missing, collapse = ", "))
  }

  if (!is.null(values) && length(covariate) > 1L && !is.list(values))
    stop("`values` must be a named list when summarizing multiple covariates")

  if (is.null(values))
  {
    if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
        any(probs < 0 | probs > 1) || probs[1] >= probs[2])
      stop("`probs` must contain two increasing probabilities in [0, 1]")
  }

  zcrit <- qnorm((1 + conf.level) / 2)

  rows <- lapply(covariate, function(nm) {
    vals <- .mlpe_response_change_values(nm, zmat, values, probs)
    contrast <- vals[2] - vals[1]
    estimate <- sm$phi_table[nm, "Estimate"] * contrast
    se <- sm$phi_table[nm, "Std. Error"] * abs(contrast)

    data.frame(
      covariate = nm,
      low = vals[1],
      high = vals[2],
      contrast = contrast,
      estimate = estimate,
      conf.low = estimate - zcrit * se,
      conf.high = estimate + zcrit * se,
      conf.level = conf.level,
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

.mlpe_response_change_values <- function(covariate, zmat, values, probs)
{
  if (is.null(values))
    return(as.numeric(stats::quantile(zmat[, covariate], probs = probs,
                                      names = FALSE, na.rm = TRUE)))

  vals <- if (is.list(values))
  {
    if (is.null(names(values)) || !covariate %in% names(values))
      stop("`values` must contain an entry named ", covariate)
    values[[covariate]]
  }
  else
  {
    values
  }

  if (!is.numeric(vals) || length(vals) != 2L || anyNA(vals))
    stop("each `values` entry must be a length-two numeric vector")
  as.numeric(vals)
}
