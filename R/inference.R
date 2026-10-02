#' Covariance and intervals for fitted conductance parameters
#'
#' Extracts observed-curvature uncertainty from a full or slim fit.
#' @param object A fitted \code{terradish} object.
#' @param parm Names or indices of parameters to include.
#' @param level Confidence level between zero and one.
#' @param ... Unused arguments for generic compatibility.
#' @details
#' Gaussian scale intervals are Wald intervals truncated at the fitted bounds.
#' Stored profile intervals take precedence when their confidence level matches.
#' The observation count is the number of focal sites (individuals or
#' populations), the sample size that AICc and BIC use in
#' \code{\link{aic_table}}. Pairwise observations are not counted, because
#' pairs that share a site are not independent.
#' @return \code{vcov()} returns a covariance matrix; \code{confint()} returns
#' a two-column interval matrix; \code{nobs()} returns the number of focal sites.
#' @name terradish_inference
#' @seealso \code{\link{terradish}}, \code{\link{slim_terradish}},
#'   \code{\link{gaussian_scale_profile}}, \code{\link{terradish_rescale_nu}}
#' @template small-wishart-fit
#' @examples
#' coef(fit)       # Relative log conductance per unit of x.
#' vcov(fit)       # Joint uncertainty in conductance parameters.
#' confint(fit)    # A model-based interval conditional on the supplied nu.
#' nobs(fit)       # Number of focal sites, the sample size for AICc and BIC.
#' @importFrom stats nobs vcov confint
NULL

.no_structure_boundary <- function(fit) {
  if (is.null(fit$no_structure_boundary)) isTRUE(fit$boundary)
  else isTRUE(fit$no_structure_boundary)
}

.graph_fingerprint <- function(graph) {
  path <- tempfile("terradish-graph-")
  on.exit(unlink(path))
  saveRDS(list(adj = graph$adj, demes = graph$demes,
               coordinates = graph$vertex_coordinates, x = graph$x), path,
          compress = FALSE, version = 2)
  unname(tools::md5sum(path))
}

.measurement_columns <- function(model) {
  pairwise <- attr(model, "pairwise_covariates", exact = TRUE)
  if (!is.null(pairwise)) return(as.matrix(pairwise))
  kernels <- attr(model, "kernel_covariates", exact = TRUE)
  if (is.null(kernels)) return(NULL)
  distances <- attr(kernels, "distances", exact = TRUE)
  out <- vapply(seq_len(dim(distances)[3]), function(k)
    distances[, , k][lower.tri(distances[, , k])], numeric(choose(dim(distances)[1], 2)))
  colnames(out) <- dimnames(distances)[[3]]
  out
}

.chibar_probability <- function(statistic, k) {
  if (statistic <= 0) return(1)
  sum(stats::dbinom(seq_len(k), k, .5) *
        stats::pchisq(statistic, seq_len(k), lower.tail = FALSE))
}

#' @rdname terradish_inference
#' @export
vcov.terradish <- function(object, ...) {
  theta <- coef(object)
  if (!length(theta)) return(matrix(numeric(), 0, 0))
  H <- object$mle$hessian
  if (is.null(H)) stop("The fit does not retain its conductance Hessian.", call. = FALSE)
  V <- solve(-H)
  dimnames(V) <- list(names(theta), names(theta))
  V
}

#' @rdname terradish_inference
#' @export
nobs.terradish <- function(object, ...) {
  # focal sites (individuals or populations), not site pairs
  object$dim[["focal"]]
}

#' @rdname terradish_inference
#' @export
confint.terradish <- function(object, parm, level = 0.95, ...) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`level` must be between zero and one.", call. = FALSE)
  theta <- coef(object)
  if (missing(parm)) parm <- names(theta)
  if (is.numeric(parm)) parm <- names(theta)[parm]
  if (anyNA(parm) || !all(parm %in% names(theta))) stop("Unknown parameter in `parm`.")
  se <- sqrt(diag(vcov.terradish(object)))
  z <- stats::qnorm((1 + level) / 2)
  out <- cbind(theta - z * se, theta + z * se)
  colnames(out) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), " %")
  info <- object$gaussian_scale_info
  if (is.null(info) && !is.null(object$submodels$f_internal))
    info <- attr(object$submodels$f_internal, "gaussian_scale_info")
  if (!is.null(info)) {
    sigma_names <- paste0("sigma.", info$scale_vars)
    out[sigma_names, 1] <- pmax(out[sigma_names, 1], info$lower[info$scale_vars])
    out[sigma_names, 2] <- pmin(out[sigma_names, 2], info$upper[info$scale_vars])
    attr(out, "note") <- "Gaussian scale Wald intervals are truncated at the fitted bounds."
    for (layer in names(object$scale_profiles)) {
      profile <- object$scale_profiles[[layer]]
      if (isTRUE(all.equal(profile$level, level))) {
        out[paste0("sigma.", layer), ] <- profile$interval
        attr(out, "note") <- "Stored Gaussian profile intervals are used where available; remaining scale intervals are truncated Wald intervals."
      }
    }
  }
  note <- attr(out, "note")
  out <- out[parm, , drop = FALSE]
  if (!is.null(note)) attr(out, "note") <- note
  out
}

#' Rescale the effective information in a Wishart fit
#'
#' Changes \code{nu} without refitting the coefficients. Wishart objectives
#' and their derivatives are proportional to this user-supplied quantity.
#' @param fit A full or slim fitted Wishart model.
#' @param nu New finite, positive value of \code{nu}.
#' @return A fitted object with rescaled likelihood, curvature, AIC, and
#' comparison metadata. Estimates are unchanged; \code{nu_original} records
#' the original value across repeated rescalings. Standard errors computed by
#' \code{summary()} and \code{vcov()} use the rescaled curvature.
#' @details This is a sensitivity calculation, not an estimate of effective
#' information. It does not remedy misspecification or dependent histories.
#' If \code{nu} is divided by 20, standard errors grow by \code{sqrt(20)};
#' the likelihood and likelihood-ratio statistics are divided by 20. AIC
#' changes because its likelihood term changes while its parameter penalty
#' stays fixed. Report conclusions over a scientifically plausible range.
#' @seealso \code{\link{terradish_inference}}, \code{\link{aic_table}},
#'   \code{\link{terradish_cv_folds}}
#' @template small-wishart-fit
#' @examples
#' sensitivity <- terradish_rescale_nu(fit, nu = 50)
#' rbind(original = coef(fit), sensitivity = coef(sensitivity))
#' rbind(original = sqrt(diag(vcov(fit))),
#'       sensitivity = sqrt(diag(vcov(sensitivity))))
#' @export
terradish_rescale_nu <- function(fit, nu) {
  family <- .terradish_fit_comparison_contract(fit)$likelihood_family
  if (length(family) != 1L || !family %in% c("wishart_covariance", "generalized_wishart_distance"))
    stop("`fit` must use a Wishart likelihood.", call. = FALSE)
  old <- fit$comparison$nu
  if (length(nu) != 1L || !is.finite(nu) || nu <= 0 ||
      length(old) != 1L || !is.finite(old) || old <= 0)
    stop("Both original and new `nu` must be finite and positive.", call. = FALSE)
  ratio <- nu / old
  scaled_names <- c("objective", "loglik", "gradient", "gradient_internal",
    "hessian", "hessian_internal", "phi_hessian", "partial_X", "partial_S",
    "num_gradient", "num_hessian", "num_partial_X", "num_partial_S")
  for (part in c("mle", "fit"))
    for (name in intersect(names(fit[[part]]), scaled_names))
      if (is.numeric(fit[[part]][[name]])) fit[[part]][[name]] <- fit[[part]][[name]] * ratio
  if (!is.null(fit$fit$phi_vcov_joint)) fit$fit$phi_vcov_joint <- fit$fit$phi_vcov_joint / ratio
  fit$loglik <- fit$loglik * ratio
  fit$aic <- -2 * fit$loglik + 2 * fit$df
  fit$comparison$nu <- nu
  fit$call$nu <- nu
  if (is.null(fit$nu_original)) fit$nu_original <- old
  if (!is.null(fit$convergence$max_abs_projected_gradient))
    fit$convergence$max_abs_projected_gradient <- fit$convergence$max_abs_projected_gradient * ratio
  # A profile interval belongs to the original objective scale; its grid can
  # be reused later, but retaining the old interval would be misleading.
  fit$scale_profiles <- NULL
  fit
}

#' Environmental effects in resistance-distance units
#'
#' Reports the resistance-distance equivalent of one unit of environmental
#' difference, using the environmental coefficient divided by the resistance
#' coefficient in the additive measurement model.
#' @param fit A fitted model from \code{mlpe_covariates()} or
#'   \code{wishart_covariates()}, or a fit with no environmental terms.
#' @return A data frame with covariate, ratio, SE, lower, upper, and family.
#'   Intervals use a 95 percent delta-method approximation based on the joint
#'   covariance of both coefficients, including fitted-conductance uncertainty.
#'   Fits without environmental terms return a zero-row data frame.
#' @details The ratios are \eqn{\gamma_k/\beta} for MLPE and
#'   \eqn{\lambda_k/\tau} for Wishart models. They require the same environmental
#'   transform and scaling for comparison. Their interpretation agrees exactly
#'   for the same covariance-derived distance and approximately for distances
#'   that scale linearly with it. Separate fits estimate their own conductance
#'   surfaces, so their resistance distances can differ. Absolute measurement
#'   coefficients are not comparable across families.
#'
#'   A zero resistance coefficient gives undefined ratios and uncertainty,
#'   returned as NA. A zero Wishart kernel coefficient receives a one-sided
#'   upper limit and no symmetric lower limit. Boundary inference remains
#'   approximate. The distance-equivalent ratio interpretation also appears in
#'   BEDASSLE, although that model uses geographic rather than resistance distance.
#' @references Bradburd GS, Ralph PL, Coop GM (2013). Disentangling the effects
#'   of geographic and ecological isolation on genetic differentiation.
#'   Evolution, 67, 3258--3273. \doi{10.1111/evo.12193}.
#' @seealso \code{\link{pairwise_covariates}}, \code{\link{mlpe_covariates}},
#'   \code{\link{wishart_covariates}}, \code{\link{terradish_rescale_nu}}
#' @template small-wishart-fit
#' @examples
#' environment <- pairwise_covariates(
#'   pairwise_endpoint_covariates(data.frame(climate = c(2, 4, 1, 7, 5, 9))))
#' measurement <- wishart_covariates(environment, model = "wishart_covariance")
#' kernel <- attr(measurement, "kernel_covariates")[, , 1]
#' response <- response + 0.1 * kernel
#' joint <- terradish(response ~ x, graph, measurement_model = measurement,
#'                    nu = 1000)
#' terradish_ibe_ratio(joint)
#' # Each ratio is the resistance-distance equivalent of one unit of difference.
#' @export
terradish_ibe_ratio <- function(fit) {
  phi <- fit$fit$phi[, 1]
  family <- .terradish_fit_comparison_contract(fit)$likelihood_family
  lambda <- grep("^lambda_", names(phi), value = TRUE)
  gamma <- if (!length(lambda)) intersect(colnames(fit$comparison$measurement_columns), names(phi)) else character()
  environmental <- c(lambda, gamma)
  out <- data.frame(covariate = environmental, ratio = rep(NA_real_, length(environmental)),
    SE = rep(NA_real_, length(environmental)), lower = rep(NA_real_, length(environmental)),
    upper = rep(NA_real_, length(environmental)),
    family = rep(if (is.null(family)) NA_character_ else family, length(environmental)))
  if (!length(environmental)) return(out)
  denominator <- if (length(lambda)) "tau" else "beta"
  if (!is.finite(phi[denominator]) || phi[denominator] == 0) {
    attr(out, "note") <- "The resistance coefficient is zero; environmental ratios are undefined."
    return(out)
  }
  V <- fit$fit$phi_vcov_joint
  if (is.null(V)) stop("This fit lacks joint nuisance covariance; refit with the current package.", call. = FALSE)
  for (i in seq_along(environmental)) {
    nm <- environmental[i]
    out$ratio[i] <- phi[nm] / phi[denominator]
    derivative <- setNames(rep(0, length(phi)), names(phi))
    derivative[nm] <- 1 / phi[denominator]
    derivative[denominator] <- -phi[nm] / phi[denominator]^2
    variance <- as.numeric(crossprod(derivative, V %*% derivative))
    out$SE[i] <- if (variance >= 0) sqrt(variance) else NA_real_
    out$lower[i] <- out$ratio[i] - qnorm(.975) * out$SE[i]
    out$upper[i] <- out$ratio[i] + qnorm(.975) * out$SE[i]
    if (nm %in% lambda && phi[nm] == 0) {
      out$lower[i] <- NA_real_
      out$upper[i] <- qnorm(.95) * out$SE[i]
    }
  }
  out
}
