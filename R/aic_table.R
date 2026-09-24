#' Rank fitted terradish models by information criterion
#'
#' Creates a model-selection table from fitted \code{terradish} models using
#' AIC, AICc, or BIC.
#'
#' @param mod_list List of fitted \code{terradish} models.
#' @param AICc Should second-order AIC (Akaike's information criterion with
#'   small-sample correction) be used instead of AIC?  When \code{TRUE}, the
#'   correction uses \eqn{n} = number of focal sampling sites (not the number
#'   of pairwise observations).  A common guideline is to prefer AICc over AIC
#'   when \eqn{n / K < 40}, where \eqn{K} is the number of estimated parameters.
#' @param BIC Should BIC be used instead of AIC? BIC uses
#'   \eqn{n = n_{sites}(n_{sites}-1)/2} as an implementation convention. Those
#'   pairs are not independent in MLPE and related designs, so report the
#'   convention and do not treat it as a uniquely determined effective sample
#'   size.
#' @param mod_names Optional model names. By default the right-hand side of
#'   each fitted formula is used. For MLPE measurement models with additional
#'   pairwise covariates, the default appends \code{[mlpe:n]}, where \code{n}
#'   is the number of added pairwise covariate columns.
#' @param verbose Should the table be printed to the console?
#'
#' @details
#' At nominal Wishart information, information criteria can favor extensions
#' that are unsupported across replicate histories. Use
#' \code{\link{terradish_cv_folds}} to select conductance terms and
#' \code{\link{terradish_rescale_nu}} to assess inferential sensitivity.
#'
#' Information-criterion comparison requires the same observed response,
#' focal sites, graph domain, and likelihood family. Wishart fits
#' must also use the same effective degrees of freedom, \code{nu}.
#' \code{leastsquares} and \code{mlpe} are Gaussian likelihoods for the same
#' pairwise-distance response and can be ranked when those conditions hold.
#' \code{generalized_wishart} and \code{wishart_covariance} use different
#' response representations and must not be ranked together. Cross-family
#' rankings, such as \code{mlpe} versus \code{generalized_wishart}, are invalid.
#'
#' The function checks response values, model family, fitted dimensions,
#' recorded \code{nu}, and retained graph signatures. Older fits without
#' signatures require the user to verify the graph domain.
#'
#' @return A data frame containing model ranks, parameter counts, information
#'   criterion values, delta values, weights, cumulative weights, and
#'   log-likelihoods. The \code{converged} column records convergence code 0;
#'   \code{boundary} marks active parameter bounds or no resistance structure.
#'   A warning identifies tables containing nonconverged or boundary fits.
#'
#' @examples
#' \donttest{
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' keep <- 1:10
#' melip.Fst_small <- melip.Fst[keep, keep]
#' # Coarsen the rasters so this example remains quick on CRAN.
#' covariates <- terra::aggregate(c(melip.altitude, melip.forestcover),
#'                                fact = 3, na.rm = TRUE)
#' covariates <- scale_covariates(covariates)
#' names(covariates) <- c("altitude", "forestcover")
#' surface_small <- conductance_surface(covariates, melip.coords[keep], directions = 8)
#'
#' fit1 <- terradish(melip.Fst_small ~ altitude, surface_small,
#'                   loglinear_conductance, leastsquares,
#'                   control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' fit2 <- terradish(melip.Fst_small ~ altitude + forestcover, surface_small,
#'                   loglinear_conductance, leastsquares,
#'                   control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#'
#' aic_table(list(fit1, fit2))
#' aic_table(list(fit1, fit2), AICc = TRUE)
#'
#' }
#' @export
aic_table <- function(mod_list, AICc = FALSE, BIC = FALSE, mod_names = NULL, verbose = FALSE)
{
  .terradish_assert_comparable_fits(mod_list, purpose = "information criterion")

  if (is.null(mod_names))
    mod_names <- vapply(mod_list, .default_model_name, character(1))

  mod_loglik <- vapply(mod_list, function(x) x$loglik, numeric(1))
  mod_AIC <- vapply(mod_list, function(x) x$aic, numeric(1))
  mod_df <- vapply(mod_list, function(x) x$df, numeric(1))

  if (!isTRUE(AICc) && !isTRUE(BIC))
  {
    delta <- mod_AIC - min(mod_AIC)
    wt <- exp(-0.5 * delta)
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      Delta_AIC = delta,
                      AIC_wt = wt / sum(wt),
                      .fit_index = seq_along(mod_list),
                      row.names = NULL)
    tab <- tab[order(tab$AIC), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$AIC_wt)
    tab$loglik <- mod_loglik[tab$.fit_index]

    tab[, 3:7] <- round(tab[, 3:7], digits = 4)
  }
  else if (isTRUE(AICc))
  {
    if (isTRUE(BIC))
      stop("Set only one of `AICc` or `BIC` to TRUE")

    mod_n <- vapply(mod_list, function(x) x$dim[["focal"]], numeric(1))
    if (any(mod_n <= mod_df + 1))
      stop("AICc is undefined when the number of focal sites is not greater than K + 1.",
           call. = FALSE)
    mod_AICc <- -2 * mod_loglik + 2 * mod_df * (mod_n / (mod_n - mod_df - 1))
    delta <- mod_AICc - min(mod_AICc)
    wt <- exp(-0.5 * delta)
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      AICc = mod_AICc,
                      Delta_AICc = delta,
                      AICc_wt = wt / sum(wt),
                      .fit_index = seq_along(mod_list),
                      row.names = NULL)
    tab <- tab[order(tab$AICc), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$AICc_wt)
    tab$loglik <- mod_loglik[tab$.fit_index]

    tab[, 3:8] <- round(tab[, 3:8], digits = 4)
  }
  else
  {
    mod_n <- vapply(mod_list, function(x) x$dim[["focal"]], numeric(1))
    mod_pairs <- mod_n * (mod_n - 1) / 2
    mod_BIC <- -2 * mod_loglik + mod_df * log(mod_pairs)
    delta <- mod_BIC - min(mod_BIC)
    wt <- exp(-0.5 * delta)
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      BIC = mod_BIC,
                      Delta_BIC = delta,
                      BIC_wt = wt / sum(wt),
                      .fit_index = seq_along(mod_list),
                      row.names = NULL)
    tab <- tab[order(tab$BIC), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$BIC_wt)
    tab$loglik <- mod_loglik[tab$.fit_index]

    tab[, 3:8] <- round(tab[, 3:8], digits = 4)
  }

  tab$converged <- vapply(mod_list, function(fit)
    isTRUE(fit$convergence$code == 0L), logical(1))[tab$.fit_index]
  tab$boundary <- vapply(mod_list, function(fit)
    isTRUE(fit$convergence$boundary) || isTRUE(fit$fit$boundary), logical(1))[tab$.fit_index]
  tab$.fit_index <- NULL
  if (any(!tab$converged | tab$boundary))
    warning("Some models have not converged or have parameters on a boundary; inspect fit diagnostics before ranking.",
            call. = FALSE)
  if (isTRUE(verbose))
    print(tab, row.names = FALSE)
  else
    tab
}

.model_rhs_label <- function(model_formula)
{
  if (!inherits(model_formula, "formula"))
    return(NULL)

  if (length(model_formula) >= 3L)
    return(paste(deparse(model_formula[[3]]), collapse = ""))
  if (length(model_formula) == 2L)
    return(paste(deparse(model_formula[[2]]), collapse = ""))

  NULL
}

.model_formula_label <- function(mod)
{
  rhs <- .model_rhs_label(mod$formula)
  if (!is.null(rhs))
    return(rhs)

  call_formula <- mod$call$formula
  rhs <- .model_rhs_label(call_formula)
  if (!is.null(rhs))
    return(rhs)

  if (is.call(call_formula) && length(call_formula) >= 3L)
    return(paste(deparse(call_formula[[3]]), collapse = ""))
  if (is.call(call_formula) && length(call_formula) == 2L)
    return(paste(deparse(call_formula[[2]]), collapse = ""))

  "<unknown>"
}

.mlpe_covariate_count <- function(mod)
{
  measurement_model <- mod$submodels$g
  if (!is.function(measurement_model))
    return(0L)

  pairwise_covariates <- attr(measurement_model, "pairwise_covariates",
                              exact = TRUE)
  if (is.null(pairwise_covariates))
    return(0L)

  pairwise_mat <- tryCatch(as.matrix(pairwise_covariates),
                           error = function(e) NULL)
  if (is.null(pairwise_mat))
    return(0L)

  n_covariates <- ncol(pairwise_mat)
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      is.na(n_covariates) || n_covariates < 1)
    return(0L)

  as.integer(n_covariates)
}

.default_model_name <- function(mod)
{
  base <- .model_formula_label(mod)
  n_mlpe <- .mlpe_covariate_count(mod)
  if (n_mlpe > 0L)
    paste0(base, " [mlpe:", n_mlpe, "]")
  else
    base
}

.terradish_measurement_model_name <- function(model)
{
  known <- c("leastsquares", "mlpe", "generalized_wishart", "wishart_covariance")
  ns_env <- asNamespace("terradish")
  for (nm in known)
  {
    obj <- get0(nm, envir = ns_env, mode = "function", inherits = FALSE)
    if (!is.null(obj) &&
        identical(formals(model), formals(obj)) &&
        identical(body(model), body(obj)))
      return(nm)
  }

  base_model <- attr(model, "base_model", exact = TRUE)
  if (is.character(base_model) && length(base_model) == 1L)
    return(base_model)

  NULL
}

.terradish_measurement_family <- function(model_name)
{
  if (!is.character(model_name) || length(model_name) != 1L)
    return(NULL)

  switch(model_name,
         leastsquares = "gaussian_distance",
         mlpe = "gaussian_distance",
         generalized_wishart = "generalized_wishart_distance",
         wishart_covariance = "wishart_covariance",
         NULL)
}

.terradish_comparison_contract <- function(measurement_model, nu = NULL,
                                           response = NULL)
{
  model_name <- .terradish_measurement_model_name(measurement_model)
  family <- .terradish_measurement_family(model_name)

  list(
    measurement_model = model_name,
    likelihood_family = family,
    nu = if (isTRUE(grepl("wishart", family, fixed = TRUE))) nu else NULL,
    response = response
  )
}

.terradish_fit_comparison_contract <- function(fit)
{
  contract <- fit$comparison
  if (!is.list(contract))
    contract <- list()

  measurement_model <- fit$submodels$g
  if (is.null(contract$measurement_model) && is.function(measurement_model))
    contract$measurement_model <- .terradish_measurement_model_name(measurement_model)
  if (is.null(contract$likelihood_family))
    contract$likelihood_family <- .terradish_measurement_family(contract$measurement_model)
  if (is.null(contract$response) && !is.null(fit$fit$response))
    contract$response <- fit$fit$response

  contract
}

.terradish_assert_comparable_fits <- function(fits, purpose = c("information criterion",
                                                                 "likelihood-ratio test",
                                                                 "cross-validation comparison"))
{
  purpose <- match.arg(purpose)
  if (length(fits) < 2L)
    stop("At least two fitted models are required for comparison.", call. = FALSE)

  dim_keys <- vapply(fits, function(x) paste(x$dim, collapse = "|"), character(1))
  if (length(unique(dim_keys)) != 1L)
    stop("Models must use the same focal sites and graph dimensions for a valid ",
         purpose, ".", call. = FALSE)

  contracts <- lapply(fits, .terradish_fit_comparison_contract)
  families <- vapply(contracts,
                     function(x) if (is.null(x$likelihood_family)) NA_character_ else x$likelihood_family,
                     character(1))
  known_families <- unique(stats::na.omit(families))
  if (anyNA(families))
    stop("Could not identify every model's likelihood family, so a valid ",
         purpose, " cannot be verified. Use a built-in measurement model or set ",
         "its `base_model` attribute to the corresponding built-in model.",
         call. = FALSE)
  if (length(known_families) > 1L)
    stop("Models use different likelihood families. Do not compare Gaussian distance, ",
         "generalized-Wishart distance, and covariance-Wishart fits by ", purpose, ".",
         call. = FALSE)

  responses <- lapply(contracts, `[[`, "response")
  have_responses <- !vapply(responses, is.null, logical(1))
  if (all(have_responses) &&
      !all(vapply(responses[-1L], identical, logical(1), responses[[1L]])))
    stop("Models must use the same response matrix for a valid ", purpose, ".",
         call. = FALSE)

  is_wishart <- grepl("wishart", families, fixed = TRUE)
  if (any(is_wishart, na.rm = TRUE))
  {
    nu <- vapply(contracts, function(x) {
      if (is.null(x$nu) || length(x$nu) != 1L) NA_real_ else as.numeric(x$nu)
    }, numeric(1))
    known_nu <- unique(nu[is.finite(nu)])
    if (length(known_nu) > 1L)
      stop("Wishart models must use the same effective degrees of freedom (`nu`) for a valid ",
           purpose, ".", call. = FALSE)
  }

  invisible(contracts)
}
