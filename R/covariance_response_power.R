#' Assess covariance-response parameter recovery by simulation
#'
#' Runs a fast simulation study for covariance-response \code{terradish} fits.
#' Candidate sampling designs are created by subsetting the focal sites in an
#' existing \code{\link{conductance_surface}} object, responses are generated
#' with \code{\link{simulate_covariance_response}}, and one or more candidate
#' conductance models are refit to each simulated covariance matrix.
#'
#' @param theta True conductance parameters. By default these are interpreted on
#'   the same external scale returned by \code{\link{coef}} for fitted
#'   \code{terradish} models. For most conductance models the internal and
#'   external scales are identical. For
#'   \code{\link{gaussian_smoothed_loglinear_conductance}}, the external scale
#'   reports \code{sigma.*} parameters in map units.
#' @param formula Formula for the true conductance model. The left hand side,
#'   if supplied, is ignored.
#' @param data A \code{\link{conductance_surface}} object containing all
#'   candidate focal sites.
#' @param sample_sizes Integer vector of focal sample sizes to assess.
#' @param strategies Sampling strategies. Supported values are
#'   \code{"spacefill"}, \code{"random"}, and \code{"sequential"}.
#' @param conductance_model Conductance-model factory used to simulate the true
#'   response.
#' @param fit_models Optional named list of candidate models to refit. Each
#'   element may contain \code{formula}, \code{conductance_model}, \code{theta},
#'   \code{optimizer}, \code{control}, \code{solver}, \code{solver_control},
#'   \code{approximation}, \code{approximation_control}, \code{nonnegative},
#'   and \code{leverage}. When omitted, the true model is refit.
#' @param tau Nonnegative scaling applied to the conductance-implied covariance
#'   matrix during simulation.
#' @param sigma Nonnegative nugget variance added to the covariance diagonal
#'   during simulation.
#' @param nu Effective Wishart degrees of freedom, passed to both simulation
#'   and \code{\link{wishart_covariance}} fitting.  For biallelic SNPs this is
#'   usually the number of retained SNPs; for microsatellites use approximately
#'   \eqn{\sum_l (K_l - 1)} where \eqn{K_l} is the number of observed alleles at
#'   locus \eqn{l}.  Because \code{nu} acts as an effective sample size, it is a
#'   primary lever for power: hold the sampling design fixed and vary \code{nu}
#'   to see how many markers are needed to recover an effect.  See
#'   \code{\link{wishart_covariance}} for details on how \code{nu} scales
#'   inference.
#' @param nsim Number of covariance-response simulations per sampling design.
#' @param n_designs Number of independent site designs per sample size for
#'   \code{strategy = "random"}. Deterministic strategies are evaluated once.
#' @param seed Optional random seed for reproducible site and response
#'   simulation.
#' @param theta_scale Is \code{theta} supplied on the \code{"external"} or
#'   \code{"internal"} conductance-model scale?
#' @param effect_parameters Optional character vector of parameters used in the
#'   all-effects recovery summary. Defaults to the nonzero true parameters.
#' @param alpha Significance level for Wald-style parameter recovery summaries.
#' @param conductance_cor_threshold Threshold used for the conductance-surface
#'   recovery summary.
#' @param optimizer,control,solver,solver_control,approximation,approximation_control,nonnegative,leverage,cores
#'   Defaults passed to \code{\link{terradish}} for candidate models unless a
#'   model-specific value is supplied in \code{fit_models}.
#' @param verbose Logical; print progress messages?
#'
#' @details
#' This helper avoids forward-time genetic simulation. It assumes a proposed
#' conductance surface and draws covariance responses from the same Wishart
#' likelihood fitted by \code{\link{wishart_covariance}}. That makes it useful
#' for approximate power and sample-size screening before investing in slower
#' individual-based simulations.
#'
#' Sampling designs are efficient because the raster graph is not rebuilt for
#' each sample size. The helper subsets \code{data$demes} and the corresponding
#' right-hand-side columns, so all candidate designs share the same graph,
#' covariate matrix, and conductance-model basis expansion.
#'
#' Fits that land on the no-signal boundary are retained as non-detections in
#' parameter-level power summaries. Their coefficient estimates and standard
#' errors are reported as \code{NA}, matching explicit fitting failures.
#'
#' For spline conductance models, individual basis coefficients are not usually
#' the most interpretable recovery target. Inspect \code{parameter_summary} for
#' basis-level behavior, but prefer \code{summary$conductance_power},
#' \code{summary$mean_conductance_cor}, and model-selection rates when judging
#' whether the non-linear conductance surface was recovered.
#'
#' @return An object of class \code{"terradish_covariance_power"} containing:
#' \describe{
#'   \item{\code{results}}{One row per fitted model and simulated response.}
#'   \item{\code{summary}}{Scenario-level fit, model-selection, and
#'     conductance-recovery summaries.}
#'   \item{\code{parameter_results}}{One row per assessed parameter and
#'     simulation replicate.}
#'   \item{\code{parameter_summary}}{Parameter-level power, bias, RMSE, and
#'     coverage summaries.}
#'   \item{\code{settings}}{Simulation settings and true parameters.}
#' }
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
#' covariates <- scale_covariates(covariates)
#' surface <- conductance_surface(covariates, melip.coords,
#'                                directions = 8, saveStack = TRUE)
#'
#' spline_model <- smooth_loglinear_conductance(
#'   ~ forestcover + s(altitude, df = 4), surface$x)
#' theta <- attr(spline_model, "default")
#' theta["forestcover"] <- -0.15
#' theta[grepl("^s\\(altitude\\)", names(theta))] <-
#'   c(-0.65, 1.1, -0.9, 0.45)
#'
#' gaussian_model <- gaussian_smoothed_loglinear_conductance(
#'   surface, scale_vars = "altitude")
#'
#' power <- covariance_response_power(
#'   theta = theta,
#'   formula = ~ forestcover + s(altitude, df = 4),
#'   data = surface,
#'   sample_sizes = c(8, 12, 16),
#'   strategies = c("spacefill", "random"),
#'   conductance_model = smooth_loglinear_conductance,
#'   fit_models = list(
#'     spline = list(formula = ~ forestcover + s(altitude, df = 4),
#'                   conductance_model = smooth_loglinear_conductance),
#'     gaussian = list(formula = ~ altitude + forestcover,
#'                     conductance_model = gaussian_model)
#'   ),
#'   tau = 0.8,
#'   sigma = 0.12,
#'   nu = 150,
#'   nsim = 50,
#'   n_designs = 5,
#'   seed = 1,
#'   control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
#' )
#' power$summary
#' power$parameter_summary
#' }
#'
#' @export
covariance_response_power <- function(theta,
                                      formula,
                                      data,
                                      sample_sizes,
                                      strategies = c("spacefill", "random"),
                                      conductance_model = loglinear_conductance,
                                      fit_models = NULL,
                                      tau = 1,
                                      sigma = 0,
                                      nu,
                                      nsim = 20L,
                                      n_designs = 1L,
                                      seed = NULL,
                                      theta_scale = c("external", "internal"),
                                      effect_parameters = NULL,
                                      alpha = 0.05,
                                      conductance_cor_threshold = 0.9,
                                      optimizer = c("auto", "newton", "bfgs"),
                                      control = NewtonRaphsonControl(maxit = 8,
                                                                     verbose = FALSE),
                                      solver = c("auto", "direct", "amg",
                                                 "pcg", "pcg_jacobi"),
                                      solver_control = NULL,
                                      approximation = c("none", "landmark",
                                                        "coarse_raster"),
                                      approximation_control = NULL,
                                      nonnegative = TRUE,
                                      leverage = FALSE,
                                      cores = 1L,
                                      verbose = FALSE)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(inherits(conductance_model,
                     c("terradish_conductance_model_factory",
                       "radish_conductance_model_factory")))
  stopifnot(length(cores) == 1L, is.numeric(cores), cores >= 1)
  stopifnot(is.numeric(nsim), length(nsim) == 1L, nsim >= 1)
  stopifnot(is.numeric(n_designs), length(n_designs) == 1L, n_designs >= 1)
  stopifnot(is.numeric(alpha), length(alpha) == 1L,
            is.finite(alpha), alpha > 0, alpha < 1)
  stopifnot(is.numeric(conductance_cor_threshold),
            length(conductance_cor_threshold) == 1L,
            is.finite(conductance_cor_threshold),
            conductance_cor_threshold >= -1,
            conductance_cor_threshold <= 1)
  theta_scale <- match.arg(theta_scale)
  optimizer <- match.arg(optimizer)
  solver <- match.arg(solver)
  approximation <- match.arg(approximation)
  strategies <- match.arg(strategies,
                          c("spacefill", "random", "sequential"),
                          several.ok = TRUE)

  n_focal <- length(data$demes)
  sample_sizes <- unique(as.integer(sample_sizes))
  if (!length(sample_sizes))
    stop("`sample_sizes` must contain at least one value.", call. = FALSE)
  if (any(is.na(sample_sizes)) || any(sample_sizes < 2L))
    stop("`sample_sizes` must contain integers of at least 2.", call. = FALSE)
  if (any(sample_sizes > n_focal))
    stop("`sample_sizes` cannot exceed the number of focal sites in `data`.",
         call. = FALSE)

  true_model <- conductance_model(.covariance_power_rhs_formula(formula),
                                  data$x)
  theta_info <- .covariance_power_theta(
    theta = theta,
    conductance_model = true_model,
    theta_scale = theta_scale
  )
  theta_external <- theta_info$external
  theta_internal <- theta_info$internal
  true_conductance <- true_model(theta_internal)$conductance

  if (is.null(effect_parameters))
    effect_parameters <- names(theta_external)[is.finite(theta_external) &
                                                 theta_external != 0]
  effect_parameters <- unique(effect_parameters)

  model_specs <- .covariance_power_model_specs(
    fit_models = fit_models,
    formula = formula,
    conductance_model = conductance_model,
    data = data
  )

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (old_seed_exists)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (old_seed_exists)
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed))
    set.seed(seed)

  results <- list()
  parameter_results <- list()
  result_i <- 0L
  parameter_i <- 0L
  scenario_id <- 0L

  for (sample_size in sample_sizes)
  {
    for (strategy in strategies)
    {
      design_count <- if (identical(strategy, "random"))
        as.integer(n_designs)
      else
        1L

      for (design in seq_len(design_count))
      {
        scenario_id <- scenario_id + 1L
        design_seed <- .covariance_power_next_seed()
        idx <- .covariance_power_indices(
          data = data,
          n_focal = n_focal,
          sample_size = sample_size,
          strategy = strategy,
          seed = design_seed
        )
        subset_data <- .covariance_power_subset_graph(data, idx)

        for (sim_id in seq_len(as.integer(nsim)))
        {
          sim_seed <- .covariance_power_next_seed()
          if (isTRUE(verbose))
            message("sample_size=", sample_size,
                    " strategy=", strategy,
                    " design=", design,
                    " sim=", sim_id)

          sim <- simulate_covariance_response(
            theta = theta_internal,
            formula = formula,
            data = subset_data,
            conductance_model = conductance_model,
            tau = tau,
            sigma = sigma,
            nu = nu,
            nsim = 1L,
            seed = sim_seed,
            cores = cores
          )
          genetic_cov <- sim$covariance

          fit_rows <- list()
          fit_objects <- vector("list", length(model_specs))
          names(fit_objects) <- names(model_specs)

          for (model_name in names(model_specs))
          {
            spec <- model_specs[[model_name]]
            fit_formula <- .covariance_power_response_formula(spec$formula)
            fit_env <- new.env(parent = parent.frame())
            fit_env$genetic_cov <- genetic_cov
            environment(fit_formula) <- fit_env

            fit_args <- list(
              formula = fit_formula,
              data = subset_data,
              conductance_model = spec$conductance_model,
              measurement_model = wishart_covariance,
              nu = nu,
              theta = .covariance_power_spec_value(spec, "theta", NULL),
              leverage = .covariance_power_spec_value(spec, "leverage",
                                                       leverage),
              nonnegative = .covariance_power_spec_value(spec, "nonnegative",
                                                         nonnegative),
              optimizer = .covariance_power_spec_value(spec, "optimizer",
                                                       optimizer),
              control = .covariance_power_spec_value(spec, "control",
                                                     control),
              cores = cores,
              solver = .covariance_power_spec_value(spec, "solver", solver),
              solver_control = .covariance_power_spec_value(
                spec, "solver_control", solver_control),
              approximation = .covariance_power_spec_value(
                spec, "approximation", approximation),
              approximation_control = .covariance_power_spec_value(
                spec, "approximation_control", approximation_control)
            )

            elapsed_start <- proc.time()[["elapsed"]]
            fit <- tryCatch(
              suppressWarnings(do.call(terradish, fit_args)),
              error = function(e) e
            )
            elapsed <- proc.time()[["elapsed"]] - elapsed_start

            if (inherits(fit, "error"))
            {
              fit_status <- paste("ERROR:", conditionMessage(fit))
              fit_row <- .covariance_power_result_row(
                scenario_id = scenario_id,
                strategy = strategy,
                sample_size = sample_size,
                design = design,
                simulation = sim_id,
                model = model_name,
                status = fit_status,
                elapsed = elapsed,
                n_focal = sample_size,
                n_pairs = sample_size * (sample_size - 1L) / 2,
                fit = NULL,
                true_conductance = true_conductance,
                conductance_cor_threshold = conductance_cor_threshold
              )
              fit_rows[[model_name]] <- fit_row
              common_parameters <- intersect(effect_parameters,
                                             spec$parameter_names)
              if (length(common_parameters))
              {
                for (parameter in common_parameters)
                {
                  parameter_i <- parameter_i + 1L
                  parameter_results[[parameter_i]] <-
                    .covariance_power_parameter_row(
                      fit_row = fit_row,
                      parameter = parameter,
                      true = theta_external[[parameter]],
                      estimate = NA_real_,
                      se = NA_real_,
                      alpha = alpha
                    )
                }
              }
              next
            }

            fit_objects[[model_name]] <- fit
            fit_row <- .covariance_power_result_row(
              scenario_id = scenario_id,
              strategy = strategy,
              sample_size = sample_size,
              design = design,
              simulation = sim_id,
              model = model_name,
              status = "OK",
              elapsed = elapsed,
              n_focal = sample_size,
              n_pairs = sample_size * (sample_size - 1L) / 2,
              fit = fit,
              true_conductance = true_conductance,
              conductance_cor_threshold = conductance_cor_threshold
            )
            fit_rows[[model_name]] <- fit_row

            if (!fit_row$fit_ok)
            {
              common_parameters <- intersect(effect_parameters,
                                             spec$parameter_names)
              if (length(common_parameters))
              {
                for (parameter in common_parameters)
                {
                  parameter_i <- parameter_i + 1L
                  parameter_results[[parameter_i]] <-
                    .covariance_power_parameter_row(
                      fit_row = fit_row,
                      parameter = parameter,
                      true = theta_external[[parameter]],
                      estimate = NA_real_,
                      se = NA_real_,
                      alpha = alpha
                    )
                }
              }
              next
            }

            coef_fit <- coef(fit)
            se_fit <- .covariance_power_se(fit)
            common_parameters <- intersect(effect_parameters, names(coef_fit))
            if (length(common_parameters))
            {
              for (parameter in common_parameters)
              {
                parameter_i <- parameter_i + 1L
                parameter_results[[parameter_i]] <-
                  .covariance_power_parameter_row(
                    fit_row = fit_row,
                    parameter = parameter,
                    true = theta_external[[parameter]],
                    estimate = coef_fit[[parameter]],
                    se = se_fit[[parameter]],
                    alpha = alpha
                  )
              }
            }
          }

          fit_table <- do.call(rbind, fit_rows)
          fit_table <- .covariance_power_mark_selection(fit_table)
          for (row in seq_len(nrow(fit_table)))
          {
            result_i <- result_i + 1L
            results[[result_i]] <- fit_table[row, , drop = FALSE]
          }
        }
      }
    }
  }

  results <- if (length(results))
    do.call(rbind, results)
  else
    data.frame()
  rownames(results) <- NULL

  parameter_results <- if (length(parameter_results))
    do.call(rbind, parameter_results)
  else
    data.frame()
  rownames(parameter_results) <- NULL

  out <- list(
    call = match.call(),
    results = results,
    summary = .covariance_power_add_parameter_power(
      .covariance_power_summary(results, conductance_cor_threshold),
      parameter_results
    ),
    parameter_results = parameter_results,
    parameter_summary = .covariance_power_parameter_summary(parameter_results),
    settings = list(
      formula = formula,
      sample_sizes = sample_sizes,
      strategies = strategies,
      tau = tau,
      sigma = sigma,
      nu = nu,
      nsim = as.integer(nsim),
      n_designs = as.integer(n_designs),
      theta = theta_external,
      theta_internal = theta_internal,
      theta_scale = theta_scale,
      effect_parameters = effect_parameters,
      alpha = alpha,
      conductance_cor_threshold = conductance_cor_threshold
    )
  )
  class(out) <- "terradish_covariance_power"
  out
}

#' @method print terradish_covariance_power
#' @export
print.terradish_covariance_power <- function(x, ...)
{
  cat("Covariance-response power assessment\n")
  cat("Simulations per design:", x$settings$nsim, "\n")
  cat("Sample sizes:", paste(x$settings$sample_sizes, collapse = ", "), "\n")
  cat("Strategies:", paste(x$settings$strategies, collapse = ", "), "\n\n")
  if (nrow(x$summary))
    print(x$summary, row.names = FALSE)
  invisible(x)
}

.covariance_power_rhs_formula <- function(formula)
{
  terms_obj <- terms(formula)
  labels <- attr(terms_obj, "term.labels")
  if (!length(labels))
    stop("`formula` must include at least one conductance covariate.",
         call. = FALSE)
  reformulate(labels)
}

.covariance_power_response_formula <- function(formula)
{
  terms_obj <- terms(formula)
  labels <- attr(terms_obj, "term.labels")
  if (!length(labels))
    return(stats::as.formula("genetic_cov ~ 1"))
  reformulate(labels, response = "genetic_cov")
}

.covariance_power_theta <- function(theta, conductance_model,
                                    theta_scale = c("external", "internal"))
{
  theta_scale <- match.arg(theta_scale)
  default <- attr(.externalize_conductance_model(conductance_model),
                  "default", exact = TRUE)
  if (is.null(default))
    stop("Conductance model does not define default parameters.",
         call. = FALSE)

  theta <- c(theta)
  if (length(theta) != length(default))
    stop("`theta` must have one value per conductance parameter.",
         call. = FALSE)
  if (is.null(names(theta)))
    names(theta) <- names(default)
  if (anyDuplicated(names(theta)))
    stop("`theta` names must be unique.", call. = FALSE)
  missing_names <- setdiff(names(default), names(theta))
  extra_names <- setdiff(names(theta), names(default))
  if (length(missing_names) || length(extra_names))
    stop("`theta` names must match the conductance-model parameters.",
         call. = FALSE)
  theta <- theta[names(default)]

  if (identical(theta_scale, "external"))
  {
    external <- theta
    internal <- .conductance_model_to_internal(theta, conductance_model)
  }
  else
  {
    internal <- theta
    external <- .conductance_model_to_external(theta, conductance_model)
  }
  list(external = external, internal = internal)
}

.covariance_power_model_specs <- function(fit_models, formula,
                                          conductance_model, data)
{
  if (is.null(fit_models))
    fit_models <- list(truth = list(formula = formula,
                                    conductance_model = conductance_model))
  if (!is.list(fit_models) || !length(fit_models))
    stop("`fit_models` must be a non-empty list.", call. = FALSE)
  if (is.null(names(fit_models)) || any(!nzchar(names(fit_models))))
    names(fit_models) <- paste0("model", seq_along(fit_models))

  out <- lapply(seq_along(fit_models), function(i)
  {
    spec <- fit_models[[i]]
    if (inherits(spec, "formula"))
      spec <- list(formula = spec, conductance_model = conductance_model)
    if (!is.list(spec))
      stop("Each `fit_models` entry must be a list or formula.",
           call. = FALSE)
    if (is.null(spec$formula))
      spec$formula <- formula
    if (is.null(spec$conductance_model))
      spec$conductance_model <- conductance_model
    if (!inherits(spec$formula, "formula"))
      stop("Each candidate model `formula` must be a formula.",
           call. = FALSE)
    if (!inherits(spec$conductance_model,
                  c("terradish_conductance_model_factory",
                    "radish_conductance_model_factory")))
      stop("Each candidate `conductance_model` must be a conductance-model factory.",
           call. = FALSE)
    model <- spec$conductance_model(
      .covariance_power_rhs_formula(spec$formula), data$x)
    external_model <- .externalize_conductance_model(model)
    spec$parameter_names <- names(attr(external_model, "default",
                                       exact = TRUE))
    spec
  })
  names(out) <- names(fit_models)
  out
}

.covariance_power_spec_value <- function(spec, name, default)
{
  if (!is.null(spec[[name]]))
    spec[[name]]
  else
    default
}

.covariance_power_next_seed <- function()
{
  sample.int(.Machine$integer.max, 1L)
}

.covariance_power_indices <- function(data, n_focal, sample_size, strategy,
                                      seed)
{
  control <- list(n_landmarks = as.integer(sample_size),
                  method = strategy,
                  seed = seed)
  .landmark_indices(data, n_focal, control)
}

.covariance_power_subset_graph <- function(data, index)
{
  subset_data <- data
  subset_data$demes <- data$demes[index]
  if (!is.null(data$rhs))
    subset_data$rhs <- data$rhs[, index, drop = FALSE]
  subset_data
}

.covariance_power_se <- function(fit)
{
  theta <- coef(fit)
  out <- rep(NA_real_, length(theta))
  names(out) <- names(theta)
  vcov <- tryCatch(.safe_hessian_inverse(fit$fit$hessian),
                   error = function(e) NULL)
  if (is.null(vcov))
    return(out)
  se <- sqrt(pmax(diag(vcov), 0))
  names(se) <- names(theta)
  out[names(se)] <- se
  out
}

.covariance_power_result_row <- function(scenario_id, strategy, sample_size,
                                         design, simulation, model, status,
                                         elapsed, n_focal, n_pairs, fit,
                                         true_conductance,
                                         conductance_cor_threshold)
{
  ok <- identical(status, "OK") && inherits(fit, c("terradish", "radish")) &&
    !is.null(fit$mle$theta)
  fitted_conductance <- if (ok)
    tryCatch(fit$submodels$f(fit$mle$theta)$conductance,
             error = function(e) rep(NA_real_, length(true_conductance)))
  else
    rep(NA_real_, length(true_conductance))

  # When the fit fails for every replicate in a cell, fitted_conductance is
  # all NA and stats::cor(use = "complete.obs") raises "no complete element
  # pairs". tryCatch demotes that to NA so a failed cell still produces a
  # complete result row instead of aborting the entire scenario.
  safe_cor <- function(x, y, method = "pearson") {
    tryCatch(
      suppressWarnings(stats::cor(x, y, use = "complete.obs",
                                  method = method)),
      error = function(e) NA_real_
    )
  }
  conductance_cor      <- safe_cor(true_conductance, fitted_conductance)
  conductance_spearman <- safe_cor(true_conductance, fitted_conductance,
                                   method = "spearman")

  aicc <- if (ok && n_focal > fit$df + 1)
    -2 * fit$loglik + 2 * fit$df * (n_focal / (n_focal - fit$df - 1))
  else
    NA_real_

  data.frame(
    scenario_id = scenario_id,
    strategy = strategy,
    sample_size = as.integer(sample_size),
    design = as.integer(design),
    simulation = as.integer(simulation),
    model = model,
    n_focal = as.integer(n_focal),
    n_pairs = as.numeric(n_pairs),
    status = status,
    fit_ok = ok,
    elapsed = elapsed,
    loglik = if (ok) fit$loglik else NA_real_,
    df = if (ok) fit$df else NA_real_,
    AIC = if (ok) fit$aic else NA_real_,
    AICc = aicc,
    Delta_AIC = NA_real_,
    Delta_AICc = NA_real_,
    selected_AIC = FALSE,
    selected_AICc = FALSE,
    conductance_cor = conductance_cor,
    conductance_spearman = conductance_spearman,
    conductance_recovered = isTRUE(is.finite(conductance_cor) &&
                                     conductance_cor >=
                                       conductance_cor_threshold),
    stringsAsFactors = FALSE
  )
}

.covariance_power_mark_selection <- function(x)
{
  ok_aic <- is.finite(x$AIC)
  if (any(ok_aic))
  {
    min_aic <- min(x$AIC[ok_aic])
    x$Delta_AIC[ok_aic] <- x$AIC[ok_aic] - min_aic
    x$selected_AIC[ok_aic] <- x$AIC[ok_aic] == min_aic
  }

  ok_aicc <- is.finite(x$AICc)
  if (any(ok_aicc))
  {
    min_aicc <- min(x$AICc[ok_aicc])
    x$Delta_AICc[ok_aicc] <- x$AICc[ok_aicc] - min_aicc
    x$selected_AICc[ok_aicc] <- x$AICc[ok_aicc] == min_aicc
  }

  x
}

.covariance_power_parameter_row <- function(fit_row, parameter, true,
                                            estimate, se, alpha)
{
  zcrit <- stats::qnorm(1 - alpha / 2)
  lower <- estimate - zcrit * se
  upper <- estimate + zcrit * se
  p_value <- if (is.finite(estimate) && is.finite(se) && se > 0)
    2 * stats::pnorm(-abs(estimate / se))
  else
    NA_real_
  significant <- is.finite(p_value) && p_value < alpha
  correct_sign <- is.finite(estimate) && sign(estimate) == sign(true)
  recovered <- isTRUE(significant && correct_sign)
  covers_true <- is.finite(lower) && is.finite(upper) &&
    lower <= true && true <= upper

  data.frame(
    scenario_id = fit_row$scenario_id,
    strategy = fit_row$strategy,
    sample_size = fit_row$sample_size,
    design = fit_row$design,
    simulation = fit_row$simulation,
    model = fit_row$model,
    parameter = parameter,
    status = fit_row$status,
    fit_ok = fit_row$fit_ok,
    true = true,
    estimate = estimate,
    se = se,
    lower = lower,
    upper = upper,
    p_value = p_value,
    significant = significant,
    correct_sign = correct_sign,
    recovered = recovered,
    covers_true = covers_true,
    stringsAsFactors = FALSE
  )
}

.covariance_power_summary <- function(results, conductance_cor_threshold)
{
  if (!nrow(results))
    return(data.frame())
  groups <- split(results,
                  interaction(results$strategy, results$sample_size,
                              results$model, drop = TRUE))
  out <- lapply(groups, function(x)
  {
    data.frame(
      strategy = x$strategy[[1]],
      sample_size = x$sample_size[[1]],
      model = x$model[[1]],
      n_replicates = nrow(x),
      fit_rate = mean(x$fit_ok),
      conductance_power = mean(x$conductance_recovered),
      mean_conductance_cor = .covariance_power_mean(x$conductance_cor),
      mean_conductance_spearman = .covariance_power_mean(
        x$conductance_spearman),
      selected_AIC_rate = mean(x$selected_AIC),
      selected_AICc_rate = mean(x$selected_AICc),
      mean_elapsed = .covariance_power_mean(x$elapsed),
      conductance_cor_threshold = conductance_cor_threshold,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$sample_size, out$strategy, out$model), , drop = FALSE]
}

.covariance_power_add_parameter_power <- function(summary, parameter_results)
{
  summary$all_parameter_power <- NA_real_
  if (!nrow(summary) || !nrow(parameter_results))
    return(summary)

  replicate_groups <- split(
    parameter_results,
    interaction(parameter_results$strategy,
                parameter_results$sample_size,
                parameter_results$model,
                parameter_results$design,
                parameter_results$simulation,
                drop = TRUE)
  )
  replicate_power <- do.call(rbind, lapply(replicate_groups, function(x)
  {
    data.frame(
      strategy = x$strategy[[1]],
      sample_size = x$sample_size[[1]],
      model = x$model[[1]],
      all_recovered = all(x$recovered),
      stringsAsFactors = FALSE
    )
  }))

  scenario_groups <- split(
    replicate_power,
    interaction(replicate_power$strategy,
                replicate_power$sample_size,
                replicate_power$model,
                drop = TRUE)
  )
  scenario_power <- do.call(rbind, lapply(scenario_groups, function(x)
  {
    data.frame(
      strategy = x$strategy[[1]],
      sample_size = x$sample_size[[1]],
      model = x$model[[1]],
      all_parameter_power = mean(x$all_recovered),
      stringsAsFactors = FALSE
    )
  }))

  key_summary <- paste(summary$strategy, summary$sample_size, summary$model,
                       sep = "\r")
  key_power <- paste(scenario_power$strategy, scenario_power$sample_size,
                     scenario_power$model, sep = "\r")
  idx <- match(key_summary, key_power)
  has_match <- !is.na(idx)
  summary$all_parameter_power[has_match] <-
    scenario_power$all_parameter_power[idx[has_match]]
  summary
}

.covariance_power_parameter_summary <- function(parameter_results)
{
  if (!nrow(parameter_results))
    return(data.frame())
  groups <- split(parameter_results,
                  interaction(parameter_results$strategy,
                              parameter_results$sample_size,
                              parameter_results$model,
                              parameter_results$parameter,
                              drop = TRUE))
  out <- lapply(groups, function(x)
  {
    fit_ok <- x$fit_ok
    estimate <- x$estimate
    true <- x$true[[1]]
    data.frame(
      strategy = x$strategy[[1]],
      sample_size = x$sample_size[[1]],
      model = x$model[[1]],
      parameter = x$parameter[[1]],
      true = true,
      n_replicates = nrow(x),
      fit_rate = mean(fit_ok),
      power = mean(x$recovered),
      correct_sign_rate = mean(x$correct_sign),
      significant_rate = mean(x$significant),
      coverage = mean(x$covers_true),
      mean_estimate = .covariance_power_mean(estimate),
      bias = .covariance_power_mean(estimate - true),
      rmse = sqrt(.covariance_power_mean((estimate - true)^2)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$sample_size, out$strategy, out$model, out$parameter),
      , drop = FALSE]
}

.covariance_power_mean <- function(x)
{
  if (!any(is.finite(x)))
    return(NA_real_)
  mean(x[is.finite(x)])
}
