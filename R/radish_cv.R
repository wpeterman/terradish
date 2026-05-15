.terradish_measurement_model <- function(model, subset = NULL)
{
  if (inherits(model, c("terradish_measurement_model",
                        "radish_measurement_model")))
  {
    subsetter <- attr(model, "subsetter", exact = TRUE)
    if (!is.null(subset) && is.function(subsetter))
      return(subsetter(subset))
    return(model)
  }

  if (!is.character(model) || length(model) != 1L)
    stop("`model` must be a measurement-model function or one of 'mlpe', 'wishart', or 'ls'")

  switch(tolower(model),
         mlpe = mlpe,
         wishart = generalized_wishart,
         generalized_wishart = generalized_wishart,
         ls = leastsquares,
         leastsquares = leastsquares,
         stop("Unknown measurement model: ", model))
}

.fit_terradish_with_fallback <- function(formula, data, measurement_model, nu = NULL,
                                         conductance_model = loglinear_conductance,
                                         theta = NULL, cores = 1L, dots = list(),
                                         response_matrix = NULL)
{
  fit_once <- function(optimizer)
  {
    args <- c(list(formula = formula,
                   data = data,
                   conductance_model = conductance_model,
                   measurement_model = measurement_model,
                   nu = nu,
                   optimizer = optimizer,
                   cores = cores),
              dots)
    if (!is.null(theta))
      args$theta <- theta
    eval_env <- list2env(list(gd_mat = response_matrix), parent = parent.frame())
    call <- as.call(c(list(as.name("terradish")), args))
    tryCatch(eval(call, envir = eval_env), error = identity)
  }

  fit <- fit_once("newton")
  if (inherits(fit, "error"))
    fit <- fit_once("bfgs")
  if (inherits(fit, "error"))
    stop("Could not optimize terradish model: ", conditionMessage(fit), call. = FALSE)

  fit
}

.cv_conductance_model_for_surface <- function(conductance_model, formula,
                                              surface, reference_model = NULL)
{
  rebuild_for_surface <- attr(conductance_model, "rebuild_for_surface",
                              exact = TRUE)
  if (isTRUE(attr(conductance_model, "requires_fixed_graph", exact = TRUE)) &&
      is.function(rebuild_for_surface))
  {
    rebuilt_factory <- function(formula, x)
      rebuild_for_surface(formula, surface, reference_model = reference_model)
    attrs <- attributes(conductance_model)
    attrs$rebuild_for_surface <- NULL
    attrs$requires_fixed_graph <- NULL
    attrs$class <- class(conductance_model)
    attributes(rebuilt_factory) <- attrs
    return(rebuilt_factory)
  }

  conductance_model
}

.result_dir_files <- function(path)
{
  stopifnot(dir.exists(path))
  all_dirs <- list.files(path, full.names = TRUE)
  all_results_path <- all_dirs[grepl("AllResults", basename(all_dirs))]
  if (length(all_results_path) != 1L)
    stop("Expected exactly one file matching 'AllResults' in `Results_dir`")

  all_results <- readRDS(all_results_path)
  out <- list(all_results = all_results,
              all_dirs = all_dirs)

  if (!is.null(all_results$pts))
    out$pts <- all_results$pts
  if (!is.null(all_results$covariates))
    out$covariates <- tryCatch(.as_spatraster(all_results$covariates),
                               error = function(e) all_results$covariates)
  if (!is.null(all_results$conductance_surface))
    out$conductance_surface <- tryCatch(.as_spatraster(all_results$conductance_surface),
                                        error = function(e) all_results$conductance_surface)
  if (!is.null(all_results$sim_genind))
    out$sim_genind <- all_results$sim_genind
  if (!is.null(all_results$effect_size))
    out$effect_size <- all_results$effect_size

  out
}

.match_radish_model_file <- function(all_dirs, radish_model)
{
  model_key <- tolower(radish_model)
  patterns <- switch(model_key,
                     wishart = c("--wishart", "--generalized_wishart"),
                     generalized_wishart = c("--generalized_wishart", "--wishart"),
                     mlpe = "--mlpe",
                     ls = c("--ls", "--leastsquares"),
                     leastsquares = c("--ls", "--leastsquares"),
                     stop("Unknown `radish_model`: ", radish_model))

  hits <- unique(unlist(lapply(patterns, function(pat) {
    all_dirs[grepl(pat, all_dirs, fixed = TRUE)]
  })))

  if (length(hits) != 1L)
    stop("Expected exactly one fitted model matching `radish_model` in `Results_dir`")

  hits
}

.align_effect_size <- function(est_eff, effect_size)
{
  out <- rep(NA_real_, length(est_eff))
  names(out) <- names(est_eff)

  if (is.null(effect_size))
    return(out)

  effect_size <- as.numeric(effect_size)
  effect_names <- names(effect_size)
  if (!is.null(effect_names))
  {
    keep <- intersect(names(est_eff), effect_names)
    out[keep] <- effect_size[keep]
  }
  else if (length(effect_size) == length(est_eff))
  {
    out[] <- effect_size
  }

  out
}

#' Rank fitted terradish models by information criterion
#'
#' Creates a model-selection table from fitted \code{terradish} models using
#' AIC, AICc, or BIC.
#'
#' @param mod_list List of fitted \code{terradish} models.
#' @param AICc Should second-order AIC be used instead of AIC?
#' @param BIC Should BIC be used instead of AIC?
#' @param mod_names Optional model names. By default the right-hand side of
#'   each fitted formula is used. For MLPE measurement models with additional
#'   pairwise covariates, the default appends \code{[mlpe:n]}, where \code{n}
#'   is the number of added pairwise covariate columns.
#' @param verbose Should the table be printed to the console?
#'
#' @return A data frame containing model ranks, parameter counts, information
#'   criterion values, delta values, weights, cumulative weights, and
#'   log-likelihoods.
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' keep <- 1:12
#' melip.Fst_small <- melip.Fst[keep, keep]
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
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
#' @export
aic_table <- function(mod_list, AICc = FALSE, BIC = FALSE, mod_names = NULL, verbose = FALSE)
{
  mod_dim_keys <- vapply(mod_list,
                         function(x) paste(x$dim, collapse = "|"),
                         character(1))
  if (length(unique(mod_dim_keys)) != 1L)
    stop("Models must be fit to the same number of focal points and graph size")

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
                      row.names = NULL)
    tab <- tab[order(tab$AIC), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$AIC_wt)
    tab$loglik <- mod_loglik[match(tab$model, mod_names)]
    tab[, 3:7] <- round(tab[, 3:7], digits = 4)
  }
  else if (isTRUE(AICc))
  {
    if (isTRUE(BIC))
      stop("Set only one of `AICc` or `BIC` to TRUE")

    mod_n <- vapply(mod_list, function(x) x$dim[["focal"]], numeric(1))
    mod_AICc <- -2 * mod_loglik + 2 * mod_df * (mod_n / (mod_n - mod_df - 1))
    delta <- mod_AICc - min(mod_AICc)
    wt <- exp(-0.5 * delta)
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      AICc = mod_AICc,
                      Delta_AICc = delta,
                      AICc_wt = wt / sum(wt),
                      row.names = NULL)
    tab <- tab[order(tab$AICc), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$AICc_wt)
    tab$loglik <- mod_loglik[match(tab$model, mod_names)]
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
                      row.names = NULL)
    tab <- tab[order(tab$BIC), , drop = FALSE]
    tab$Cum.wt <- cumsum(tab$BIC_wt)
    tab$loglik <- mod_loglik[match(tab$model, mod_names)]
    tab[, 3:8] <- round(tab[, 3:8], digits = 4)
  }

  if (isTRUE(verbose))
    print(tab, row.names = FALSE)
  else
    tab
}

.cv_model_name <- function(mod)
{
  .default_model_name(mod)
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

#' Inspect a saved terradish results directory
#'
#' Reads a terradish-style results directory and returns the discovered files
#' and any recognized saved objects, such as \code{AllResults} metadata,
#' covariates, focal points, or simulation effect sizes.
#'
#' @param Results_dir Full path to the top-level results directory.
#'
#' @return A named list containing the directory contents plus any recognized
#'   saved objects recovered from the \code{AllResults} file.
#'
#' @details
#' This helper is intended for downstream simulation or workflow code that
#' needs to inspect saved terradish outputs before extracting parameter tables
#' with \code{\link{terradish_parameters}}.
#'
#' @examples
#' \dontrun{
#' tmp <- tempfile()
#' dir.create(tmp)
#' saveRDS(list(effect_size = c(altitude = 0.2)),
#'         file.path(tmp, "AllResults_list.rds"))
#' terradish_results(tmp)
#' }
#'
#' @seealso \code{\link{terradish_parameters}}
#'
#' @export
terradish_results <- function(Results_dir)
{
  .result_dir_files(Results_dir)
}

#' Extract parameter estimates from saved terradish results
#'
#' Reads a saved fitted \code{terradish} model from a simulation results
#' directory
#' and returns a table of estimated coefficients and standard errors.
#'
#' @param Results_dir Full path to the top-level results directory.
#' @param model Which saved \code{terradish} model to read. One of
#'   \code{"wishart"}, \code{"generalized_wishart"}, \code{"mlpe"},
#'   \code{"ls"}, or \code{"leastsquares"}.
#' @param radish_model Deprecated alias for \code{model}.
#' @param save_table Should the parameter table be written to
#'   \code{Results_dir} as a CSV file?
#' @param conv Optional convergence flag or vector to append to the output.
#' @param ... Reserved for future use.
#'
#' @return A data frame of fitted and, when available, true effect sizes.
#'
#' @details
#' This helper is intended for saved-results workflows rather than core model
#' fitting. It assumes \code{Results_dir} contains a single fitted terradish
#' model matching \code{model} and a single file matching \code{"AllResults"}.
#' When present, named \code{effect_size} values from the \code{AllResults}
#' object are aligned to the fitted coefficient names. The exported
#' \code{\link{terradish_results}} helper can be used to inspect the recovered
#' directory contents before extracting parameter summaries.
#'
#' The returned standard errors are taken from
#' \code{sqrt(diag(summary(mod)$vcov))} when available, so they describe the
#' conductance coefficients of the saved fitted model rather than any
#' simulation-wide uncertainty summary.
#'
#' @seealso \code{\link{terradish_results}}
#'
#' @examples
#' \dontrun{
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#' covariates <- c(terra::scale(melip.altitude), terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#' fit <- terradish(melip.Fst ~ altitude + forestcover, surface,
#'               loglinear_conductance, leastsquares)
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' saveRDS(fit, file.path(tmp, "fit--ls.rds"))
#' saveRDS(list(effect_size = coef(fit)), file.path(tmp, "AllResults_list.rds"))
#' terradish_parameters(tmp, model = "ls", save_table = FALSE)
#' }
#'
#' @export
terradish_parameters <- function(Results_dir,
                                 model = "wishart",
                                 save_table = TRUE,
                                 conv = NULL,
                                 radish_model,
                                 ...)
{
  if (!missing(radish_model))
  {
    warning("`radish_model` is deprecated; use `model` instead.", call. = FALSE)
    model <- radish_model
  }
  params <- .result_dir_files(Results_dir)
  model_path <- .match_radish_model_file(params$all_dirs, model)
  mod <- readRDS(model_path)

  est_eff <- coef(mod)
  vcov_mat <- tryCatch(summary(mod)$vcov, error = function(e) NULL)
  se <- rep(NA_real_, length(est_eff))
  if (!is.null(vcov_mat))
    se <- sqrt(diag(vcov_mat))

  df <- data.frame(parameter = names(est_eff),
                   true_eff = .align_effect_size(est_eff, params$effect_size),
                   est_eff = unname(est_eff),
                   est_SE = unname(se),
                   stringsAsFactors = FALSE)

  if (!is.null(conv))
    df$converged <- conv

  if (isTRUE(save_table))
    write.csv(df,
              file = file.path(Results_dir, "parameter_estimates.csv"),
              row.names = FALSE)

  df
}

#' Legacy radish parameter-extraction wrapper
#'
#' Deprecated compatibility wrapper retained for older code that still calls
#' \code{radish_parameters()}.
#'
#' @param ... Arguments passed through the deprecated
#'   \code{radish_parameters()} compatibility wrapper to
#'   \code{\link{terradish_parameters}}.
#' @name legacy_radish_parameters_wrapper
#' @keywords internal
NULL

#' @rdname legacy_radish_parameters_wrapper
#' @export
radish_parameters <- function(...)
{
  .terradish_deprecate("radish_parameters", "terradish_parameters")
  .terradish_forward_call(match.call(), "terradish_parameters")
}

#' Cross validation for terradish models
#'
#' Randomly splits focal points into training and test sets, fits a
#' \code{terradish} model on the training set, and evaluates the fitted
#' coefficients on the test set.
#'
#' @param pts Focal point coordinates as a \code{terra::SpatVector}, matrix, or
#'   data frame with \code{x}/\code{y} columns.
#' @param covariates Spatial covariates as a \code{terra::SpatRaster}.
#' @param fmla Formula describing the model to assess. The left-hand side must
#'   evaluate to a genetic distance matrix in the calling environment.
#' @param model Measurement model, either as a function or one of
#'   \code{"mlpe"}, \code{"wishart"}, or \code{"ls"}.
#' @param nu Number of genetic markers, passed to the measurement model.
#' @param prop_train Proportion of focal points assigned to the training set.
#' @param seed Optional random seed used for the split. If \code{NULL}, one is
#'   generated and returned.
#' @param fit_full Should the model also be fit to all focal points?
#' @param directions Neighborhood definition passed to
#'   \code{\link{conductance_surface}}.
#' @param conductance_model Conductance-model factory used for training and
#'   held-out evaluation. Defaults to \code{\link{loglinear_conductance}}.
#'   Factories that require a fixed graph, such as
#'   \code{\link{gaussian_smoothed_loglinear_conductance}}, are rebuilt on each
#'   cross-validation surface when they provide a \code{rebuild_for_surface}
#'   method.
#' @param cores Number of worker processes to use in downstream \code{terradish}
#'   fits and grid evaluation.
#' @param ... Additional arguments passed to \code{\link{terradish}}, such as
#'   \code{control}.
#'
#' @return A named list containing the training fit, test-set log-likelihood,
#'   train/test indices, and optionally the full-data fit.
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' keep <- 1:12
#' melip.Fst_small <- melip.Fst[keep, keep]
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#'
#' cv_fit <- terradish_cv(melip.coords[keep], covariates,
#'                     melip.Fst_small ~ altitude + forestcover,
#'                     model = "ls",
#'                     prop_train = 0.75,
#'                     seed = 1,
#'                     fit_full = FALSE,
#'                     control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' cv_fit$cv_loglik
#'
#' @export
terradish_cv <- function(pts,
                      covariates,
                      fmla,
                      model = "mlpe",
                      nu = NULL,
                      prop_train = 0.8,
                      seed = NULL,
                      fit_full = TRUE,
                      directions = 8,
                      conductance_model = loglinear_conductance,
                      cores = 1L,
                      ...)
{
  stopifnot(length(prop_train) == 1, is.numeric(prop_train), prop_train > 0, prop_train < 1)
  stopifnot(length(cores) == 1, is.numeric(cores), cores >= 1)
  stopifnot(inherits(conductance_model,
                    c("terradish_conductance_model_factory",
                      "radish_conductance_model_factory")))

  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  set.seed(seed)

  covariates <- .as_spatraster(covariates)
  pts <- .coords_matrix(pts, covariates)
  total_ <- nrow(pts)
  train_n <- floor(total_ * prop_train)
  if (train_n < 2 || (total_ - train_n) < 2)
    stop("`prop_train` must leave at least two focal points in both train and test sets")

  train_ <- sort(sample.int(total_, train_n))
  test_ <- setdiff(seq_len(total_), train_)

  fmla <- as.formula(fmla)
  trm <- terms(fmla)
  response <- attr(trm, "response")
  if (response != 1L)
    stop("'fmla' must have a genetic distance matrix on the left-hand side")

  vars <- attr(trm, "variables")
  gd_mat <- eval(vars[[response + 1L]], envir = parent.frame(), enclos = environment(fmla))
  rhs <- attr(trm, "term.labels")
  fmla_radish <- if (length(rhs))
    reformulate(rhs, response = "gd_mat")
  else
    as.formula("gd_mat ~ 1")
  environment(fmla_radish) <- environment()

  measurement_model_train <- .terradish_measurement_model(model, subset = train_)
  dots <- list(...)
  training_dots <- dots
  approximation <- if ("approximation" %in% names(dots))
    dots$approximation
  else
    "none"
  approximation_control <- dots$approximation_control
  if (identical(approximation, "coarse_raster"))
  {
    training_dots$approximation <- NULL
    training_dots$approximation_control <- NULL
  }
  gd_mat_full <- as.matrix(gd_mat)

  gd_mat <- gd_mat_full[train_, train_, drop = FALSE]
  training_surface <- conductance_surface(covariates, pts[train_, , drop = FALSE],
                                          directions = directions)
  training_conductance_model <- .cv_conductance_model_for_surface(
    conductance_model, fmla_radish, training_surface
  )
  fit <- .fit_terradish_with_fallback(fmla_radish, training_surface,
                                      measurement_model = measurement_model_train,
                                      nu = nu,
                                      conductance_model = training_conductance_model,
                                      cores = cores,
                                      dots = training_dots,
                                      response_matrix = gd_mat)

  if (length(coef(fit)) < 1L)
    stop("Training fit produced no coefficients, so cross-validation cannot continue")

  gd_mat <- gd_mat_full[test_, test_, drop = FALSE]
  test_surface <- conductance_surface(covariates, pts[test_, , drop = FALSE],
                                      directions = directions)
  test_conductance_model <- .cv_conductance_model_for_surface(
    conductance_model, fmla_radish, test_surface,
    reference_model = fit$submodels$f_internal
  )
  measurement_model_test <- .terradish_measurement_model(model, subset = test_)
  ll <- terradish_grid(theta = matrix(coef(fit), nrow = 1),
                    formula = fmla_radish,
                    data = test_surface,
                    conductance_model = test_conductance_model,
                    measurement_model = measurement_model_test,
                    nu = nu,
                    cores = cores,
                    approximation = approximation,
                    approximation_control = approximation_control)

  out <- list(train_mod = fit,
              cv_loglik = ll$loglik,
              seed = seed,
              train_index = train_,
              test_index = test_)

  if (isTRUE(fit_full))
  {
    gd_mat <- gd_mat_full
    full_surface <- conductance_surface(covariates, pts, directions = directions)
    full_conductance_model <- .cv_conductance_model_for_surface(
      conductance_model, fmla_radish, full_surface,
      reference_model = fit$submodels$f_internal
    )
    measurement_model_full <- .terradish_measurement_model(model)
    out$full_mod <- .fit_terradish_with_fallback(fmla_radish, full_surface,
                                                 measurement_model = measurement_model_full,
                                                 nu = nu,
                                                 conductance_model = full_conductance_model,
                                                 theta = fit$mle$theta,
                                                 cores = cores,
                                                 dots = training_dots,
                                                 response_matrix = gd_mat)
  }

  out
}

#' Repeated cross validation for terradish models
#'
#' Repeats \code{\link{terradish_cv}} across multiple random splits and
#' summarizes held-out log-likelihood across replicates.
#'
#' @param pts Focal point coordinates as a \code{terra::SpatVector}, matrix, or
#'   data frame with \code{x}/\code{y} columns.
#' @param covariates Spatial covariates as a \code{terra::SpatRaster}.
#' @param fmla Formula describing the model to assess. The left-hand side must
#'   evaluate to a genetic distance matrix in the calling environment.
#' @param model Measurement model, either as a function or one of
#'   \code{"mlpe"}, \code{"wishart"}, or \code{"ls"}.
#' @param nu Number of genetic markers, passed to the measurement model.
#' @param prop_train Proportion of focal points assigned to the training set.
#' @param n_reps Number of repeated train/test splits to evaluate.
#' @param seeds Optional integer vector of seeds to use for each replicate. If
#'   supplied, \code{n_reps} is inferred from \code{length(seeds)}.
#' @param fit_full Should each replicate also fit a model to all focal points?
#' @param keep_fits Should the per-replicate \code{terradish_cv()} outputs be
#'   retained in the returned object?
#' @param directions Neighborhood definition passed to
#'   \code{\link{conductance_surface}}.
#' @param conductance_model Conductance-model factory used for each replicate.
#'   See \code{\link{terradish_cv}}.
#' @param cores Number of worker processes to use in downstream \code{terradish}
#'   fits and grid evaluation.
#' @param ... Additional arguments passed to \code{\link{terradish}}, such as
#'   \code{control}.
#'
#' @return A list containing a per-replicate summary table, mean and standard
#'   deviation of held-out log-likelihood, and optionally the individual
#'   cross-validation fits.
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' keep <- 1:12
#' melip.Fst_small <- melip.Fst[keep, keep]
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#'
#' cv_rep <- terradish_cv_replicates(melip.coords[keep], covariates,
#'                                   melip.Fst_small ~ altitude + forestcover,
#'                                   model = "ls",
#'                                   n_reps = 2,
#'                                   fit_full = FALSE,
#'                                   control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' cv_rep$mean_loglik
#'
#' @export
terradish_cv_replicates <- function(pts,
                                    covariates,
                                    fmla,
                                    model = "mlpe",
                                    nu = NULL,
                                    prop_train = 0.8,
                                    n_reps = 5L,
                                    seeds = NULL,
                                    fit_full = FALSE,
                                    keep_fits = FALSE,
                                    directions = 8,
                                    conductance_model = loglinear_conductance,
                                    cores = 1L,
                                    ...)
{
  eval_env <- parent.frame()

  if (is.null(seeds))
  {
    stopifnot(length(n_reps) == 1L, is.numeric(n_reps), n_reps >= 1)
    n_reps <- as.integer(n_reps)
    seeds <- sample.int(.Machine$integer.max, n_reps)
  }
  else
  {
    if (!is.numeric(seeds) || anyNA(seeds))
      stop("`seeds` must be a numeric vector without missing values")
    seeds <- as.integer(seeds)
    n_reps <- length(seeds)
  }

  fits <- vector("list", n_reps)
  summary_tab <- data.frame(
    replicate = seq_len(n_reps),
    seed = seeds,
    cv_loglik = NA_real_
  )

  for (i in seq_len(n_reps))
  {
    fit_i <- eval(
      as.call(c(
        list(as.name("terradish_cv")),
        list(
          pts = pts,
          covariates = covariates,
          fmla = fmla,
          model = model,
          nu = nu,
          prop_train = prop_train,
          seed = seeds[i],
          fit_full = fit_full,
          directions = directions,
          conductance_model = conductance_model,
          cores = cores
        ),
        list(...)
      )),
      envir = eval_env
    )

    summary_tab$cv_loglik[i] <- fit_i$cv_loglik
    if (isTRUE(keep_fits))
      fits[[i]] <- fit_i
  }

  out <- list(
    summary = summary_tab,
    mean_loglik = mean(summary_tab$cv_loglik),
    sd_loglik = sd(summary_tab$cv_loglik),
    seeds = seeds
  )

  if (isTRUE(keep_fits))
    out$fits <- fits

  class(out) <- "terradish_cv_replicates"
  out
}

#' Methods for repeated terradish cross-validation results
#'
#' S3 methods for working with objects returned by
#' \code{\link{terradish_cv_replicates}}.
#'
#' @name terradish_cv_replicates_methods
#' @title Methods for repeated terradish cross-validation results
#'
#' @param x A repeated cross-validation object returned by
#'   \code{\link{terradish_cv_replicates}}.
#' @param object A repeated cross-validation object returned by
#'   \code{\link{terradish_cv_replicates}}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments passed through to generic methods.
#'
#' @return
#' \itemize{
#'   \item \code{print()} returns its input invisibly.
#'   \item \code{summary()} returns an object of class
#'     \code{"summary.terradish_cv_replicates"}.
#' }
#'
#' @export
print.terradish_cv_replicates <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Repeated terradish cross-validation\n")
  cat("Replicates:", nrow(x$summary), "\n")
  cat("Mean held-out loglikelihood:", format(x$mean_loglik, digits = digits), "\n")
  cat("SD held-out loglikelihood:", format(x$sd_loglik, digits = digits), "\n\n")
  print.default(x$summary, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname terradish_cv_replicates_methods
#' @export
summary.terradish_cv_replicates <- function(object, ...)
{
  tab <- object$summary
  out <- list(
    replicates = nrow(tab),
    mean_loglik = object$mean_loglik,
    sd_loglik = object$sd_loglik,
    min_loglik = min(tab$cv_loglik),
    max_loglik = max(tab$cv_loglik),
    summary = tab,
    seeds = object$seeds
  )
  class(out) <- "summary.terradish_cv_replicates"
  out
}

#' @rdname terradish_cv_replicates_methods
#' @export
print.summary.terradish_cv_replicates <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Repeated terradish cross-validation summary\n")
  cat("Replicates:", x$replicates, "\n")
  cat("Mean held-out loglikelihood:", format(x$mean_loglik, digits = digits), "\n")
  cat("SD held-out loglikelihood:", format(x$sd_loglik, digits = digits), "\n")
  cat("Range:", format(x$min_loglik, digits = digits), "to",
      format(x$max_loglik, digits = digits), "\n\n")
  print.default(x$summary, row.names = FALSE, ...)
  invisible(x)
}

#' Legacy radish cross-validation wrapper
#'
#' Deprecated compatibility wrapper retained for older code that still calls
#' \code{radish_cv()}.
#'
#' @param ... Additional arguments passed to \code{\link{terradish}}.
#' @name legacy_radish_cv_wrapper
#' @keywords internal
NULL

#' @rdname legacy_radish_cv_wrapper
#' @export
radish_cv <- function(...)
{
  .terradish_deprecate("radish_cv", "terradish_cv")
  .terradish_forward_call(match.call(), "terradish_cv")
}

# Internal compatibility aliases retained for housekeeping-sized refactors.
.radish_measurement_model <- .terradish_measurement_model
.fit_radish_with_fallback <- .fit_terradish_with_fallback

#' Summarize terradish cross-validation models
#'
#' Creates a ranked log-likelihood table from fitted cross-validation results
#' and, optionally, an information-criterion table from the full-data fits.
#'
#' @param cv_list A list of objects returned by \code{\link{terradish_cv}}.
#' @param cv_names Optional model names. By default the right-hand side of the
#'   fitted training formula is used.
#' @param aic Should an information-criterion table also be computed from
#'   \code{full_mod}?
#' @param AICc Should \code{aic = TRUE} use second-order Akaike's Information
#'   Criterion instead of AIC?
#' @param BIC Should \code{aic = TRUE} use BIC instead of AIC?
#' @param ... Reserved for future use.
#'
#' @return Either a ranked cross-validation table or a list containing the
#'   cross-validation table and information-criterion table.
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' keep <- 1:12
#' melip.Fst_small <- melip.Fst[keep, keep]
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#' surface_small <- conductance_surface(covariates, melip.coords[keep], directions = 8)
#'
#' fit1 <- terradish(melip.Fst_small ~ altitude, surface_small,
#'                loglinear_conductance, leastsquares,
#'                control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' fit2 <- terradish(melip.Fst_small ~ altitude + forestcover, surface_small,
#'                loglinear_conductance, leastsquares,
#'                control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#'
#' cv_model_selection(list(
#'   list(train_mod = fit1, cv_loglik = fit1$loglik, full_mod = fit1),
#'   list(train_mod = fit2, cv_loglik = fit2$loglik, full_mod = fit2)
#' ), aic = TRUE, AICc = TRUE)
#'
#' @export
cv_model_selection <- function(cv_list,
                               cv_names = NULL,
                               aic = FALSE,
                               AICc = FALSE,
                               BIC = FALSE,
                               ...)
{
  if (is.null(cv_names))
  {
    cv_names <- vapply(cv_list,
                       function(x) .cv_model_name(x$train_mod),
                       character(1))
  }

  mod_df <- vapply(cv_list, function(x) x$train_mod$df, numeric(1))
  mod_loglik <- vapply(cv_list, function(x) x$cv_loglik, numeric(1))
  delta_ll <- mod_loglik - max(mod_loglik)

  cv_tab <- data.frame(model = cv_names,
                       K = mod_df,
                       loglik = mod_loglik,
                       Delta_LL = delta_ll,
                       row.names = NULL)
  cv_tab <- cv_tab[order(-cv_tab$loglik), , drop = FALSE]

  if (!isTRUE(aic))
    return(cv_tab)

  mod_list <- lapply(cv_list, `[[`, "full_mod")
  if (any(vapply(mod_list, is.null, logical(1))))
    stop("All elements of `cv_list` must contain `full_mod` when `aic = TRUE`")

  mod_names <- vapply(mod_list, .cv_model_name, character(1))
  list(loglik_tab = cv_tab,
       AIC_tab = aic_table(mod_list,
                           AICc = AICc,
                           BIC = BIC,
                           mod_names = mod_names))
}
