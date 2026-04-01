.radish_measurement_model <- function(model)
{
  if (inherits(model, "radish_measurement_model"))
    return(model)

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

.fit_radish_with_fallback <- function(formula, data, measurement_model, nu = NULL,
                                      theta = NULL, cores = 1L, dots = list(),
                                      response_matrix = NULL)
{
  fit_once <- function(optimizer)
  {
    args <- c(list(formula = formula,
                   data = data,
                   conductance_model = loglinear_conductance,
                   measurement_model = measurement_model,
                   nu = nu,
                   optimizer = optimizer,
                   cores = cores),
              dots)
    if (!is.null(theta))
      args$theta <- theta
    eval_env <- list2env(list(gd_mat = response_matrix), parent = parent.frame())
    call <- as.call(c(list(as.name("radish")), args))
    tryCatch(eval(call, envir = eval_env), error = identity)
  }

  fit <- fit_once("newton")
  if (inherits(fit, "error"))
    fit <- fit_once("bfgs")
  if (inherits(fit, "error"))
    stop("Could not optimize radish model: ", conditionMessage(fit), call. = FALSE)

  fit
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

.aic_table <- function(mod_list, AICc = FALSE, BIC = FALSE, mod_names = NULL, verbose = FALSE)
{
  mod_dim_keys <- vapply(mod_list,
                         function(x) paste(x$dim, collapse = "|"),
                         character(1))
  if (length(unique(mod_dim_keys)) != 1L)
    stop("Models must be fit to the same number of focal points and graph size")

  if (is.null(mod_names))
    mod_names <- vapply(mod_list, function(x) as.character(x$formula)[2], character(1))

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
  if (!is.null(mod$call$formula))
    return(paste(deparse(mod$call$formula[[3]]), collapse = ""))
  as.character(mod$formula)[2]
}

#' Extract parameter estimates from saved radish results
#'
#' Reads a saved fitted \code{radish} model from a simulation results directory
#' and returns a table of estimated coefficients and standard errors.
#'
#' @param Results_dir Full path to the top-level results directory.
#' @param radish_model Which saved \code{radish} model to read. One of
#'   \code{"wishart"}, \code{"generalized_wishart"}, \code{"mlpe"},
#'   \code{"ls"}, or \code{"leastsquares"}.
#' @param save_table Should the parameter table be written to
#'   \code{Results_dir} as a CSV file?
#' @param conv Optional convergence flag or vector to append to the output.
#' @param ... Reserved for future use.
#'
#' @return A data frame of fitted and, when available, true effect sizes.
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
#' fit <- radish(melip.Fst ~ altitude + forestcover, surface,
#'               loglinear_conductance, leastsquares)
#'
#' tmp <- tempfile()
#' dir.create(tmp)
#' saveRDS(fit, file.path(tmp, "fit--ls.rds"))
#' saveRDS(list(effect_size = coef(fit)), file.path(tmp, "AllResults_list.rds"))
#' radish_parameters(tmp, radish_model = "ls", save_table = FALSE)
#' }
#'
#' @export
radish_parameters <- function(Results_dir,
                              radish_model = "wishart",
                              save_table = TRUE,
                              conv = NULL,
                              ...)
{
  params <- .result_dir_files(Results_dir)
  model_path <- .match_radish_model_file(params$all_dirs, radish_model)
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

#' Cross validation for radish models
#'
#' Randomly splits focal points into training and test sets, fits a
#' \code{radish} model on the training set, and evaluates the fitted
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
#' @param cores Number of worker processes to use in downstream \code{radish}
#'   fits and grid evaluation.
#' @param ... Additional arguments passed to \code{\link{radish}}, such as
#'   \code{control}. Landmark approximation arguments
#'   (\code{approximation} and \code{approximation_control}) are also forwarded
#'   to the held-out \code{\link{radish_grid}} evaluation.
#'   \code{approximation = "coarse_raster"} is currently applied only to the
#'   held-out grid evaluation, while the training fit remains exact.
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
#' cv_fit <- radish_cv(melip.coords[keep], covariates,
#'                     melip.Fst_small ~ altitude + forestcover,
#'                     model = "ls",
#'                     prop_train = 0.75,
#'                     seed = 1,
#'                     fit_full = FALSE,
#'                     control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' cv_fit$cv_loglik
#'
#' @export
radish_cv <- function(pts,
                      covariates,
                      fmla,
                      model = "mlpe",
                      nu = NULL,
                      prop_train = 0.8,
                      seed = NULL,
                      fit_full = TRUE,
                      directions = 8,
                      cores = 1L,
                      ...)
{
  stopifnot(length(prop_train) == 1, is.numeric(prop_train), prop_train > 0, prop_train < 1)
  stopifnot(length(cores) == 1, is.numeric(cores), cores >= 1)

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

  measurement_model <- .radish_measurement_model(model)
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
  fit <- .fit_radish_with_fallback(fmla_radish, training_surface,
                                   measurement_model = measurement_model,
                                   nu = nu,
                                   cores = cores,
                                   dots = training_dots,
                                   response_matrix = gd_mat)

  if (length(coef(fit)) < 1L)
    stop("Training fit produced no coefficients, so cross-validation cannot continue")

  gd_mat <- gd_mat_full[test_, test_, drop = FALSE]
  test_surface <- conductance_surface(covariates, pts[test_, , drop = FALSE],
                                      directions = directions)
  ll <- radish_grid(theta = matrix(coef(fit), nrow = 1),
                    formula = fmla_radish,
                    data = test_surface,
                    conductance_model = loglinear_conductance,
                    measurement_model = measurement_model,
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
    out$full_mod <- .fit_radish_with_fallback(fmla_radish, full_surface,
                                              measurement_model = measurement_model,
                                              nu = nu,
                                              theta = fit$mle$theta,
                                              cores = cores,
                                              dots = training_dots,
                                              response_matrix = gd_mat)
  }

  out
}

#' Summarize radish cross-validation models
#'
#' Creates a ranked log-likelihood table from fitted cross-validation results
#' and, optionally, an AIC table from the full-data fits.
#'
#' @param cv_list A list of objects returned by \code{\link{radish_cv}}.
#' @param cv_names Optional model names. By default the right-hand side of the
#'   fitted training formula is used.
#' @param aic Should an AIC table also be computed from \code{full_mod}?
#' @param ... Reserved for future use.
#'
#' @return Either a ranked cross-validation table or a list containing the
#'   cross-validation table and AIC table.
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
#' fit1 <- radish(melip.Fst_small ~ altitude, surface_small,
#'                loglinear_conductance, leastsquares,
#'                control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#' fit2 <- radish(melip.Fst_small ~ altitude + forestcover, surface_small,
#'                loglinear_conductance, leastsquares,
#'                control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
#'
#' cv_model_selection(list(
#'   list(train_mod = fit1, cv_loglik = fit1$loglik, full_mod = fit1),
#'   list(train_mod = fit2, cv_loglik = fit2$loglik, full_mod = fit2)
#' ))
#'
#' @export
cv_model_selection <- function(cv_list,
                               cv_names = NULL,
                               aic = FALSE,
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
       AIC_tab = .aic_table(mod_list, mod_names = mod_names))
}
