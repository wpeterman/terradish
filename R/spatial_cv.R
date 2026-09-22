.terradish_restore_seed <- function(seed)
{
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (old_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  if (!is.null(seed))
    set.seed(seed)
  function()
  {
    if (old_exists)
      assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}

.terradish_fold_coordinates <- function(pts, projected, require_projected = TRUE)
{
  if (inherits(pts, "SpatVector"))
  {
    if (is.null(projected) && isTRUE(require_projected))
    {
      if (!nzchar(terra::crs(pts)))
        stop("`pts` has no coordinate reference system. Set `projected` explicitly.",
             call. = FALSE)
      projected <- !terra::is.lonlat(pts)
    }
    coords <- terra::crds(pts)
  }
  else
  {
    coords <- as.matrix(pts)
    if (is.null(projected) && isTRUE(require_projected))
      stop("Set `projected = TRUE` or `FALSE` when `pts` is not a SpatVector.",
           call. = FALSE)
  }

  if (isTRUE(require_projected) && !isTRUE(projected))
    stop("Spatial folds require projected planar coordinates. Project `pts` before clustering.",
         call. = FALSE)
  if (!is.numeric(coords) || ncol(coords) < 2L || anyNA(coords[, 1:2, drop = FALSE]))
    stop("`pts` must provide two finite coordinate columns.", call. = FALSE)
  coords[, 1:2, drop = FALSE]
}

#' Construct reproducible cross-validation folds
#'
#' Creates either spatially blocked folds by k-means clustering projected site
#' coordinates or balanced random folds. Fold labels can also be supplied
#' directly to \code{\link{terradish_cv_folds}} when study-specific blocks are
#' more defensible than an automated partition.
#'
#' @param pts Focal-site coordinates as a \code{terra::SpatVector}, matrix, or
#'   data frame.
#' @param k Number of folds. Must be at least 2 and smaller than the number of
#'   sites.
#' @param method Fold construction method. \code{"spatial_kmeans"} applies
#'   k-means to the unscaled projected coordinates. \code{"random"} assigns
#'   sites to balanced random folds.
#' @param seed Optional random seed.
#' @param projected Logical indicating whether the coordinates are projected.
#'   For spatial k-means with a \code{SpatVector}, \code{NULL} uses its
#'   coordinate reference system. For spatial k-means with a matrix or data
#'   frame, this argument must be supplied explicitly. It is ignored for random
#'   folds.
#' @param nstart Number of random starts used by k-means.
#'
#' @details
#' Coordinates are deliberately not standardized before spatial clustering.
#' Consequently, Euclidean separation is measured in the map units of the
#' projected coordinate reference system. Automated k-means blocks are a
#' reproducible diagnostic, not a substitute for ecological blocking based on
#' watersheds, barriers, sampling campaigns, or other design information.
#'
#' @return An integer vector of fold labels with attributes recording the
#'   construction method, seed, and number of folds.
#'
#' @examples
#' data(melip)
#' coords <- terra::unwrap(melip.coords)
#' coords_projected <- terra::project(coords, "EPSG:5070")
#' folds <- terradish_folds(coords_projected, k = 3, seed = 1)
#' table(folds)
#'
#' @export
terradish_folds <- function(pts, k = 5L,
                            method = c("spatial_kmeans", "random"),
                            seed = 1L, projected = NULL, nstart = 25L)
{
  method <- match.arg(method)
  coords <- .terradish_fold_coordinates(
    pts, projected, require_projected = identical(method, "spatial_kmeans"))
  k <- as.integer(k)[1]
  if (is.na(k) || k < 2L || k >= nrow(coords))
    stop("`k` must be at least 2 and smaller than the number of sites.", call. = FALSE)
  nstart <- as.integer(nstart)[1]
  if (is.na(nstart) || nstart < 1L)
    stop("`nstart` must be a positive integer.", call. = FALSE)

  restore <- .terradish_restore_seed(seed)
  on.exit(restore(), add = TRUE)
  fold <- if (identical(method, "spatial_kmeans"))
    stats::kmeans(coords, centers = k, nstart = nstart)$cluster
  else
    sample(rep(seq_len(k), length.out = nrow(coords)))

  fold <- as.integer(fold)
  attr(fold, "method") <- method
  attr(fold, "seed") <- seed
  attr(fold, "k") <- k
  fold
}

.terradish_subset_graph_focal <- function(data, index)
{
  out <- data
  out$demes <- data$demes[index]
  if (!is.null(data$rhs))
    out$rhs <- data$rhs[, index, drop = FALSE]
  out
}

.terradish_cv_formula <- function(formula)
{
  trm <- terms(formula)
  if (attr(trm, "response") != 1L)
    stop("Each formula must have a response matrix on the left-hand side.", call. = FALSE)
  rhs <- attr(trm, "term.labels")
  out <- if (length(rhs)) reformulate(rhs, response = "gd_mat") else gd_mat ~ 1
  environment(out) <- environment()
  out
}

.terradish_cv_score <- function(theta, formula, data, response,
                                conductance_model, measurement_model,
                                nu, nonnegative, measurement_control)
{
  term_labels <- attr(terms(formula), "term.labels")
  model_formula <- if (length(term_labels)) reformulate(term_labels) else ~1
  f <- conductance_model(model_formula, data$x)
  args <- list(f = f, g = measurement_model, s = data, S = response,
               theta = theta, nu = nu, gradient = FALSE, hessian = FALSE,
               partial = FALSE, nonnegative = nonnegative, cores = 1L)
  if (!is.null(measurement_control))
    args$measurement_control <- measurement_control
  obj <- do.call(terradish_algorithm, args)
  list(loglik = -obj$objective, phi = obj$phi,
       convergence = obj$subproblem$convergence,
       iterations = obj$subproblem$iters)
}

#' Fixed-domain cross-validation for terradish models
#'
#' Fits one or more conductance formulas to each training fold while retaining
#' the complete landscape graph, then evaluates the trained coefficients on the
#' held-out focal sites. Nuisance parameters are reprofiled on every test fold.
#'
#' @param data A prebuilt \code{terradish_graph}, usually returned by
#'   \code{\link{conductance_surface}}.
#' @param formulas A formula or named list of formulas. All formulas must use
#'   the same response matrix and measurement model.
#' @param folds Fold membership for each focal site. Use
#'   \code{\link{terradish_folds}}, or supply study-specific labels.
#' @param model Measurement model function or one of \code{"mlpe"},
#'   \code{"wishart"}, or \code{"ls"}.
#' @param nu Effective Wishart degrees of freedom when required by the selected
#'   measurement model.
#' @param conductance_model Conductance-model factory.
#' @param baseline Should each test fold also be scored at zero conductance
#'   coefficients, the uniform-conductance baseline on the same graph?
#' @param keep_fits Whether to retain no training fits, slim fits, or full fits.
#' @param checkpoint Optional file path. Completed fold results are written
#'   atomically after every successful or failed fit so a long run can be
#'   inspected or resumed.
#' @param resume Logical. If \code{TRUE} and \code{checkpoint} exists, continue
#'   a compatible run without reevaluating completed formula-fold combinations.
#' @param cores Number of workers used inside each model fit. Fold evaluation is
#'   serial to avoid nested parallel clusters.
#' @param measurement_control Optional \code{\link{NewtonRaphsonControl}} object
#'   for profiling nuisance parameters.
#' @param ... Additional arguments passed to \code{\link{terradish}}.
#'
#' @details
#' The graph domain and raster covariates remain fixed. Only the focal-site
#' columns in \code{demes} and \code{rhs}, the response matrix, and any
#' measurement-model site attributes are subset. This evaluates transfer to
#' held-out sites within the observed landscape rather than rebuilding a new
#' landscape for every fold.
#'
#' \code{loglik_gain} is the trained model's held-out log-likelihood minus the
#' uniform-conductance baseline log-likelihood. Larger values indicate improved
#' prediction over that baseline. Raw log-likelihood and gain are comparable
#' only among models using the same response representation, measurement
#' likelihood, folds, graph, and effective degrees of freedom. They are not a
#' common scoring rule for comparing MLPE, generalized-Wishart distance, and
#' covariance-Wishart models.
#'
#' @return An object of class \code{"terradish_cv_folds"}. Its \code{folds}
#'   element records the site assignments, \code{results} contains one row per
#'   formula and fold, and \code{summary} reports summed and mean held-out scores
#'   by formula. Optional retained fits are stored in \code{fits}.
#'
#' @examples
#' \donttest{
#' data(melip)
#' altitude <- terra::unwrap(melip.altitude)
#' forestcover <- terra::unwrap(melip.forestcover)
#' coords <- terra::unwrap(melip.coords)
#' keep <- 1:8
#' response <- melip.Fst[keep, keep]
#' covariates <- terra::aggregate(c(altitude, forestcover), fact = 6,
#'                                na.rm = TRUE)
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, coords[keep])
#' coords_projected <- terra::project(coords[keep], "EPSG:5070")
#' folds <- terradish_folds(coords_projected, k = 2, seed = 1)
#' cv <- terradish_cv_folds(
#'   surface,
#'   list(altitude = response ~ altitude,
#'        both = response ~ altitude + forestcover),
#'   folds = folds,
#'   model = "ls",
#'   control = NewtonRaphsonControl(maxit = 1))
#' cv$summary
#' }
#'
#' @export
terradish_cv_folds <- function(data, formulas, folds, model = "mlpe", nu = NULL,
                               conductance_model = loglinear_conductance,
                               baseline = TRUE,
                               keep_fits = c("none", "slim", "full"),
                               checkpoint = NULL, resume = TRUE, cores = 1L,
                               measurement_control = NULL, ...)
{
  if (!inherits(data, c("terradish_graph", "radish_graph")))
    stop("`data` must be a terradish graph.", call. = FALSE)
  keep_fits <- match.arg(keep_fits)
  if (inherits(formulas, "formula"))
    formulas <- list(model = formulas)
  if (!is.list(formulas) || !length(formulas) ||
      !all(vapply(formulas, inherits, logical(1), "formula")))
    stop("`formulas` must be a formula or a list of formulas.", call. = FALSE)
  if (is.null(names(formulas)) || any(!nzchar(names(formulas))))
    names(formulas) <- paste0("model_", seq_along(formulas))
  if (anyDuplicated(names(formulas)))
    stop("Formula names must be unique.", call. = FALSE)
  response_labels <- vapply(formulas, function(x) {
    paste(deparse(attr(terms(x), "variables")[[2L]]), collapse = "")
  }, character(1))
  if (length(unique(response_labels)) != 1L)
    stop("All formulas must use the same response expression.", call. = FALSE)

  folds <- as.vector(folds)
  if (length(folds) != length(data$demes) || anyNA(folds))
    stop("`folds` must contain one non-missing label per focal site.", call. = FALSE)
  fold_levels <- unique(folds)
  if (length(fold_levels) < 2L || any(tabulate(match(folds, fold_levels)) < 2L))
    stop("At least two folds are required, with at least two focal sites per fold.",
         call. = FALSE)

  response_expr <- attr(terms(formulas[[1L]]), "variables")[[2L]]
  gd_full <- as.matrix(eval(response_expr, envir = parent.frame(),
                            enclos = environment(formulas[[1L]])))
  if (!identical(dim(gd_full), c(length(folds), length(folds))))
    stop("The response matrix dimensions must match the number of focal sites.",
         call. = FALSE)
  model_name <- .terradish_measurement_model_name(.terradish_measurement_model(model))
  if (is.null(model_name))
    model_name <- "custom"
  family <- .terradish_measurement_family(model_name)
  dots <- list(...)
  nonnegative <- if (is.null(dots$nonnegative)) TRUE else isTRUE(dots$nonnegative)
  if (!is.null(measurement_control))
    dots$measurement_control <- measurement_control

  formula_labels <- vapply(formulas, function(x) paste(deparse(x), collapse = ""),
                           character(1))
  checkpoint_meta <- list(
    folds = folds,
    formulas = formula_labels,
    measurement_model = model_name,
    nu = nu,
    baseline = isTRUE(baseline),
    keep_fits = keep_fits
  )
  rows <- list()
  fits <- list()
  completed <- character()
  if (!is.null(checkpoint) && isTRUE(resume) && file.exists(checkpoint))
  {
    saved <- readRDS(checkpoint)
    if (!is.list(saved) || is.null(saved$metadata) || is.null(saved$results) ||
        !identical(saved$metadata, checkpoint_meta))
      stop("The existing checkpoint is not compatible with this cross-validation run.",
           call. = FALSE)
    rows <- split(saved$results, seq_len(nrow(saved$results)))
    fits <- saved$fits
    if (is.null(fits))
      fits <- list()
    completed <- paste(saved$results$model, saved$results$fold, sep = "::")
  }
  row_i <- length(rows)
  for (model_i in seq_along(formulas))
  {
    cv_formula <- .terradish_cv_formula(formulas[[model_i]])
    for (fold_i in seq_along(fold_levels))
    {
      label <- fold_levels[[fold_i]]
      result_key <- paste(names(formulas)[model_i], label, sep = "::")
      if (result_key %in% completed)
        next
      test <- which(folds == label)
      train <- setdiff(seq_along(folds), test)
      train_data <- .terradish_subset_graph_focal(data, train)
      test_data <- .terradish_subset_graph_focal(data, test)
      train_model <- .terradish_measurement_model(model, subset = train)
      test_model <- .terradish_measurement_model(model, subset = test)
      gd_mat <- gd_full[train, train, drop = FALSE]
      started <- proc.time()[["elapsed"]]

      fit <- tryCatch(
        .fit_terradish_with_fallback(
          cv_formula, train_data, train_model, nu = nu,
          conductance_model = conductance_model, cores = cores, dots = dots,
          response_matrix = gd_mat),
        error = identity
      )
      if (!inherits(fit, "error") && length(fit$mle$theta_internal) < 1L)
        fit <- simpleError("The training fit produced no conductance coefficients.")
      row_i <- row_i + 1L
      if (inherits(fit, "error"))
      {
        rows[[row_i]] <- data.frame(
          model = names(formulas)[model_i], fold = as.character(label),
          n_train = length(train), n_test = length(test),
          loglik = NA_real_, baseline_loglik = NA_real_, loglik_gain = NA_real_,
          nuisance_iterations = NA_integer_, nuisance_convergence = NA_integer_,
          elapsed_seconds = proc.time()[["elapsed"]] - started,
          error = conditionMessage(fit), stringsAsFactors = FALSE)
      }
      else
      {
        score <- tryCatch(
          .terradish_cv_score(
            fit$mle$theta_internal, cv_formula, test_data,
            gd_full[test, test, drop = FALSE], conductance_model, test_model,
            nu, nonnegative = nonnegative, measurement_control = measurement_control),
          error = identity
        )
        base <- NULL
        if (!inherits(score, "error") && isTRUE(baseline))
          base <- tryCatch(
            .terradish_cv_score(
              rep(0, length(fit$mle$theta_internal)), cv_formula, test_data,
              gd_full[test, test, drop = FALSE], conductance_model, test_model,
              nu, nonnegative = nonnegative, measurement_control = measurement_control),
            error = identity)
        score_error <- if (inherits(score, "error")) conditionMessage(score) else
          if (inherits(base, "error")) conditionMessage(base) else ""
        ll <- if (inherits(score, "error")) NA_real_ else score$loglik
        base_ll <- if (is.null(base) || inherits(base, "error")) NA_real_ else base$loglik
        rows[[row_i]] <- data.frame(
          model = names(formulas)[model_i], fold = as.character(label),
          n_train = length(train), n_test = length(test),
          loglik = ll, baseline_loglik = base_ll,
          loglik_gain = if (is.finite(ll) && is.finite(base_ll)) ll - base_ll else NA_real_,
          nuisance_iterations = if (inherits(score, "error")) NA_integer_ else score$iterations,
          nuisance_convergence = if (inherits(score, "error")) NA_integer_ else score$convergence,
          elapsed_seconds = proc.time()[["elapsed"]] - started,
          error = score_error, stringsAsFactors = FALSE)
        if (!identical(keep_fits, "none"))
          fits[[result_key]] <-
            if (identical(keep_fits, "slim")) slim_terradish(fit) else fit
      }

      if (!is.null(checkpoint))
      {
        checkpoint_dir <- dirname(checkpoint)
        if (!dir.exists(checkpoint_dir))
          dir.create(checkpoint_dir, recursive = TRUE)
        tmp <- paste0(checkpoint, ".tmp")
        saveRDS(list(metadata = checkpoint_meta,
                     results = do.call(rbind, rows),
                     fits = fits), tmp)
        if (file.exists(checkpoint))
          unlink(checkpoint)
        if (!file.rename(tmp, checkpoint))
          stop("Could not replace the cross-validation checkpoint.", call. = FALSE)
      }
    }
  }

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  summary <- do.call(rbind, lapply(split(results, results$model), function(x) {
    sum_or_na <- function(value) if (all(is.na(value))) NA_real_ else sum(value, na.rm = TRUE)
    mean_or_na <- function(value) if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
    data.frame(model = x$model[1L], folds = nrow(x), failed = sum(nzchar(x$error)),
               total_loglik = sum_or_na(x$loglik),
               mean_loglik = mean_or_na(x$loglik),
               total_loglik_gain = sum_or_na(x$loglik_gain),
               mean_loglik_gain = mean_or_na(x$loglik_gain),
               stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL
  out <- list(folds = folds, results = results, summary = summary,
              measurement_model = model_name, likelihood_family = family,
              nu = nu, baseline = isTRUE(baseline))
  if (length(fits)) out$fits <- fits
  class(out) <- "terradish_cv_folds"
  out
}

#' @export
print.terradish_cv_folds <- function(x, ...)
{
  cat("Fixed-domain terradish cross-validation\n")
  cat("Measurement model:", x$measurement_model, "\n")
  cat("Folds:", length(unique(x$folds)), "\n\n")
  print(x$summary, row.names = FALSE, ...)
  invisible(x)
}
