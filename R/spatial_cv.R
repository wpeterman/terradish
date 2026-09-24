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
                                nu, nonnegative, measurement_control, phi = NULL)
{
  term_labels <- attr(terms(formula), "term.labels")
  model_formula <- if (length(term_labels)) reformulate(term_labels) else ~1
  f <- conductance_model(model_formula, data$x)
  fixed <- !is.null(phi)
  if (fixed) {
    original <- measurement_model
    fixed_phi <- c(phi)
    measurement_model <- function(E, S, phi, ...) {
      if (missing(phi)) return(list(phi = fixed_phi, lower = fixed_phi, upper = fixed_phi))
      original(E, S, phi = fixed_phi, ...)
    }
    attributes(measurement_model) <- attributes(original)
  }
  args <- list(f = f, g = measurement_model, s = data, S = response,
               theta = theta, nu = nu, gradient = FALSE, hessian = FALSE,
               partial = FALSE, nonnegative = nonnegative, cores = 1L)
  if (!is.null(measurement_control))
    args$measurement_control <- measurement_control
  obj <- do.call(terradish_algorithm, args)
  list(loglik = -obj$objective, phi = obj$phi,
       convergence = obj$subproblem$convergence,
       iterations = if (fixed) 0L else obj$subproblem$iters)
}

.cv_named_models <- function(value, names, label) {
  if (is.function(value) || is.character(value)) return(setNames(rep(list(value), length(names)), names))
  if (!is.list(value) || is.null(names(value)) || anyDuplicated(names(value)) ||
      !setequal(names(value), names))
    stop("`", label, "` must be one model or a named list matching `formulas`.", call. = FALSE)
  value[names]
}

.terradish_cv_prediction_theta <- function(fit, nuisance) {
  theta <- fit$mle$theta_internal
  if (length(theta)) return(theta)
  # With beta/tau fixed at zero, predictions do not depend on conductance.
  # A valid factory default supplies the otherwise unidentified coefficients.
  if (nuisance == "fixed" && .no_structure_boundary(fit$fit))
    return(attr(fit$submodels$f_internal, "default", exact = TRUE))
  stop("The training fit produced no conductance coefficients.", call. = FALSE)
}

.cv_prediction_metadata <- function(x) {
  # Fitted terms retain a formula environment for prediction. That transient
  # call frame is not part of a model definition or a checkpoint signature.
  attrs <- attributes(x)
  attrs$.Environment <- NULL
  if (is.list(x)) x <- lapply(x, .cv_prediction_metadata)
  if (length(attrs)) attrs <- lapply(attrs, .cv_prediction_metadata)
  attributes(x) <- attrs
  x
}

.cv_function_identity <- function(fn) {
  attrs <- attributes(fn)
  attrs <- attrs[!vapply(attrs, is.function, logical(1))]
  references <- unique(all.names(body(fn)))
  captured <- lapply(references, function(name) get0(name, envir = environment(fn), inherits = FALSE))
  names(captured) <- references
  captured <- captured[vapply(captured, function(x)
    !is.null(x) && (is.atomic(x) || is.data.frame(x)), logical(1))]
  list(formals = formals(fn), body = body(fn),
       attributes = .cv_prediction_metadata(attrs),
       captured = .cv_prediction_metadata(captured))
}

.cv_signature <- function(value) {
  path <- tempfile("terradish-cv-")
  on.exit(unlink(path))
  con <- file(path, "wb")
  tryCatch(serialize(value, con, version = 2), finally = close(con))
  unname(tools::md5sum(path))
}

#' Fixed-domain cross-validation for terradish models
#'
#' Fits one or more conductance formulas to each training fold while retaining
#' the complete landscape graph, then evaluates the trained coefficients on the
#' held-out focal sites. Measurement parameters can be reprofiled on each test
#' fold or held at their training estimates for predictive scoring.
#'
#' @param data A prebuilt \code{terradish_graph}, usually returned by
#'   \code{\link{conductance_surface}}.
#' @param formulas A formula or named list of formulas. All formulas must use
#'   the same response matrix and likelihood family.
#' @param folds Fold membership for each focal site. Use
#'   \code{\link{terradish_folds}}, or supply study-specific labels. A list of
#'   fold vectors requests repeated cross-validation.
#' @param model Measurement model function or one of \code{"mlpe"},
#'   \code{"wishart"}, or \code{"ls"}. With fixed nuisance parameters, a named
#'   list matching the formula names can compare models within one family.
#' @param nu Effective Wishart degrees of freedom when required by the selected
#'   measurement model.
#' @param conductance_model One conductance-model factory or a named list
#'   matching the formula names, allowing log-linear, Gaussian, and spline
#'   models in the same comparison.
#' @param nuisance \code{"reprofile"} (default) estimates measurement parameters
#'   separately on each held-out block. \code{"fixed"} holds every measurement
#'   parameter, including environmental coefficients, at its training estimate.
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
#' Reprofiled scores compare conductance formulas under one measurement model;
#' they must not be compared across calls with different measurement models or
#' effective degrees of freedom. Fixed-parameter scores are predictive marginal
#' densities and can compare measurement models within one likelihood family.
#' Neither mode scores pairs spanning training and test sites. Spatial folds
#' score only within-cluster pairs. Wishart score gains scale with \code{nu},
#' although ideal fitted rankings are invariant to its scale.
#'
#' Checkpoints include a signature of the response, graph, focal cells,
#' covariates, model definitions, formulas, \code{nu}, and fold assignments.
#' A mismatched signature stops resumption. Nonconverged fits and failed scores
#' count as failed folds. With fixed nuisance parameters, no-structure boundary
#' fits remain scoreable because predictions do not depend on conductance.
#' A failed uniform-baseline score is recorded separately: model scores remain
#' usable, but baseline gains are unavailable for that fold.
#' Models with failed folds receive
#' NA totals and an explicit success count. Rankings use only folds successful
#' for every model. Paired differences and their standard errors are descriptive;
#' repeated folds share observations and are not independent replicates.
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
#'   by formula, common-fold totals, paired differences from the best model,
#'   the standard error of their mean, and a one-standard-error flag.
#'   Repeated runs also report the SD of repeat totals over common folds.
#'   Optional retained fits are stored in \code{fits}.
#'
#' @examples
#' \donttest{
#' r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8,
#'                  ymin = 0, ymax = 8, crs = "EPSG:32617")
#' terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
#' names(r) <- "x"
#' sites <- terra::xyFromCell(r, c(1, 5, 8, 12, 20, 27, 33, 37, 45, 53, 57, 64))
#' surface <- conductance_surface(r, sites)
#' response <- simulate_covariance_response(c(x = 0.8), ~x, surface,
#'   nu = 1000, seed = 93)$covariance
#' folds <- terradish_folds(sites, k = 3, method = "random", seed = 1)
#' cv <- terradish_cv_folds(
#'   surface,
#'   list(uniform = response ~ 1, gradient = response ~ x),
#'   folds = folds,
#'   model = wishart_covariance, nu = 1000, nuisance = "fixed",
#'   control = NewtonRaphsonControl(maxit = 100))
#' cv$summary
#' # Check fold success counts before comparing common-fold predictive scores.
#' # Random folds illustrate the call; choose ecological spatial blocks in practice.
#' }
#'
#' @export
terradish_cv_folds <- function(data, formulas, folds, model = "mlpe", nu = NULL,
                               conductance_model = loglinear_conductance,
                               baseline = TRUE,
                               keep_fits = c("none", "slim", "full"),
                               checkpoint = NULL, resume = TRUE, cores = 1L,
                               measurement_control = NULL,
                               nuisance = c("reprofile", "fixed"), ...)
{
  if (!inherits(data, c("terradish_graph", "radish_graph")))
    stop("`data` must be a terradish graph.", call. = FALSE)
  keep_fits <- match.arg(keep_fits)
  nuisance <- match.arg(nuisance)
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

  repeats <- if (is.list(folds)) folds else list(folds)
  if (!length(repeats)) stop("`folds` must contain at least one repeat.", call. = FALSE)
  repeats <- lapply(repeats, as.vector)
  tasks <- do.call(rbind, lapply(seq_along(repeats), function(i) {
    values <- repeats[[i]]
    if (length(values) != length(data$demes) || anyNA(values))
      stop("`folds` must contain one non-missing label per focal site in each repeat.", call. = FALSE)
    levels <- unique(values)
    if (length(levels) < 2L || any(tabulate(match(values, levels)) < 2L))
      stop("At least two folds are required, with at least two focal sites per fold.", call. = FALSE)
    data.frame(repetition = i, fold = as.character(levels))
  }))

  response_expr <- attr(terms(formulas[[1L]]), "variables")[[2L]]
  gd_full <- as.matrix(eval(response_expr, envir = parent.frame(),
                            enclos = environment(formulas[[1L]])))
  if (!identical(dim(gd_full), rep(length(data$demes), 2)))
    stop("The response matrix dimensions must match the number of focal sites.",
         call. = FALSE)
  factories <- .cv_named_models(conductance_model, names(formulas), "conductance_model")
  models <- lapply(.cv_named_models(model, names(formulas), "model"), .terradish_measurement_model)
  if (nuisance == "reprofile" && !all(vapply(models, identical, logical(1), models[[1]])))
    stop("Reprofiled scores require one measurement model per call; use `nuisance = \"fixed\"` for measurement-model comparisons.", call. = FALSE)
  model_name <- vapply(models, function(g) {
    name <- .terradish_measurement_model_name(g)
    if (is.null(name)) "custom" else name
  }, character(1))
  families <- vapply(model_name, function(name) {
    family <- .terradish_measurement_family(name)
    if (is.null(family)) "unknown" else family
  }, character(1))
  if (length(unique(families)) != 1L || any(families == "unknown"))
    stop("Cross-validation models must share one known likelihood family and response.", call. = FALSE)
  family <- families[1]
  if (length(nu) > 1L) stop("One `nu` value is required per call.", call. = FALSE)
  dots <- list(...)
  nonnegative <- if (is.null(dots$nonnegative)) TRUE else isTRUE(dots$nonnegative)
  if (!is.null(measurement_control))
    dots$measurement_control <- measurement_control

  formula_labels <- vapply(formulas, function(x) paste(deparse(x), collapse = ""),
                           character(1))
  baseline_variable <- all.vars(delete.response(terms(formulas[[1L]])))[1]
  if (is.na(baseline_variable)) baseline_variable <- colnames(data$x)[1]
  baseline_formula <- reformulate(baseline_variable, response = "gd_mat")
  signature <- .cv_signature(list(response = gd_full,
    graph = c(vertices = nrow(data$x), edges = ncol(data$adj)), focal = data$demes,
    covariates = colnames(data$x), sums = colSums(data$x), formulas = formula_labels,
    conductance = lapply(seq_along(factories), function(i) {
      factory <- factories[[i]]
      rhs <- attr(terms(formulas[[i]]), "term.labels")
      fitted_model <- factory(if (length(rhs)) reformulate(rhs) else ~1, data$x)
      list(factory = .cv_function_identity(factory), model = .cv_function_identity(fitted_model))
    }), measurement = lapply(models, .cv_function_identity), nu = nu, folds = repeats))
  checkpoint_meta <- list(
    signature = signature,
    folds = repeats,
    formulas = formula_labels,
    measurement_model = model_name,
    nu = nu,
    baseline = isTRUE(baseline),
    keep_fits = keep_fits,
    nuisance = nuisance
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
    completed <- paste(saved$results$model, saved$results$repetition, saved$results$fold, sep = "::")
  }
  row_i <- length(rows)
  for (model_i in seq_along(formulas))
  {
    cv_formula <- .terradish_cv_formula(formulas[[model_i]])
    for (fold_i in seq_len(nrow(tasks)))
    {
      repetition <- tasks$repetition[fold_i]
      memberships <- repeats[[repetition]]
      label <- tasks$fold[fold_i]
      result_key <- paste(names(formulas)[model_i], repetition, label, sep = "::")
      if (result_key %in% completed)
        next
      test <- which(as.character(memberships) == label)
      train <- setdiff(seq_along(memberships), test)
      train_data <- .terradish_subset_graph_focal(data, train)
      test_data <- .terradish_subset_graph_focal(data, test)
      train_model <- .terradish_measurement_model(models[[model_i]], subset = train)
      test_model <- .terradish_measurement_model(models[[model_i]], subset = test)
      factory <- factories[[model_i]]
      gd_mat <- gd_full[train, train, drop = FALSE]
      started <- proc.time()[["elapsed"]]

      fit <- tryCatch(
        .fit_terradish_with_fallback(
          cv_formula, train_data, train_model, nu = nu,
          conductance_model = factory, cores = cores, dots = dots,
          response_matrix = gd_mat),
        error = identity
      )
      if (!inherits(fit, "error") && !isTRUE(fit$convergence$code == 0))
        fit <- simpleError("The training fit did not converge.")
      if (!inherits(fit, "error")) {
        prediction_theta <- tryCatch(.terradish_cv_prediction_theta(fit, nuisance),
                                     error = identity)
        if (inherits(prediction_theta, "error")) fit <- prediction_theta
      }
      row_i <- row_i + 1L
      if (inherits(fit, "error"))
      {
        rows[[row_i]] <- data.frame(
          model = names(formulas)[model_i], fold = as.character(label),
          repetition = repetition,
          n_train = length(train), n_test = length(test),
          loglik = NA_real_, baseline_loglik = NA_real_, loglik_gain = NA_real_,
          nuisance_iterations = NA_integer_, nuisance_convergence = NA_integer_,
          elapsed_seconds = proc.time()[["elapsed"]] - started,
          error = conditionMessage(fit), baseline_error = "", stringsAsFactors = FALSE)
      }
      else
      {
        score <- tryCatch(
          .terradish_cv_score(
            prediction_theta, cv_formula, test_data,
            gd_full[test, test, drop = FALSE], fit$submodels$f_factory, test_model,
            nu, nonnegative = nonnegative, measurement_control = measurement_control,
            phi = if (nuisance == "fixed") fit$fit$phi else NULL),
          error = identity
        )
        if (!inherits(score, "error") && !isTRUE(score$convergence == 0))
          score <- simpleError("Held-out nuisance profiling did not converge.")
        base <- NULL
        if (!inherits(score, "error") && isTRUE(baseline)) {
          base <- tryCatch({
            baseline_phi <- NULL
            if (nuisance == "fixed") {
              baseline_training <- .terradish_cv_score(0, baseline_formula, train_data, gd_mat,
                loglinear_conductance, train_model, nu, nonnegative,
                measurement_control)
              if (!isTRUE(baseline_training$convergence == 0))
                stop("Uniform-baseline training nuisance profiling did not converge.")
              baseline_phi <- baseline_training$phi
            }
            .terradish_cv_score(
              0, baseline_formula, test_data,
              gd_full[test, test, drop = FALSE], loglinear_conductance, test_model,
              nu, nonnegative = nonnegative, measurement_control = measurement_control,
              phi = baseline_phi)
            },
            error = identity)
          if (!inherits(base, "error") && !isTRUE(base$convergence == 0))
            base <- simpleError("Uniform-baseline held-out nuisance profiling did not converge.")
        }
        score_error <- if (inherits(score, "error")) conditionMessage(score) else ""
        baseline_error <- if (inherits(base, "error")) conditionMessage(base) else ""
        ll <- if (inherits(score, "error")) NA_real_ else score$loglik
        base_ll <- if (is.null(base) || inherits(base, "error")) NA_real_ else base$loglik
        rows[[row_i]] <- data.frame(
          model = names(formulas)[model_i], fold = as.character(label),
          repetition = repetition,
          n_train = length(train), n_test = length(test),
          loglik = ll, baseline_loglik = base_ll,
          loglik_gain = if (is.finite(ll) && is.finite(base_ll)) ll - base_ll else NA_real_,
          nuisance_iterations = if (inherits(score, "error")) NA_integer_ else score$iterations,
          nuisance_convergence = if (inherits(score, "error")) NA_integer_ else score$convergence,
          elapsed_seconds = proc.time()[["elapsed"]] - started,
          error = score_error, baseline_error = baseline_error, stringsAsFactors = FALSE)
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
  if (any(nzchar(results$baseline_error))) {
    bad <- nzchar(results$baseline_error)
    warning("Failed uniform-baseline scores: ", paste(paste(results$model[bad],
      results$repetition[bad], results$fold[bad], sep = "/"), collapse = ", "),
      ". Model scores remain available; affected gains are NA.", call. = FALSE)
  }
  key <- paste(results$repetition, results$fold, sep = "::")
  successful <- is.finite(results$loglik) & !nzchar(results$error)
  common_keys <- names(which(vapply(split(successful, key), function(x)
    length(x) == length(formulas) && all(x), logical(1))))
  results$common_fold <- key %in% common_keys
  summary <- do.call(rbind, lapply(split(results, results$model), function(x) {
    ok <- is.finite(x$loglik) & !nzchar(x$error)
    sum_or_na <- function(value) if (!all(ok) || any(!is.finite(value))) NA_real_ else sum(value)
    mean_or_na <- function(value) if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
    repeated <- tapply(x$loglik[x$common_fold], x$repetition[x$common_fold], sum)
    data.frame(model = x$model[1L], folds = nrow(x), failed = sum(!ok), n_folds_ok = sum(ok),
               total_loglik = sum_or_na(x$loglik),
               common_folds_total = if (any(x$common_fold)) sum(x$loglik[x$common_fold]) else NA_real_,
               n_common_folds = sum(x$common_fold),
               between_repeat_sd = if (length(repeated) > 1L) sd(repeated) else NA_real_,
               mean_loglik = mean_or_na(x$loglik),
               total_loglik_gain = sum_or_na(x$loglik_gain),
               mean_loglik_gain = mean_or_na(x$loglik_gain),
               stringsAsFactors = FALSE)
  }))
  summary <- summary[order(-summary$common_folds_total, summary$model, na.last = TRUE), , drop = FALSE]
  best <- summary$model[1L]
  best_rows <- results[results$model == best, ]
  best_scores <- setNames(best_rows$loglik, paste(best_rows$repetition, best_rows$fold, sep = "::"))
  results$paired_difference <- ifelse(results$common_fold,
    results$loglik - best_scores[key], NA_real_)
  summary$mean_difference <- summary$se_difference <- NA_real_
  summary$within_one_se <- FALSE
  for (i in seq_len(nrow(summary))) {
    differences <- results$paired_difference[results$model == summary$model[i] & results$common_fold]
    if (length(differences)) summary$mean_difference[i] <- mean(differences)
    if (length(differences) > 1L) {
      summary$se_difference[i] <- sd(differences) / sqrt(length(differences))
      summary$within_one_se[i] <- abs(mean(differences)) <= summary$se_difference[i]
    }
  }
  if (any(!successful))
    warning("Failed cross-validation folds: ", paste(paste(results$model[!successful],
      results$repetition[!successful], results$fold[!successful], sep = "/"), collapse = ", "),
      ". Totals are NA; ranking uses folds common to every model.", call. = FALSE)
  rownames(summary) <- NULL
  out <- list(folds = folds, results = results, summary = summary,
              measurement_model = model_name, likelihood_family = family,
              nu = nu, baseline = isTRUE(baseline), nuisance = nuisance,
              repeats = length(repeats), checkpoint_signature = signature)
  if (length(fits)) out$fits <- fits
  class(out) <- "terradish_cv_folds"
  out
}

#' @export
print.terradish_cv_folds <- function(x, ...)
{
  cat("Fixed-domain terradish cross-validation\n")
  cat("Measurement model:", x$measurement_model, "\n")
  cat("Repeats:", x$repeats, "| nuisance:", x$nuisance, "\n")
  cat("Ranks use common successful folds; differences are paired against the best model.\n\n")
  print(x$summary, row.names = FALSE, ...)
  invisible(x)
}
