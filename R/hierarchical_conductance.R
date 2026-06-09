# =====================================================================
# Tier 2: hierarchical conductance surface  log c = X theta + Z_u u
# A spatially smooth Gaussian residual field u is added to the covariate
# log-conductance, so the model is mechanistic where covariates explain the
# data and flexible where they do not. Nests terradish (tau^2 -> 0) and the
# FEEMS-style free surface (no covariates). Fit by penalized L-BFGS over
# (theta, u), gradient-only, reusing terradish_algorithm's adjoint, which
# avoids forming the dense Hessian over the field. The field variance tau^2 is
# selected by a Laplace (empirical-Bayes) marginal likelihood.
# =====================================================================

# Build a loglinear conductance model directly from a combined design matrix
# D = [X | Z_u]. conductance = exp(D %*% par); par = c(theta, u).
.design_loglinear_model <- function(D)
{
  D <- as.matrix(D)
  cm <- function(theta)
  {
    cond <- as.vector(exp(D %*% theta))
    cond <- .validate_conductance_values(cond, context = "hierarchical conductance")
    dfm  <- cond * D
    list(conductance        = cond,
         confint            = function(...) NULL,
         df__dx             = function(k) cond * theta[k],
         df__dtheta         = function(k) dfm[, k],
         df__dtheta_matrix  = dfm,
         d2f__dtheta_dtheta = function(k, l) cond * D[, k] * D[, l],
         d2f__dtheta_dx     = function(k, l) cond * ((k == l) + D[, k] * theta[l]))
  }
  class(cm) <- c("terradish_conductance_model", "radish_conductance_model")
  attr(cm, "default") <- rep(0, ncol(D))
  cm
}

# Build the coarse piecewise-constant field basis Z (n_active x m) and the
# GMRF precision Q = L_coarse + eps I on a G x G grid over the graph extent.
.build_conductance_field <- function(coords, G = 6L, eps = 1e-3)
{
  coords <- as.matrix(coords)
  G <- as.integer(G)
  if (length(G) != 1L || is.na(G) || G < 2L)
    stop("`field_resolution` (G) must be an integer >= 2.", call. = FALSE)
  xr <- range(coords[, 1]); yr <- range(coords[, 2])
  cx <- pmin(G, 1L + floor(G * (coords[, 1] - xr[1]) / (diff(xr) + 1e-9)))
  cy <- pmin(G, 1L + floor(G * (coords[, 2] - yr[1]) / (diff(yr) + 1e-9)))
  cell <- (cy - 1L) * G + cx
  occ  <- sort(unique(cell))
  m    <- length(occ)
  if (m < 2L)
    stop("The conductance field has fewer than two occupied coarse cells; ",
         "increase `field_resolution`.", call. = FALSE)
  remap <- match(cell, occ)
  Z <- matrix(0, nrow(coords), m)
  Z[cbind(seq_len(nrow(coords)), remap)] <- 1
  gx <- ((occ - 1L) %% G) + 1L
  gy <- ((occ - 1L) %/% G) + 1L
  A <- matrix(0, m, m)
  for (a in seq_len(m))
    for (b in seq_len(m))
      if (a != b && abs(gx[a] - gx[b]) + abs(gy[a] - gy[b]) == 1L)
        A[a, b] <- 1
  Q <- diag(rowSums(A)) - A + eps * diag(m)
  list(Z = Z, Q = Q, m = m, G = G, occ = occ,
       coarse_xy = cbind(gx, gy), eps = eps)
}

# Penalized objective (+ optional gradient) over par = c(theta, u) at fixed tau2.
# phi_state is an environment holding a warm-started nuisance vector.
.hierarchical_penalized <- function(par, tau2, cm, measurement_model, data, S,
                                    nu, uidx, Q, nonnegative, solver,
                                    solver_control, phi_state, want_grad,
                                    want_hess = FALSE)
{
  res <- terradish_algorithm(f = cm, g = measurement_model, s = data, S = S,
                             theta = par, nu = nu,
                             phi = if (is.null(phi_state)) NULL else phi_state$value,
                             gradient = want_grad, hessian = want_hess,
                             partial = FALSE, nonnegative = nonnegative,
                             solver = solver, solver_control = solver_control)
  if (!is.null(phi_state)) phi_state$value <- res$phi
  u <- par[uidx]
  Qu <- as.vector(Q %*% u)
  pen <- as.numeric(crossprod(u, Qu)) / (2 * tau2)
  out <- list(objective = res$objective + pen,
              loglik_objective = res$objective,
              phi = res$phi, boundary = res$boundary)
  if (want_grad) {
    g <- res$gradient
    g[uidx] <- g[uidx] + Qu / tau2
    out$gradient <- g
  }
  if (want_hess) out$hessian <- res$hessian   # likelihood Hessian (penalty added by caller)
  out
}

# Inner fit: minimize the penalized objective over (theta, u) at fixed tau2.
.hierarchical_fit_fixed <- function(par0, tau2, cm, measurement_model, data, S,
                                    nu, uidx, Q, nonnegative, solver,
                                    solver_control, maxit, factr)
{
  phi_state <- new.env(parent = emptyenv())
  phi_state$value <- NULL
  fn <- function(p)
    .hierarchical_penalized(p, tau2, cm, measurement_model, data, S, nu, uidx,
                            Q, nonnegative, solver, solver_control, phi_state,
                            want_grad = FALSE)$objective
  gr <- function(p)
    .hierarchical_penalized(p, tau2, cm, measurement_model, data, S, nu, uidx,
                            Q, nonnegative, solver, solver_control, phi_state,
                            want_grad = TRUE)$gradient
  opt <- optim(par0, fn, gr, method = "L-BFGS-B",
               control = list(maxit = maxit, factr = factr))
  list(par = opt$par, value = opt$value, convergence = opt$convergence,
       phi = phi_state$value)
}

# Laplace (empirical-Bayes) marginal log-likelihood of tau2 at the fitted optimum.
.hierarchical_logml <- function(fit, tau2, cm, measurement_model, data, S, nu,
                                uidx, Q, nonnegative, solver, solver_control)
{
  m <- length(uidx)
  u <- fit$par[uidx]
  res <- .hierarchical_penalized(fit$par, tau2, cm, measurement_model, data, S,
                                 nu, uidx, Q, nonnegative, solver,
                                 solver_control, phi_state = NULL,
                                 want_grad = FALSE, want_hess = TRUE)
  Huu <- res$hessian[uidx, uidx, drop = FALSE] + Q / tau2
  ld_H <- tryCatch(2 * sum(log(diag(chol((Huu + t(Huu)) / 2)))),
                   error = function(e) {
                     ev <- eigen(Huu, symmetric = TRUE, only.values = TRUE)$values
                     sum(log(pmax(ev, .Machine$double.eps)))
                   })
  ld_Q <- determinant(Q, logarithm = TRUE)$modulus
  nll  <- res$loglik_objective
  Qu   <- as.vector(Q %*% u)
  pen  <- as.numeric(crossprod(u, Qu)) / (2 * tau2)
  as.numeric(-nll - (m / 2) * log(tau2) + 0.5 * ld_Q - pen - 0.5 * ld_H)
}

#' Hierarchical (covariate + smooth field) conductance surface
#'
#' Fits a conductance surface whose log-conductance is the sum of a mechanistic
#' covariate term and a spatially smooth Gaussian random field,
#' \eqn{\log c = X\theta + u}.  This nests the standard \code{\link{terradish}}
#' model (recovered as the field variance \eqn{\tau^2 \to 0}) and a
#' FEEMS-style free surface (recovered when no covariates are supplied), and it
#' protects the covariate effects \eqn{\theta} from omitted-variable bias by
#' letting the field \eqn{u} absorb spatial structure the covariates do not
#' explain.  The fitted field is itself a deliverable: a map of conductance
#' variation your covariates failed to capture.
#'
#' @param formula Model formula with the genetic response matrix on the
#'   left-hand side and conductance covariates on the right, exactly as in
#'   \code{\link{terradish}}.  With no right-hand-side covariates
#'   (\code{response ~ 1}) the model is a pure smooth surface.
#' @param data A \code{terradish_graph} from \code{\link{conductance_surface}}.
#' @param conductance_model Conductance-model factory for the covariate term.
#'   Only log-link models are supported (the field is additive on the log scale);
#'   the default \code{\link{loglinear_conductance}} is recommended.
#' @param measurement_model A \code{terradish} measurement model, e.g.
#'   \code{\link{wishart_covariance}} (default), \code{\link{generalized_wishart}},
#'   \code{\link{mlpe}}, or \code{\link{leastsquares}}.
#' @param nu Effective Wishart degrees of freedom, required by the Wishart
#'   measurement models.
#' @param field_resolution Integer \eqn{G}: the field lives on a coarse
#'   \eqn{G \times G} grid over the graph extent (piecewise-constant within each
#'   occupied coarse cell).  A coarser field than the covariate grid is the main
#'   identifiability safeguard.  Default \code{6}.
#' @param tau2 Field variance.  Either a positive number (fit at that fixed
#'   value) or \code{"reml"} (default), which selects \eqn{\tau^2} by maximizing
#'   a Laplace empirical-Bayes marginal likelihood over \code{tau2_grid}.
#' @param tau2_grid Candidate \eqn{\tau^2} values used when \code{tau2 = "reml"}.
#'   Defaults to a log-spaced grid spanning \code{1e-2} to \code{1e2}.
#' @param eps Ridge added to the field GMRF precision (\eqn{L_{coarse} + \epsilon I})
#'   so the prior is proper.  Default \code{1e-3}.
#' @param theta Optional starting values for the covariate parameters.
#' @param nonnegative Force nonnegative measurement-model slope where applicable.
#' @param solver,solver_control Passed to \code{\link{terradish_algorithm}}.
#' @param maxit,factr \code{L-BFGS-B} controls for the inner penalized fit.
#' @param verbose Print progress during \eqn{\tau^2} selection?
#'
#' @details
#' The penalized objective is
#' \deqn{J(\theta, u) = -\ell(\theta, u \mid \phi) + \frac{1}{2\tau^2} u^\top (L_{coarse} + \epsilon I)\, u,}
#' where \eqn{\ell} is the measurement-model log-likelihood (with nuisance
#' parameters \eqn{\phi} profiled out) and \eqn{L_{coarse}} is the Laplacian of
#' the coarse-field adjacency graph, so the penalty is exactly the FEEMS spatial
#' smoothness prior written as a proper GMRF.  Optimization is by \code{L-BFGS-B}
#' over \eqn{(\theta, u)} using only the gradient, which flows through the same
#' single sparse-Cholesky solve that \code{\link{terradish_algorithm}} already
#' uses; the dense Hessian over the field is never formed during fitting.
#'
#' \eqn{\tau^2} is chosen (when \code{tau2 = "reml"}) by a Laplace approximation
#' to the marginal likelihood that integrates out the field, evaluated at each
#' candidate value's penalized optimum.
#'
#' @return An object of class \code{"terradish_hierarchical"} containing the
#'   covariate estimates \code{theta} (with conditional standard errors), the
#'   fitted field \code{u}, the selected \code{tau2}, the profiled nuisance
#'   parameters \code{phi}, the maximized log-likelihood and marginal
#'   log-likelihood, and the field construction needed by
#'   \code{\link{conductance_field}}.
#'
#' @seealso \code{\link{terradish}}, \code{\link{conductance_field}},
#'   \code{\link{conductance_surface}}, \code{\link{wishart_covariance}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data(melip)
#' covs <- c(terra::scale(terra::unwrap(melip.altitude)),
#'           terra::scale(terra::unwrap(melip.forestcover)))
#' names(covs) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covs, terra::unwrap(melip.coords),
#'                                directions = 8)
#'
#' fit <- terradish_hierarchical(
#'   melip.Fst ~ altitude + forestcover,
#'   data              = surface,
#'   measurement_model = generalized_wishart,
#'   nu                = 1000,
#'   field_resolution  = 6
#' )
#' summary(fit)
#' plot(conductance_field(fit, surface))   # map of unexplained conductance
#' }
#'
#' @export
terradish_hierarchical <- function(formula, data,
                                   conductance_model = loglinear_conductance,
                                   measurement_model = wishart_covariance,
                                   nu = NULL,
                                   field_resolution = 6L,
                                   tau2 = "reml",
                                   tau2_grid = NULL,
                                   eps = 1e-3,
                                   theta = NULL,
                                   nonnegative = TRUE,
                                   solver = c("direct", "auto", "amg", "pcg", "pcg_jacobi"),
                                   solver_control = NULL,
                                   maxit = 500L,
                                   factr = 1e7,
                                   verbose = TRUE)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(inherits(measurement_model, c("terradish_measurement_model",
                                          "radish_measurement_model")))
  solver <- match.arg(solver)

  # response + RHS formula (mirrors terradish())
  tm       <- terms(formula)
  vars     <- as.character(attr(tm, "variables"))[-1]
  response <- attr(tm, "response")
  if (!response)
    stop("'formula' must have the genetic response matrix on the left-hand side")
  S        <- eval(attr(tm, "variables")[[response + 1L]], parent.frame())
  S        <- as.matrix(S)
  is_ibd   <- length(vars) == 1
  rhs      <- if (!is_ibd) reformulate(attr(tm, "term.labels")) else formula(~1)

  stopifnot(nrow(S) == ncol(S), length(data$demes) == nrow(S))

  # covariate design X (no intercept; same convention as loglinear_conductance)
  X <- if (!is_ibd) assemble_model_matrix(rhs, data$x) else
    matrix(numeric(0), nrow = nrow(data$x), ncol = 0L)
  p <- ncol(X)
  theta_names <- if (p) colnames(X) else character(0)

  # coarse smooth field
  field <- .build_conductance_field(data$vertex_coordinates,
                                    G = field_resolution, eps = eps)
  m <- field$m
  D <- cbind(X, field$Z)
  cm <- .design_loglinear_model(D)
  uidx <- (p + 1L):(p + m)

  # starting values
  par0 <- rep(0, p + m)
  if (!is.null(theta)) {
    stopifnot(length(theta) == p)
    par0[seq_len(p)] <- theta
  }

  fit_at <- function(t2, par_start)
    .hierarchical_fit_fixed(par_start, t2, cm, measurement_model, data, S, nu,
                            uidx, field$Q, nonnegative, solver, solver_control,
                            maxit, factr)

  # ---- select tau2 ----
  tau2_selection <- NULL
  if (is.character(tau2) && identical(tau2, "reml")) {
    if (is.null(tau2_grid))
      tau2_grid <- 10 ^ seq(-2, 2, length.out = 9)
    logml <- rep(NA_real_, length(tau2_grid))
    par_start <- par0
    for (i in seq_along(tau2_grid)) {
      fi <- fit_at(tau2_grid[i], par_start)
      par_start <- fi$par                       # warm start across the grid
      logml[i] <- .hierarchical_logml(fi, tau2_grid[i], cm, measurement_model,
                                      data, S, nu, uidx, field$Q, nonnegative,
                                      solver, solver_control)
      if (verbose)
        message(sprintf("  tau2 = %10.4g   logML = %.3f", tau2_grid[i], logml[i]))
    }
    tau2_hat <- tau2_grid[which.max(logml)]
    tau2_selection <- data.frame(tau2 = tau2_grid, logML = logml)
  } else {
    stopifnot(is.numeric(tau2), length(tau2) == 1L, tau2 > 0)
    tau2_hat <- tau2
  }

  # ---- final fit + conditional theta covariance ----
  fit <- fit_at(tau2_hat, par0)
  hres <- .hierarchical_penalized(fit$par, tau2_hat, cm, measurement_model, data,
                                  S, nu, uidx, field$Q, nonnegative, solver,
                                  solver_control, phi_state = NULL,
                                  want_grad = FALSE, want_hess = TRUE)
  H_pen <- hres$hessian
  H_pen[uidx, uidx] <- H_pen[uidx, uidx] + field$Q / tau2_hat
  vcov_full <- .safe_invert((H_pen + t(H_pen)) / 2)
  theta_hat <- fit$par[seq_len(p)]
  u_hat <- fit$par[uidx]
  names(theta_hat) <- theta_names
  theta_se <- if (p) sqrt(pmax(diag(vcov_full)[seq_len(p)], 0)) else numeric(0)
  names(theta_se) <- theta_names

  loglik <- -hres$loglik_objective
  logml_final <- .hierarchical_logml(fit, tau2_hat, cm, measurement_model, data,
                                     S, nu, uidx, field$Q, nonnegative, solver,
                                     solver_control)

  # Effective degrees of freedom of the (penalized) field: trace of the
  # smoother "hat" on the u-block, edf = tr[(H_lik + Q/tau2)^-1 H_lik]. Total df
  # counts covariates + nuisance phi + tau2 + the field's effective df (so AIC
  # treats the smooth field honestly, as in a GAM / mixed model).
  Huu_lik <- hres$hessian[uidx, uidx, drop = FALSE]
  edf_field <- tryCatch(
    sum(diag(solve(H_pen[uidx, uidx, drop = FALSE], Huu_lik))),
    error = function(e) NA_real_)
  n_phi <- length(fit$phi)
  df <- p + n_phi + 1 + edf_field          # +1 for tau2
  vcov_theta <- if (p) vcov_full[seq_len(p), seq_len(p), drop = FALSE] else
    matrix(numeric(0), 0, 0)
  if (p) dimnames(vcov_theta) <- list(theta_names, theta_names)

  out <- list(
    call = match.call(),
    formula = formula,
    theta = theta_hat,
    theta_se = theta_se,
    u = u_hat,
    field = field,
    tau2 = tau2_hat,
    tau2_selection = tau2_selection,
    phi = fit$phi,
    loglik = loglik,
    logML = logml_final,
    df = df,
    edf_field = edf_field,
    aic = -2 * loglik + 2 * df,
    vcov = vcov_theta,
    vcov_full = vcov_full,
    convergence = fit$convergence,
    npar_covariate = p,
    npar_field = m,
    measurement_model = .terradish_measurement_model_name(measurement_model),
    nu = nu)
  class(out) <- c("terradish_hierarchical", "terradish")
  out
}

#' Extract the fitted conductance field from a hierarchical model
#'
#' Returns the smooth residual field \eqn{u} (the conductance variation the
#' covariates did not explain) from a \code{\link{terradish_hierarchical}} fit,
#' as a \code{terra::SpatRaster} when the graph stored its raster stack, or as a
#' per-cell numeric vector otherwise.
#'
#' @param fit A \code{"terradish_hierarchical"} object.
#' @param data The \code{terradish_graph} used to fit the model.
#' @param type \code{"field"} (default) returns the smooth field \eqn{u};
#'   \code{"logconductance"} returns the full fitted \eqn{\log c = X\theta + u};
#'   \code{"conductance"} returns \eqn{\exp(\log c)}.
#'
#' @return A single-layer \code{SpatRaster} (if \code{data$stack} is available)
#'   or a numeric vector over active cells.
#'
#' @seealso \code{\link{terradish_hierarchical}}
#' @export
conductance_field <- function(fit, data, type = c("field", "logconductance", "conductance"))
{
  stopifnot(inherits(fit, "terradish_hierarchical"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  type <- match.arg(type)

  u_cells <- as.vector(fit$field$Z %*% fit$u)
  values_out <- u_cells
  if (type != "field") {
    p <- fit$npar_covariate
    Xpart <- if (p) {
      rhs <- if (length(all.vars(fit$formula)) > 1)
        reformulate(attr(terms(fit$formula), "term.labels")) else formula(~1)
      X <- assemble_model_matrix(rhs, data$x)
      as.vector(X %*% fit$theta)
    } else rep(0, length(u_cells))
    logc <- Xpart + u_cells
    values_out <- if (identical(type, "conductance")) exp(logc) else logc
  }

  if (!is.null(data$stack)) {
    template <- data$stack[[1]]
    tv <- terra::values(template, dataframe = FALSE)[, 1]
    missing <- is.na(tv)
    tv[!missing] <- values_out
    out <- terra::setValues(template, tv)
    names(out) <- type
    return(out)
  }
  warning("No raster stored in `data`; returning a per-cell vector.")
  values_out
}

#' @export
print.terradish_hierarchical <- function(x, ...)
{
  cat("Hierarchical conductance surface  log c = X theta + u\n")
  cat(sprintf("  measurement model: %s%s\n", x$measurement_model,
              if (is.null(x$nu)) "" else sprintf("  (nu = %g)", x$nu)))
  cat(sprintf("  covariate params: %d   field cells: %d (%dx%d coarse grid)\n",
              x$npar_covariate, x$npar_field, x$field$G, x$field$G))
  cat(sprintf("  tau^2 = %.4g   field sd = %.3f   field edf = %.2f\n",
              x$tau2, sd(x$u), x$edf_field))
  cat(sprintf("  loglik = %.3f   marginal loglik = %.3f\n", x$loglik, x$logML))
  cat(sprintf("  df = %.2f   AIC = %.2f\n", x$df, x$aic))
  if (x$npar_covariate) {
    cat("\nConductance coefficients:\n")
    z <- x$theta / x$theta_se
    tab <- cbind(Estimate = x$theta, `Std. Error` = x$theta_se,
                 `z value` = z, `Pr(>|z|)` = pmin(2 * (1 - pnorm(abs(z))), 1))
    printCoefmat(tab, digits = 4, signif.stars = getOption("show.signif.stars"))
  }
  invisible(x)
}

#' @export
logLik.terradish_hierarchical <- function(object, ...) {
  val <- object$loglik; attr(val, "df") <- object$df; class(val) <- "logLik"; val
}

#' @export
AIC.terradish_hierarchical <- function(object, ..., k = 2) {
  if (identical(k, 2)) return(object$aic)
  -2 * object$loglik + k * object$df
}

#' @export
vcov.terradish_hierarchical <- function(object, ...) object$vcov

#' @export
confint.terradish_hierarchical <- function(object, parm, level = 0.95, ...) {
  if (!object$npar_covariate)
    return(matrix(numeric(0), 0, 2))
  z <- qnorm(1 - (1 - level) / 2)
  ci <- cbind(object$theta - z * object$theta_se, object$theta + z * object$theta_se)
  colnames(ci) <- paste0(round(100 * c((1 - level) / 2, 1 - (1 - level) / 2), 1), "%")
  rownames(ci) <- names(object$theta)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
plot.terradish_hierarchical <- function(x, data, type = c("field", "conductance", "logconductance"), ...) {
  type <- match.arg(type)
  r <- conductance_field(x, data, type = type)
  if (inherits(r, "SpatRaster"))
    terra::plot(r, main = paste0("Hierarchical conductance: ", type), ...)
  else
    stop("No raster stack in `data`; conductance_field() returned a vector.", call. = FALSE)
  invisible(r)
}

#' @export
summary.terradish_hierarchical <- function(object, ...)
{
  print(object, ...)
  if (!is.null(object$tau2_selection)) {
    cat("\ntau^2 selection (Laplace marginal likelihood):\n")
    print(round(object$tau2_selection, 3))
  }
  invisible(object)
}

#' @export
coef.terradish_hierarchical <- function(object, ...) object$theta
