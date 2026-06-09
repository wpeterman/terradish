# =====================================================================
# Tier 3: directional, non-reversible conductance (directed generator engine).
#
# log G_{a->b} = (eta_a + eta_b)/2 + d_{ab}' gamma_dir,  eta_k = X_k theta,
# where d_{ab} = -d_{ba} are directional edge covariates (elevation drop, flow,
# wind). gamma_dir = 0 gives a reversible generator (symmetric rates) -> the
# usual isolation-by-resistance model. The model distance is the symmetric
# COMMUTE time R^dir_ij = h(i->j) + h(j->i) of the directed Markov chain, so the
# (symmetric) genetic response and all existing measurement models are used
# unchanged. Fitting is by penalized/profiled likelihood with a reverse-mode
# (transpose-solve) adjoint; see dev/TIER3_DESIGN.md.
#
# Engine note: this R implementation factorizes one sub-generator per focal
# absorber per evaluation (forward + transpose solve). Correct and adequate for
# moderate graphs; an Eigen::SparseLU C++ backend (one factorization reused for
# forward + transpose) is the planned optimization for large rasters.
# =====================================================================

#' Directional edge covariates from a spatial layer
#'
#' Builds the antisymmetric per-directed-edge covariate
#' \eqn{d_{ab} = x_a - x_b} used by \code{\link{terradish_directed}} to drive
#' directional gene flow (e.g. elevation drop, flow accumulation, wind
#' potential).  \eqn{d_{ab} = -d_{ba}}, so a positive coefficient makes movement
#' from high to low \code{x} faster.
#'
#' @param x A \code{terra::SpatRaster} layer (or a numeric vector with one value
#'   per active graph cell) giving the potential whose gradient drives direction.
#' @param data A \code{terradish_graph} from \code{\link{conductance_surface}}.
#' @return A list with \code{edges} (an integer matrix of directed edges, columns
#'   \code{a}, \code{b}, 1-based node indices) and \code{d} (the matching numeric
#'   vector \eqn{d_{ab} = x_a - x_b}).
#' @seealso \code{\link{terradish_directed}}
#' @export
edge_gradient <- function(x, data)
{
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  if (inherits(x, "PackedSpatRaster")) x <- terra::unwrap(x)
  if (inherits(x, "SpatRaster")) {
    vals <- terra::values(x, dataframe = FALSE)[, 1]
    # active cells in vertex order: extract at vertex coordinates
    vals <- terra::extract(x, data$vertex_coordinates)[, ncol(terra::extract(x, data$vertex_coordinates))]
  } else {
    vals <- as.numeric(x)
  }
  n <- nrow(data$vertex_coordinates)
  if (length(vals) != n)
    stop("`x` must give one value per active graph cell (", n, ").", call. = FALSE)
  ed <- .directed_edges(data)
  list(edges = ed, d = vals[ed[, 1]] - vals[ed[, 2]])
}

# Directed edge list (both directions) from the graph's undirected edge_pairs.
.directed_edges <- function(data)
{
  ep <- data$edge_pairs
  if (is.null(ep)) {
    ep <- t(data$adj) + 1L                      # adj is 0-based upper-tri
  }
  ep <- as.matrix(ep); storage.mode(ep) <- "integer"
  rbind(cbind(a = ep[, 1], b = ep[, 2]),
        cbind(a = ep[, 2], b = ep[, 1]))
}

# Generator model: par = c(theta (p), gamma (q)) -> directed rates + design.
.directed_generator <- function(formula, data, directional)
{
  # Covariates absent from the symmetric formula may still drive direction
  # (passed via `directional`), so muffle the "unused covariate" warning here.
  X <- withCallingHandlers(
    assemble_model_matrix(formula, data$x),     # n_cells x p (no intercept)
    warning = function(w)
      if (grepl("Removed unused spatial covariates", conditionMessage(w)))
        invokeRestart("muffleWarning"))
  p <- ncol(X)
  ed <- .directed_edges(data)
  if (is.null(directional$edges) || nrow(directional$edges) != nrow(ed))
    stop("`directional` edges do not match the graph's directed edges; build it with edge_gradient(., data).",
         call. = FALSE)
  D <- as.matrix(directional$d)
  q <- ncol(D)
  sab <- (X[ed[, 1], , drop = FALSE] + X[ed[, 2], , drop = FALSE]) / 2  # n_edges x p
  n <- nrow(X)

  list(
    n = n, p = p, q = q, edges = ed, sab = sab, D = D,
    theta_names = colnames(X),
    gamma_names = if (is.null(colnames(D))) paste0("gamma", seq_len(q)) else
      paste0("gamma_", colnames(D)),
    rates = function(par) {
      theta <- par[seq_len(p)]; gamma <- par[p + seq_len(q)]
      as.numeric(exp(sab %*% theta + D %*% gamma))
    })
}

#' Likelihood of a directed (non-reversible) conductance surface
#'
#' Forward map and reverse-mode gradient for the directional generator engine.
#' Builds the directed Markov generator from edge rates, computes the symmetric
#' commute-time covariance \code{E} among focal sites, evaluates a measurement
#' model, and returns the profiled objective and its gradient with respect to
#' \code{(theta, gamma)} via a transpose-solve adjoint.  Users normally call
#' \code{\link{terradish_directed}} rather than this function.
#'
#' @param gen A directed generator model (internal; from
#'   \code{terradish_directed}).
#' @param g A \code{terradish} measurement model.
#' @param data A \code{terradish_graph}.
#' @param S Observed genetic response matrix.
#' @param par Numeric vector \code{c(theta, gamma)}.
#' @param nu Effective Wishart degrees of freedom (for Wishart models).
#' @param gradient Compute the gradient with respect to \code{par}?
#' @param phi Optional warm-start nuisance parameters.
#' @param nonnegative Force nonnegative measurement-model slope where applicable.
#' @return A list with \code{objective}, \code{covariance} (the commute-time
#'   \code{E}), \code{phi}, and (if \code{gradient}) \code{gradient}.
#' @keywords internal
#' @export
terradish_directed_algorithm <- function(gen, g, data, S, par, nu = NULL,
                                         gradient = TRUE, phi = NULL,
                                         nonnegative = TRUE)
{
  n <- gen$n
  focal <- data$demes
  nf <- length(focal)
  ed <- gen$edges; ea <- ed[, 1]; eb <- ed[, 2]
  rate <- gen$rates(par)

  G <- Matrix::sparseMatrix(i = ea, j = eb, x = rate, dims = c(n, n))
  Matrix::diag(G) <- 0
  G <- G - Matrix::Diagonal(n, x = Matrix::rowSums(G))

  Jc <- diag(nf) - matrix(1 / nf, nf, nf)
  H <- matrix(0, n, nf)
  hcache <- vector("list", nf)
  for (fj in seq_len(nf)) {
    j <- focal[fj]
    idx <- seq_len(n)[-j]
    Q <- G[idx, idx, drop = FALSE]
    hred <- tryCatch(as.numeric(Matrix::solve(Q, rep(-1, n - 1L))),
                     error = function(e) NULL)
    if (is.null(hred))
      stop("Directed generator is numerically singular at these parameters ",
           "(very strong asymmetry). Reduce the directional effect or bounds.",
           call. = FALSE)
    hfull <- numeric(n); hfull[idx] <- hred
    H[, fj] <- hfull
    hcache[[fj]] <- list(j = j, idx = idx, Q = Q, hfull = hfull)
  }
  Hf <- H[focal, , drop = FALSE]
  R <- Hf + t(Hf)
  E <- as.matrix(-0.5 * Jc %*% R %*% Jc)

  if (is.null(g))                                  # forward only: return E
    return(list(covariance = E, hcache = hcache, Hf = Hf))

  sub <- radish_subproblem(g, E, S, nu = nu, phi = phi, nonnegative = nonnegative,
                           control = NewtonRaphsonControl(verbose = FALSE,
                                                          ftol = 1e-10, ctol = 1e-10))
  out <- list(objective = sub$loglikelihood, covariance = E, phi = sub$phi,
              boundary = sub$boundary)
  if (!gradient) return(out)

  dL_dE <- as.matrix(sub$gradient)
  dL_dR <- -0.5 * (Jc %*% dL_dE %*% Jc)
  dL_dH <- dL_dR + t(dL_dR)                       # R = Hf + Hf^T
  dL_drate <- numeric(length(ea))
  for (fj in seq_len(nf)) {
    hc <- hcache[[fj]]; j <- hc$j; idx <- hc$idx
    b <- numeric(n); b[focal] <- dL_dH[, fj]; b[j] <- 0
    ared <- as.numeric(Matrix::solve(Matrix::t(hc$Q), b[idx]))   # transpose solve
    afull <- numeric(n); afull[idx] <- ared
    keep <- ea != j
    dL_drate[keep] <- dL_drate[keep] +
      afull[ea[keep]] * (hc$hfull[ea[keep]] - hc$hfull[eb[keep]])
  }
  # chain to (theta, gamma)
  grad_theta <- as.numeric(crossprod(gen$sab, dL_drate * rate))
  grad_gamma <- as.numeric(crossprod(gen$D,   dL_drate * rate))
  out$gradient <- c(grad_theta, grad_gamma)
  names(out$gradient) <- c(gen$theta_names, gen$gamma_names)
  out
}

#' Directional (non-reversible) conductance surface
#'
#' Fits a conductance surface in which gene flow can be \strong{directional}:
#' the movement rate from cell \eqn{a} to neighbor \eqn{b} need not equal the
#' reverse rate.  Movement is modeled as a covariate-parameterized continuous-
#' time Markov generator
#' \deqn{\log G_{a\to b} = \tfrac12(\eta_a + \eta_b) + d_{ab}^\top \gamma,\quad \eta_k = x_k^\top\theta,}
#' where \eqn{\theta} are the usual (symmetric) conductance effects and
#' \eqn{\gamma} are directional effects on antisymmetric edge covariates
#' \eqn{d_{ab} = -d_{ba}} (e.g. elevation drop, flow, wind; see
#' \code{\link{edge_gradient}}).  The genetic response is linked through the
#' symmetric \strong{commute time} of the chain, so any \code{terradish}
#' measurement model applies unchanged.  With \eqn{\gamma = 0} the generator is
#' reversible and the model reduces to isolation by resistance.
#'
#' @param formula Genetic response matrix \code{~} symmetric conductance
#'   covariates (as in \code{\link{terradish}}).
#' @param data A \code{terradish_graph} from \code{\link{conductance_surface}}.
#' @param directional Directional edge covariates from \code{\link{edge_gradient}}
#'   (or a list with \code{edges} and \code{d}).
#' @param measurement_model A \code{terradish} measurement model
#'   (e.g. \code{\link{generalized_wishart}}, \code{\link{leastsquares}},
#'   \code{\link{mlpe}}).
#' @param nu Effective Wishart degrees of freedom, if required by the model.
#' @param gamma_bound Symmetric bound on the directional coefficients during
#'   optimization (guards the ill-conditioned strong-asymmetry regime).
#'   Default \code{2}.
#' @param theta_bound Symmetric bound on the conductance coefficients. Default
#'   \code{5}.
#' @param nonnegative Force nonnegative measurement-model slope where applicable.
#' @param control Passed to the inner nuisance-parameter fits.
#'
#' @details
#' \strong{Scope.} This is the tractable directional model (directed commute
#' time): a single-lineage hitting-time formulation, reducing to resistance
#' distance when reversible.  It is not the full structured-coalescent model
#' (intractable on fine rasters).  Whether directional effects are identifiable
#' depends on the data: asymmetry leaves a signature in the symmetric commute
#' time on bounded landscapes, but it can be weak; report \eqn{\gamma} with its
#' uncertainty and check sensitivity.
#'
#' @return An object of class \code{"terradish_directed"} with the symmetric
#'   conductance estimates \code{theta}, the directional estimates \code{gamma}
#'   (with standard errors), the profiled nuisance parameters \code{phi}, the
#'   maximized log-likelihood, and the fit details.
#'
#' @seealso \code{\link{edge_gradient}}, \code{\link{terradish}},
#'   \code{\link{conductance_surface}}
#' @export
terradish_directed <- function(formula, data, directional,
                               measurement_model = generalized_wishart,
                               nu = NULL, gamma_bound = 2, theta_bound = 5,
                               nonnegative = TRUE,
                               control = NewtonRaphsonControl(verbose = FALSE))
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(inherits(measurement_model, c("terradish_measurement_model",
                                          "radish_measurement_model")))

  tm <- terms(formula)
  response <- attr(tm, "response")
  if (!response) stop("'formula' must have the genetic response matrix on the LHS")
  S <- as.matrix(eval(attr(tm, "variables")[[response + 1L]], parent.frame()))
  rhs <- if (length(attr(tm, "term.labels"))) reformulate(attr(tm, "term.labels")) else formula(~1)

  gen <- .directed_generator(rhs, data, directional)
  p <- gen$p; q <- gen$q
  par0 <- rep(0, p + q)
  lower <- c(rep(-theta_bound, p), rep(-gamma_bound, q))
  upper <- c(rep(theta_bound, p),  rep(gamma_bound, q))

  phi_state <- new.env(parent = emptyenv()); phi_state$value <- NULL
  fn <- function(par) {
    r <- tryCatch(terradish_directed_algorithm(gen, measurement_model, data, S, par,
                                               nu = nu, gradient = FALSE,
                                               phi = phi_state$value,
                                               nonnegative = nonnegative),
                  error = function(e) NULL)
    if (is.null(r)) return(1e12)
    phi_state$value <- r$phi
    r$objective
  }
  gr <- function(par) {
    r <- tryCatch(terradish_directed_algorithm(gen, measurement_model, data, S, par,
                                               nu = nu, gradient = TRUE,
                                               phi = phi_state$value,
                                               nonnegative = nonnegative),
                  error = function(e) NULL)
    if (is.null(r)) return(rep(0, p + q))
    phi_state$value <- r$phi
    r$gradient
  }

  opt <- optim(par0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = 300))

  final <- terradish_directed_algorithm(gen, measurement_model, data, S, opt$par,
                                        nu = nu, gradient = TRUE,
                                        nonnegative = nonnegative)
  # asymptotic vcov via numerical Hessian of the profiled objective (low-dim)
  # Hessian of the negative log-likelihood via the Jacobian of the analytic
  # gradient (cheap + accurate: reuses the transpose-solve adjoint rather than
  # finite-differencing the objective).
  pnames <- c(gen$theta_names, gen$gamma_names)
  H <- if (!requireNamespace("numDeriv", quietly = TRUE)) NULL else
    tryCatch({
      Hm <- numDeriv::jacobian(gr, opt$par)
      (Hm + t(Hm)) / 2
    }, error = function(e) NULL)
  vcov <- if (is.null(H)) matrix(NA_real_, p + q, p + q) else .safe_invert(H)
  dimnames(vcov) <- list(pnames, pnames)
  se <- sqrt(pmax(diag(vcov), 0)); names(se) <- pnames

  theta <- opt$par[seq_len(p)]; names(theta) <- gen$theta_names
  gamma <- opt$par[p + seq_len(q)]; names(gamma) <- gen$gamma_names
  loglik <- -final$objective
  n_phi <- length(final$phi)
  df <- p + q + n_phi
  eta <- as.numeric(assemble_model_matrix(rhs, data$x) %*% theta)  # symmetric log-conductance

  out <- list(call = match.call(), formula = formula,
              theta = theta, gamma = gamma, se = se, vcov = vcov,
              phi = final$phi, loglik = loglik, df = df,
              aic = -2 * loglik + 2 * df,
              covariance = final$covariance, fitted = final$covariance,
              response = S, convergence = opt$convergence,
              measurement_model = .terradish_measurement_model_name(measurement_model),
              nu = nu, npar = c(theta = p, gamma = q),
              logconductance = eta,
              dim = c(vertices = gen$n, focal = length(data$demes)))
  class(out) <- c("terradish_directed", "terradish")
  out
}

#' @export
print.terradish_directed <- function(x, ...)
{
  cat("Directional (non-reversible) conductance surface\n")
  cat(sprintf("  measurement model: %s%s\n", x$measurement_model,
              if (is.null(x$nu)) "" else sprintf("  (nu = %g)", x$nu)))
  cat(sprintf("  loglik = %.3f\n", x$loglik))
  tab <- cbind(Estimate = c(x$theta, x$gamma),
               `Std. Error` = x$se)
  cat("\nSymmetric conductance (theta) and directional (gamma) effects:\n")
  print(round(tab, 4))
  cat("\n(gamma = 0 => reversible / isolation-by-resistance)\n")
  invisible(x)
}

#' @export
coef.terradish_directed <- function(object, ...) c(object$theta, object$gamma)

#' @export
logLik.terradish_directed <- function(object, ...) {
  val <- object$loglik; attr(val, "df") <- object$df; class(val) <- "logLik"; val
}

#' @export
AIC.terradish_directed <- function(object, ..., k = 2) {
  if (identical(k, 2)) return(object$aic)
  -2 * object$loglik + k * object$df
}

#' @export
vcov.terradish_directed <- function(object, ...) object$vcov

#' @export
confint.terradish_directed <- function(object, parm, level = 0.95, ...) {
  est <- c(object$theta, object$gamma); se <- object$se
  z <- qnorm(1 - (1 - level) / 2)
  ci <- cbind(est - z * se, est + z * se)
  colnames(ci) <- paste0(round(100 * c((1 - level) / 2, 1 - (1 - level) / 2), 1), "%")
  rownames(ci) <- names(est)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
summary.terradish_directed <- function(object, conf.level = 0.95, ...) {
  est <- c(object$theta, object$gamma); se <- object$se
  z <- est / se
  zt <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
              `Pr(>|z|)` = pmin(2 * (1 - pnorm(abs(z))), 1))
  out <- list(ztable = zt, loglik = object$loglik, df = object$df, aic = object$aic,
              phi = object$phi, npar = object$npar,
              measurement_model = object$measurement_model, nu = object$nu,
              dim = object$dim, call = object$call)
  class(out) <- "summary.terradish_directed"
  out
}

#' @export
print.summary.terradish_directed <- function(x, digits = max(3L, getOption("digits") - 3L),
                                             signif.stars = getOption("show.signif.stars"), ...) {
  cat("Directional (non-reversible) conductance surface\n")
  cat("Call:   ", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat(sprintf("Measurement model: %s%s\n", x$measurement_model,
              if (is.null(x$nu)) "" else sprintf(" (nu = %g)", x$nu)))
  cat(sprintf("Loglik: %.3f  (df = %d)   AIC: %.2f\n\n", x$loglik, x$df, x$aic))
  cat("Coefficients (theta = symmetric conductance; gamma = directional):\n")
  printCoefmat(x$ztable, digits = digits, signif.stars = signif.stars, na.print = "NA")
  cat("\n(gamma = 0 => reversible / isolation-by-resistance)\n")
  invisible(x)
}

#' @export
plot.terradish_directed <- function(x, data, type = c("conductance", "logconductance"), ...) {
  type <- match.arg(type)
  if (missing(data) || is.null(data$stack))
    stop("Pass the `terradish_graph` (with a stored raster stack) as `data` to plot.", call. = FALSE)
  vals <- if (identical(type, "conductance")) exp(x$logconductance) else x$logconductance
  template <- data$stack[[1]]
  tv <- terra::values(template, dataframe = FALSE)[, 1]
  missing <- is.na(tv); tv[!missing] <- vals
  r <- terra::setValues(template, tv); names(r) <- paste0("symmetric_", type)
  terra::plot(r, main = "Symmetric conductance (reversible part)", ...)
  invisible(r)
}
