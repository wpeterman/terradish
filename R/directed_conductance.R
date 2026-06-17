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
# Engine note: the default Matrix backend is correct and adequate for moderate
# graphs. The optional Eigen::SparseLU C++ backend caches one factorization per
# focal absorber and reuses it for the forward and transpose adjoint solves.
# =====================================================================

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
#' @param solver Linear solver backend. \code{"matrix"} uses the reference
#'   \pkg{Matrix} implementation; \code{"sparse_lu_cpp"} uses the Eigen
#'   SparseLU C++ backend and reuses each absorber factorization for the
#'   forward and transpose adjoint solves.
#' @return A list with \code{objective}, \code{covariance} (the commute-time
#'   \code{E}), \code{phi}, and (if \code{gradient}) \code{gradient}.
#' @keywords internal
#' @export
terradish_directed_algorithm <- function(gen, g, data, S, par, nu = NULL,
                                         gradient = TRUE, phi = NULL,
                                         nonnegative = TRUE,
                                         solver = c("matrix", "sparse_lu_cpp"))
{
  solver <- match.arg(solver)
  n <- gen$n
  focal <- data$demes
  nf <- length(focal)
  ed <- gen$edges; ea <- ed[, 1]; eb <- ed[, 2]
  rate <- gen$rates(par)

  Jc <- diag(nf) - matrix(1 / nf, nf, nf)

  if (identical(solver, "sparse_lu_cpp")) {
    fw <- tryCatch(directed_sparse_lu_forward_cpp(ed, rate, as.integer(focal), n),
                   error = function(e) e)
    if (inherits(fw, "error"))
      stop("Directed SparseLU backend failed at these parameters ",
           "(very strong asymmetry can make the generator singular): ",
           conditionMessage(fw), call. = FALSE)
    Hf <- fw$Hf
    hcache <- fw$hcache
  } else {
    G <- Matrix::sparseMatrix(i = ea, j = eb, x = rate, dims = c(n, n))
    Matrix::diag(G) <- 0
    G <- G - Matrix::Diagonal(n, x = Matrix::rowSums(G))

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
  }
  R <- Hf + t(Hf)
  E <- as.matrix(-0.5 * Jc %*% R %*% Jc)

  if (is.null(g))                                  # forward only: return E
    return(list(covariance = E, hcache = hcache, Hf = Hf, solver = solver))

  sub <- radish_subproblem(g, E, S, nu = nu, phi = phi, nonnegative = nonnegative,
                           control = NewtonRaphsonControl(verbose = FALSE,
                                                          ftol = 1e-10, ctol = 1e-10))
  out <- list(objective = sub$loglikelihood, covariance = E, phi = sub$phi,
              boundary = sub$boundary)
  if (!gradient) return(out)

  dL_dE <- as.matrix(sub$gradient)
  dL_dR <- -0.5 * (Jc %*% dL_dE %*% Jc)
  dL_dH <- dL_dR + t(dL_dR)                       # R = Hf + Hf^T
  if (identical(solver, "sparse_lu_cpp")) {
    dL_drate <- tryCatch(
      directed_sparse_lu_adjoint_cpp(hcache, ed, as.integer(focal), dL_dH, n),
      error = function(e) e)
    if (inherits(dL_drate, "error"))
      stop("Directed SparseLU adjoint failed: ", conditionMessage(dL_drate),
           call. = FALSE)
  } else {
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
#' @param solver Linear solver backend for the directed hitting-time systems.
#'   \code{"matrix"} uses the reference \pkg{Matrix} implementation.
#'   \code{"sparse_lu_cpp"} uses Eigen SparseLU through C++ and reuses each
#'   absorber factorization for the forward and transpose adjoint solves.
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
                               control = NewtonRaphsonControl(verbose = FALSE),
                               solver = c("matrix", "sparse_lu_cpp"))
{
  solver <- match.arg(solver)
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
                                               nonnegative = nonnegative,
                                               solver = solver),
                  error = function(e) NULL)
    if (is.null(r)) return(1e12)
    phi_state$value <- r$phi
    r$objective
  }
  gr <- function(par) {
    r <- tryCatch(terradish_directed_algorithm(gen, measurement_model, data, S, par,
                                               nu = nu, gradient = TRUE,
                                               phi = phi_state$value,
                                               nonnegative = nonnegative,
                                               solver = solver),
                  error = function(e) NULL)
    if (is.null(r)) return(rep(0, p + q))
    phi_state$value <- r$phi
    r$gradient
  }

  opt <- optim(par0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = 300))

  final <- terradish_directed_algorithm(gen, measurement_model, data, S, opt$par,
                                        nu = nu, gradient = TRUE,
                                        nonnegative = nonnegative,
                                        solver = solver)
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
  eta <- withCallingHandlers(
    as.numeric(assemble_model_matrix(rhs, data$x) %*% theta),  # symmetric log-conductance
    warning = function(w)
      if (grepl("Removed unused spatial covariates", conditionMessage(w)))
        invokeRestart("muffleWarning"))

  out <- list(call = match.call(), formula = formula,
              theta = theta, gamma = gamma, se = se, vcov = vcov,
              phi = final$phi, loglik = loglik, df = df,
              aic = -2 * loglik + 2 * df,
              covariance = final$covariance, fitted = final$covariance,
              response = S, convergence = opt$convergence,
              measurement_model = .terradish_measurement_model_name(measurement_model),
              nu = nu, npar = c(theta = p, gamma = q),
              solver = solver,
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
  if (!is.null(x$solver)) cat(sprintf("  solver: %s\n", x$solver))
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
              solver = object$solver,
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
  if (!is.null(x$solver)) cat(sprintf("Solver: %s\n", x$solver))
  cat(sprintf("Loglik: %.3f  (df = %d)   AIC: %.2f\n\n", x$loglik, x$df, x$aic))
  cat("Coefficients (theta = symmetric conductance; gamma = directional):\n")
  printCoefmat(x$ztable, digits = digits, signif.stars = signif.stars, na.print = "NA")
  cat("\n(gamma = 0 => reversible / isolation-by-resistance)\n")
  invisible(x)
}

#' Extract fitted directed edge rates
#'
#' Computes fitted forward and reverse edge rates from a
#' \code{\link{terradish_directed}} model. The returned table is useful for
#' diagnosing and plotting the directional component of the model: the
#' symmetric part describes where conductance is high or low, while
#' \code{log_rate_ratio} describes how strongly the fitted generator favors one
#' direction over the reverse direction along each graph edge.
#'
#' @param object A fitted \code{terradish_directed} object.
#' @param data The \code{terradish_graph} used to fit \code{object}.
#' @param directional Directional edge covariates from
#'   \code{\link{edge_gradient}}. This must be the same directional object used
#'   when fitting \code{object}.
#' @param level Confidence level used for the optional edge-level Wald summary
#'   of \code{log_rate_ratio}. The default is \code{0.95}.
#'
#' @return A data frame with one row per undirected graph edge. Columns include
#'   the original edge endpoints (\code{a}, \code{b}), endpoint coordinates,
#'   fitted rates in both directions, the geometric mean symmetric rate,
#'   arithmetic mean bidirectional rate, \code{log_rate_ratio =
#'   log(rate_ab / rate_ba)}, and arrow-ready coordinates from the lower-rate
#'   endpoint toward the higher-rate endpoint. If \code{object} contains a
#'   finite covariance matrix for \code{gamma}, standard errors, z statistics,
#'   p-values, and Wald significance indicators are also returned for
#'   \code{log_rate_ratio}.
#'
#' @seealso \code{\link{terradish_directed}}, \code{\link{edge_gradient}}
#' @export
directed_rates <- function(object, data, directional, level = 0.95)
{
  if (!inherits(object, "terradish_directed"))
    stop("`object` must be a fitted terradish_directed object.", call. = FALSE)
  if (!inherits(data, c("terradish_graph", "radish_graph")))
    stop("`data` must be the terradish_graph used to fit `object`.", call. = FALSE)
  if (missing(directional))
    stop("`directional` must be supplied; use the edge covariates passed to terradish_directed().",
         call. = FALSE)
  if (is.null(data$vertex_coordinates))
    stop("`data` must contain vertex coordinates to extract directed rates.",
         call. = FALSE)

  tm <- stats::terms(object$formula)
  labels <- attr(tm, "term.labels")
  rhs <- if (length(labels)) stats::reformulate(labels) else stats::as.formula("~1")
  gen <- .directed_generator(rhs, data, directional)
  q <- gen$q
  m <- nrow(gen$edges) / 2L
  if (m != floor(m))
    stop("Directed edge list must contain paired forward and reverse edges.",
         call. = FALSE)
  if (q != length(object$gamma))
    stop("`directional` does not match the fitted directional coefficients.",
         call. = FALSE)

  par <- c(object$theta, object$gamma)
  rates <- gen$rates(par)
  rate_ab <- rates[seq_len(m)]
  rate_ba <- rates[m + seq_len(m)]

  edge_pairs <- gen$edges[seq_len(m), , drop = FALSE]
  coords <- as.matrix(data$vertex_coordinates)
  xy_a <- coords[edge_pairs[, 1], , drop = FALSE]
  xy_b <- coords[edge_pairs[, 2], , drop = FALSE]
  colnames(xy_a) <- colnames(xy_b) <- c("x", "y")

  log_rate_ratio <- log(rate_ab) - log(rate_ba)
  forward_favored <- log_rate_ratio >= 0
  x <- ifelse(forward_favored, xy_a[, 1], xy_b[, 1])
  y <- ifelse(forward_favored, xy_a[, 2], xy_b[, 2])
  xend <- ifelse(forward_favored, xy_b[, 1], xy_a[, 1])
  yend <- ifelse(forward_favored, xy_b[, 2], xy_a[, 2])

  out <- data.frame(
    a = edge_pairs[, 1],
    b = edge_pairs[, 2],
    x_a = xy_a[, 1],
    y_a = xy_a[, 2],
    x_b = xy_b[, 1],
    y_b = xy_b[, 2],
    rate_ab = rate_ab,
    rate_ba = rate_ba,
    symmetric_rate = sqrt(rate_ab * rate_ba),
    average_rate = (rate_ab + rate_ba) / 2,
    log_rate_ratio = log_rate_ratio,
    abs_log_rate_ratio = abs(log_rate_ratio),
    favored_from = ifelse(forward_favored, edge_pairs[, 1], edge_pairs[, 2]),
    favored_to = ifelse(forward_favored, edge_pairs[, 2], edge_pairs[, 1]),
    x = x,
    y = y,
    xend = xend,
    yend = yend
  )

  gamma_start <- length(object$theta) + 1L
  gamma_idx <- gamma_start:(gamma_start + q - 1L)
  vc <- object$vcov
  if (!is.null(vc) && all(dim(vc) >= c(max(gamma_idx), max(gamma_idx))) &&
      all(is.finite(vc[gamma_idx, gamma_idx, drop = FALSE])))
  {
    D <- as.matrix(gen$D)
    grad <- D[seq_len(m), , drop = FALSE] - D[m + seq_len(m), , drop = FALSE]
    vc_g <- vc[gamma_idx, gamma_idx, drop = FALSE]
    se <- sqrt(pmax(rowSums((grad %*% vc_g) * grad), 0))
    z <- log_rate_ratio / se
    p <- pmin(2 * (1 - stats::pnorm(abs(z))), 1)
    crit <- stats::qnorm(1 - (1 - level) / 2)
    out$log_rate_ratio_se <- se
    out$z <- z
    out$p_value <- p
    out$significant <- is.finite(se) & se > 0 & abs(log_rate_ratio) > crit * se
  }
  else
  {
    out$log_rate_ratio_se <- NA_real_
    out$z <- NA_real_
    out$p_value <- NA_real_
    out$significant <- NA
  }

  out
}

#' @export
plot.terradish_directed <- function(x, data,
                                    type = c("conductance", "logconductance",
                                             "directional", "combined"),
                                    directional,
                                    min_abs_log_ratio = 0,
                                    significant_only = FALSE,
                                    level = 0.95,
                                    ...) {
  type <- match.arg(type)
  if (missing(data))
    stop("Pass the `terradish_graph` used for fitting as `data` to plot.",
         call. = FALSE)
  if (!inherits(data, c("terradish_graph", "radish_graph")))
    stop("`data` must be the terradish_graph used to fit `x`.", call. = FALSE)
  if (!identical(type, "directional") && is.null(data$stack))
    stop("Pass a `terradish_graph` with a stored raster stack for this plot type.",
         call. = FALSE)

  if (type %in% c("conductance", "logconductance"))
  {
    vals <- if (identical(type, "conductance")) exp(x$logconductance) else x$logconductance
    template <- data$stack[[1]]
    tv <- terra::values(template, dataframe = FALSE)[, 1]
    missing <- is.na(tv); tv[!missing] <- vals
    r <- terra::setValues(template, tv); names(r) <- paste0("symmetric_", type)
    terra::plot(r, main = "Symmetric conductance (reversible part)", ...)
    return(invisible(r))
  }

  if (missing(directional))
    stop('plot(type = "directional" or "combined") requires `directional`: ',
         'supply the edge covariates used for fitting.',
         call. = FALSE)

  edge_data <- directed_rates(x, data = data, directional = directional,
                              level = level)
  edge_data <- edge_data[is.finite(edge_data$abs_log_rate_ratio) &
                           edge_data$abs_log_rate_ratio >= min_abs_log_ratio, ,
                         drop = FALSE]
  if (isTRUE(significant_only))
    edge_data <- edge_data[edge_data$significant %in% TRUE, , drop = FALSE]

  if (nrow(edge_data) == 0)
    stop("No directed edges remain after filtering.", call. = FALSE)

  if (identical(type, "directional"))
    return(.plot_directed_edges(edge_data))

  .plot_directed_combined(x, data, edge_data)
}

.directed_symmetric_raster_data <- function(x, data)
{
  template <- data$stack[[1]]
  tv <- terra::values(template, dataframe = FALSE)[, 1]
  missing <- is.na(tv); tv[!missing] <- exp(x$logconductance)
  r <- terra::setValues(template, tv)
  out <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(out)[3] <- "conductance"
  out
}

.directed_focal_data <- function(data)
{
  coords <- as.matrix(data$vertex_coordinates)
  demes <- data$demes
  data.frame(x = coords[demes, 1], y = coords[demes, 2])
}

.plot_directed_edges <- function(edge_data)
{
  edge_data <- .shorten_directed_segments(edge_data)
  ggplot2::ggplot(edge_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_segment(
      ggplot2::aes(colour = abs_log_rate_ratio, linewidth = abs_log_rate_ratio),
      arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed"),
      lineend = "round"
    ) +
    ggplot2::coord_equal(expand = TRUE) +
    ggplot2::scale_colour_gradient(low = "grey82", high = "#9e0142",
                                   name = "|log rate ratio|") +
    ggplot2::scale_linewidth_continuous(range = c(0.15, 1.6),
                                        guide = "none") +
    ggplot2::scale_x_continuous(labels = .terradish_plot_number) +
    ggplot2::scale_y_continuous(labels = .terradish_plot_number) +
    ggplot2::labs(title = "Directional edge-rate bias", x = NULL, y = NULL) +
    .terradish_plot_theme()
}

.plot_directed_combined <- function(x, data, edge_data)
{
  background <- .directed_symmetric_raster_data(x, data)
  focal <- .directed_focal_data(data)
  edge_data <- .shorten_directed_segments(edge_data)

  ggplot2::ggplot() +
    ggplot2::geom_raster(data = background,
                         ggplot2::aes(x = x, y = y, fill = conductance)) +
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   colour = abs_log_rate_ratio,
                   linewidth = abs_log_rate_ratio),
      alpha = 0.8,
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed"),
      lineend = "round"
    ) +
    ggplot2::geom_point(data = focal, ggplot2::aes(x = x, y = y),
                        shape = 21, size = 1.8, stroke = 0.35,
                        fill = "white", colour = "black") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_fill_gradientn(colours = terrain.colors(100),
                                  name = "Conductance") +
    ggplot2::scale_colour_gradient(low = "grey20", high = "#9e0142",
                                   name = "|log rate ratio|") +
    ggplot2::scale_linewidth_continuous(range = c(0.12, 1.2),
                                        guide = "none") +
    ggplot2::scale_x_continuous(labels = .terradish_plot_number) +
    ggplot2::scale_y_continuous(labels = .terradish_plot_number) +
    ggplot2::labs(title = "Symmetric conductance and directional bias",
                  x = NULL, y = NULL) +
    .terradish_plot_theme()
}

.shorten_directed_segments <- function(edge_data, inset = 0.28)
{
  dx <- edge_data$xend - edge_data$x
  dy <- edge_data$yend - edge_data$y
  edge_data$x <- edge_data$x + inset * dx
  edge_data$y <- edge_data$y + inset * dy
  edge_data$xend <- edge_data$xend - inset * dx
  edge_data$yend <- edge_data$yend - inset * dy
  edge_data
}
