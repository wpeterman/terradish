# =====================================================================
# DRAGON: Directional Resistance and Asymmetric Gene-flow On Networks.
#
# Pure-R reference engine (Phase 1) for estimating ASYMMETRIC gene flow as a
# function of directional spatial covariates, fit to a deme-level genetic
# covariance by the structured coalescent (FRAME's Strobeck forward map) with a
# covariate-parameterized directed generator and an analytic adjoint gradient.
#
# This is the coalescent-correct successor to terradish_directed(): the directed
# generator is the same idea, but the model distance is the structured-coalescent
# expected pairwise coalescence time T (solved from Strobeck's equation), not the
# symmetric commute time. See DRAGON_DESIGN.md.
#
# Model:  M_{a->b} = exp(s_ab' theta + d_ab' gamma)            (directed migration)
#         L = D - M,  D = diag(rowSums M)
#         gamma_i (coalescence rate):
#            "uniform" : exp(c0)
#            "drift"   : exp(c0 + kappa * z_i)                 (Tier-1 drift surface)
#            "coupled" : exp(c0) * pi_i^(-alpha)               (FRAME stationary tie)
#         Strobeck:  L T + T L^T + diag(gamma_i T_ii) = J
#         Gower:     Sigma_c = -1/2 tau (C T_obs C') + nugget I
#         Wishart:   l = 1/2 p [ log|Sigma_c| + tr(Sigma_c^{-1} Shat_c) ]
#
# Performance note: this reference engine assembles and factorizes the
# d^2 x d^2 sparse Strobeck operator (and its transpose for the adjoint) per
# evaluation. Correct and adequate to ~50 demes; the O(d^4) Lyapunov backend
# (Phase 2, C++/RcppEigen) is the planned speed path.
# =====================================================================

# ---- orthonormal contrast basis: (o-1) x o, rows span 1^perp, C C' = I -------
.dragon_contrast <- function(o)
{
  H <- diag(o) - 1 / o
  e <- eigen(H, symmetric = TRUE)
  ord <- order(e$values, decreasing = TRUE)[seq_len(o - 1L)]
  t(e$vectors[, ord, drop = FALSE])
}

# ---- directed edge list (both directions) from undirected edge_pairs ---------
.dragon_directed_edges <- function(edge_pairs)
{
  ep <- as.matrix(edge_pairs)
  storage.mode(ep) <- "integer"
  list(a = c(ep[, 1], ep[, 2]),
       b = c(ep[, 2], ep[, 1]))
}

# ---- directed migration matrix M and Laplacian L from edge designs -----------
# Sdesign, Ddesign: (n_directed_edge x p) and (x q) covariate matrices; either may
# have zero columns. par packs (theta[p], gamma[q]) for the migration linear pred.
.dragon_ML <- function(theta, gamma, model)
{
  d <- model$d
  lin <- numeric(length(model$ea))
  if (model$q > 0) lin <- lin + as.numeric(model$Ddesign %*% gamma)
  if (model$p > 0) lin <- lin + as.numeric(model$Sdesign %*% theta)
  w <- exp(lin)
  M <- Matrix::sparseMatrix(i = model$ea, j = model$eb, x = w,
                            dims = c(d, d))
  L <- Matrix::Diagonal(x = Matrix::rowSums(M)) - M
  list(M = M, L = L, w = w)
}

# ---- Strobeck operator A (column-major vec of L T + T L^T + diag inject) ------
.dragon_strobeck_op <- function(L, gamma)
{
  d <- nrow(L)
  Ld <- as(L, "CsparseMatrix")
  Id <- as(Matrix::Diagonal(d), "CsparseMatrix")
  Aop <- Matrix::kronecker(Id, Ld) + Matrix::kronecker(Ld, Id)  # (I (x) L)+(L (x) I)
  pidx <- (seq_len(d) - 1L) * d + seq_len(d)                    # vec index of (i,i)
  inj <- numeric(d * d); inj[pidx] <- gamma
  Aop + Matrix::Diagonal(x = inj)
}

# ---- stationary distribution of the generator -L: pi = B^{-1} 1, B = L'+11' ---
.dragon_stationary <- function(L)
{
  d <- nrow(L)
  B <- t(as.matrix(L)) + matrix(1, d, d)
  as.numeric(solve(B, rep(1, d)))
}

# ---- unpack the parameter vector by coalescence mode -------------------------
.dragon_unpack <- function(par, model)
{
  p <- model$p; q <- model$q; mode <- model$mode
  theta <- if (p > 0) par[seq_len(p)] else numeric(0)
  gamma <- par[p + seq_len(q)]
  c0    <- par[p + q + 1L]
  off   <- p + q + 1L
  if (mode == "uniform") {
    extra <- NA_real_; ltau <- par[off + 1L]; lnug <- par[off + 2L]
  } else {
    extra <- par[off + 1L]; ltau <- par[off + 2L]; lnug <- par[off + 3L]
  }
  list(theta = theta, gamma = gamma, c0 = c0, extra = extra,
       ltau = ltau, lnug = lnug)
}

# ---- forward map: returns objective (neg loglik) + cache for the gradient -----
.dragon_forward <- function(par, model)
{
  d <- model$d; pr <- .dragon_unpack(par, model)
  tau <- exp(pr$ltau); nug <- exp(pr$lnug)
  ml <- .dragon_ML(pr$theta, pr$gamma, model)
  pivec <- NULL
  gam <- switch(model$mode,
    uniform = rep(exp(pr$c0), d),
    drift   = exp(pr$c0 + pr$extra * model$z),
    coupled = {
      pivec <- .dragon_stationary(ml$L)
      if (any(pivec <= 0) || any(!is.finite(pivec)))
        return(list(objective = 1e12, ok = FALSE))
      exp(pr$c0) * pivec^(-pr$extra)
    })
  Aop <- .dragon_strobeck_op(ml$L, gam)
  Tvec <- tryCatch(as.numeric(Matrix::solve(Aop, rep(1, d * d))),
                   error = function(e) NULL)
  if (is.null(Tvec)) return(list(objective = 1e12, ok = FALSE))
  Tm <- matrix(Tvec, d, d); Tm <- 0.5 * (Tm + t(Tm))
  Tobs <- Tm[model$obs, model$obs, drop = FALSE]
  Co <- model$Co
  Sig <- -0.5 * tau * (Co %*% Tobs %*% t(Co)) + nug * diag(length(model$obs) - 1L)
  ld <- determinant(Sig, logarithm = TRUE)
  if (ld$sign <= 0) return(list(objective = 1e12, ok = FALSE))
  Sinv <- solve(Sig)
  nll <- 0.5 * model$nu * (as.numeric(ld$modulus) + sum(Sinv * model$Shat_c))
  list(objective = nll, ok = TRUE, Aop = Aop, L = ml$L, w = ml$w, Tm = Tm,
       Tobs = Tobs, Sig = Sig, Sinv = Sinv, gam = gam, pivec = pivec,
       tau = tau, nug = nug, pr = pr)
}

# ---- analytic adjoint gradient ------------------------------------------------
.dragon_grad <- function(par, model)
{
  f <- .dragon_forward(par, model)
  if (!isTRUE(f$ok)) return(list(objective = f$objective,
                                 gradient = rep(0, length(par))))
  d <- model$d; Co <- model$Co; pr <- f$pr
  G <- f$Sinv - f$Sinv %*% model$Shat_c %*% f$Sinv
  Tbar <- matrix(0, d, d)
  Tbar[model$obs, model$obs] <- -0.25 * model$nu * f$tau * (t(Co) %*% G %*% Co)
  Lamvec <- as.numeric(Matrix::solve(Matrix::t(f$Aop), as.numeric(Tbar)))
  Lam <- matrix(Lamvec, d, d)
  Phi <- (Lam + t(Lam)) %*% f$Tm
  gbar <- -diag(Lam) * diag(f$Tm)                 # dl/dgamma_i
  Lbar <- -Phi
  dc0 <- sum(gbar * f$gam)
  dextra <- NULL
  if (model$mode == "drift") {
    dextra <- sum(gbar * f$gam * model$z)
  } else if (model$mode == "coupled") {
    pivec <- f$pivec
    dextra <- sum(gbar * f$gam * (-log(pivec)))
    pbar <- gbar * f$gam * (-pr$extra) / pivec
    q <- solve(as.matrix(f$L) + matrix(1, d, d), pbar)   # (L + 11') q = pbar
    Lbar <- Lbar - outer(pivec, q)
  }
  ge <- f$w * (Lbar[cbind(model$ea, model$ea)] - Lbar[cbind(model$ea, model$eb)])
  dgamma <- if (model$q > 0) as.numeric(t(model$Ddesign) %*% ge) else numeric(0)
  dtheta <- if (model$p > 0) as.numeric(t(model$Sdesign) %*% ge) else numeric(0)
  dtau <- -0.25 * model$nu * sum(G * (Co %*% f$Tobs %*% t(Co)))
  dnug <- 0.5 * model$nu * sum(diag(G))
  dltau <- dtau * f$tau; dlnug <- dnug * f$nug
  grad <- c(dtheta, dgamma, dc0,
            if (model$mode != "uniform") dextra else NULL,
            dltau, dlnug)
  list(objective = f$objective, gradient = grad)
}

# ---- multi-start L-BFGS-B fit -------------------------------------------------
.dragon_fit_engine <- function(model, n_start = 8, alpha_bounds = c(0, 1),
                               control = list(maxit = 400), seed = 1L)
{
  np <- model$p + model$q + 3L + (model$mode != "uniform")
  lower <- rep(-Inf, np); upper <- rep(Inf, np)
  if (model$mode == "coupled") {                 # alpha sits at p+q+2
    ai <- model$p + model$q + 2L
    lower[ai] <- alpha_bounds[1]; upper[ai] <- alpha_bounds[2]
  }
  set.seed(seed)
  best <- NULL
  for (s in seq_len(n_start)) {
    x0 <- stats::rnorm(np, 0, 0.6)
    if (model$mode == "coupled") x0[model$p + model$q + 2L] <-
      stats::runif(1, alpha_bounds[1], alpha_bounds[2])
    fit <- tryCatch(stats::optim(x0,
        fn = function(p) .dragon_grad(p, model)$objective,
        gr = function(p) .dragon_grad(p, model)$gradient,
        method = "L-BFGS-B", lower = lower, upper = upper, control = control),
      error = function(e) NULL)
    if (!is.null(fit) && fit$value < 1e11 &&
        (is.null(best) || fit$value < best$value)) best <- fit
  }
  best
}

# =====================================================================
# Internal model builder + the user-facing dragon()
# =====================================================================

# Build a low-level DRAGON model object from explicit pieces (testable without
# terra). edge_pairs: undirected integer matrix; elev: node directional potential;
# Sx: optional node symmetric covariate matrix; z: optional node drift covariate;
# Shat: node-level genetic covariance (d x d); obs: 1-based observed-deme indices.
.dragon_model <- function(edge_pairs, elev, Shat, nu, obs = NULL,
                          z = NULL, Sx = NULL, dedge = NULL, mode = "uniform")
{
  d <- nrow(Shat)
  de <- .dragon_directed_edges(edge_pairs)
  ea <- de$a; eb <- de$b
  # directional design: the gradient of a scalar node potential (curl-free,
  # reversible) and/or supplied antisymmetric edge covariates `dedge` (which can
  # carry a circulation/curl component, making the generator non-reversible).
  Dd <- if (!is.null(elev)) matrix(elev[ea] - elev[eb], ncol = 1L) else
        matrix(0, length(ea), 0L)
  if (!is.null(dedge))
    Dd <- cbind(Dd, as.matrix(dedge))
  Sd <- if (is.null(Sx)) matrix(0, length(ea), 0L) else
        0.5 * (as.matrix(Sx)[ea, , drop = FALSE] + as.matrix(Sx)[eb, , drop = FALSE])
  if (is.null(obs)) obs <- seq_len(d)
  Co <- .dragon_contrast(length(obs))
  Shat_c <- Co %*% Shat[obs, obs, drop = FALSE] %*% t(Co)
  list(d = d, ea = ea, eb = eb, Ddesign = Dd, Sdesign = Sd,
       p = ncol(Sd), q = ncol(Dd), z = if (is.null(z)) rep(0, d) else as.numeric(z),
       obs = obs, Co = Co, Shat_c = Shat_c, nu = nu, mode = mode)
}

#' Estimate asymmetric gene flow as a function of directional covariates (DRAGON)
#'
#' Fits the covariate-parameterized directed structured-coalescent model: a
#' directed migration generator whose rates are a log-linear function of
#' directional spatial covariates (elevation drop, downstream flow, prevailing
#' wind), with the genetic data modeled through the structured-coalescent expected
#' pairwise coalescence time and a Wishart likelihood. This is the
#' coalescent-correct, directional successor to \code{\link{terradish_directed}}.
#'
#' @param directional Optional numeric node-level potential (one value per deme,
#'   ordered as the graph's active demes) whose gradient across edges drives
#'   directional gene flow, for example elevation, stream order, or a wind
#'   potential. The per-edge directional covariate is \eqn{d_{ab} = elev_a - elev_b}.
#'   A scalar potential is curl-free, so the generator it induces is reversible and
#'   its stationary distribution is collinear with the potential; supply
#'   \code{circulation} for a non-reversible (rotational) component. At least one of
#'   \code{directional} or \code{circulation} must be given.
#' @param data A \code{terradish_graph} from \code{\link{conductance_surface}}.
#' @param S Observed genetic covariance matrix among demes (e.g. from
#'   \code{\link{cov_from_biallelic}}), with rows/columns ordered as the graph's
#'   active demes.
#' @param nu Effective Wishart degrees of freedom (the marker count).
#' @param coalescence Coalescence-rate model: \code{"uniform"} (default),
#'   \code{"drift"} (a covariate effective-size surface; supply \code{drift}), or
#'   \code{"coupled"} (FRAME stationary tie \code{gamma_i = c * pi_i^(-alpha)}).
#' @param drift Optional node-level covariate (length = number of demes) for the
#'   \code{"drift"} coalescence model.
#' @param circulation Optional antisymmetric edge covariate(s) carrying a
#'   rotational/curl component of directed gene flow that no scalar potential can
#'   express (so the directed generator becomes non-reversible and its stationary
#'   distribution need not be collinear with \code{directional}). A numeric vector
#'   with one entry per undirected edge in \code{data$edge_pairs}, or a matrix with
#'   one such column per circulation covariate; the value for edge
#'   \code{(a, b) = edge_pairs[k, ]} applies to the directed edge \eqn{a \to b} and
#'   its negative to \eqn{b \to a}.
#' @param obs Optional 1-based indices of observed demes (defaults to all).
#' @param n_start Number of random restarts for the optimizer.
#' @return An object of class \code{"dragon"} with the fitted coefficients, the
#'   maximized log-likelihood, and the model object. Compare nested models (e.g.
#'   \code{coalescence = "uniform"} vs \code{"drift"}) by \code{\link{AIC}} or a
#'   likelihood-ratio test, and check direction-vs-drift identifiability with
#'   \code{\link{dragon_collinearity}}.
#' @seealso \code{\link{terradish_directed}}, \code{\link{dragon_collinearity}}
#' @export
dragon <- function(directional = NULL, data, S, nu,
                   coalescence = c("uniform", "drift", "coupled"),
                   drift = NULL, circulation = NULL, obs = NULL, n_start = 8L)
{
  coalescence <- match.arg(coalescence)
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  ep <- data$edge_pairs
  if (is.null(ep)) ep <- t(data$adj) + 1L
  d <- nrow(as.matrix(S))
  # directional potential at active demes (a node-level covariate whose gradient
  # across edges drives directional gene flow, e.g. elevation or stream order).
  elev <- NULL
  if (!is.null(directional)) {
    elev <- as.numeric(directional)
    if (length(elev) != d)
      stop("`directional` must give one value per deme (", d, ").", call. = FALSE)
  }
  # circulation: antisymmetric edge covariate(s) over undirected edges -> directed
  # design c(forward, -reverse). Carries the non-gradient (curl) part of flow.
  dedge <- NULL; ncirc <- 0L
  if (!is.null(circulation)) {
    cm <- as.matrix(circulation)
    if (nrow(cm) != nrow(as.matrix(ep)))
      stop("`circulation` must have one row per undirected edge in `data$edge_pairs` (",
           nrow(as.matrix(ep)), ").", call. = FALSE)
    dedge <- rbind(cm, -cm)
    ncirc <- ncol(cm)
  }
  if (is.null(elev) && is.null(dedge))
    stop("supply `directional` (node potential) and/or `circulation` (edge covariate).",
         call. = FALSE)
  z <- if (!is.null(drift)) as.numeric(drift) else NULL
  model <- .dragon_model(ep, elev, as.matrix(S), nu = nu, obs = obs,
                         z = z, dedge = dedge, mode = coalescence)
  fit <- .dragon_fit_engine(model, n_start = n_start)
  if (is.null(fit)) stop("DRAGON optimization failed to converge.", call. = FALSE)
  pr <- .dragon_unpack(fit$par, model)
  gnm <- character(0)
  if (!is.null(elev)) gnm <- "gamma_dir"
  if (ncirc > 0L)
    gnm <- c(gnm, if (ncirc == 1L) "gamma_circ" else paste0("gamma_circ", seq_len(ncirc)))
  nm <- c(gnm,
          "c0",
          if (coalescence == "drift") "kappa" else if (coalescence == "coupled") "alpha",
          "log_tau", "log_nugget")
  out <- list(coefficients = stats::setNames(fit$par, nm),
              loglik = -fit$value, npar = length(fit$par),
              coalescence = coalescence, model = model, optim = fit)
  class(out) <- "dragon"
  out
}

#' @export
print.dragon <- function(x, ...)
{
  cat("DRAGON structured-coalescent fit (coalescence =", x$coalescence, ")\n")
  cat("  log-likelihood:", format(x$loglik, digits = 6),
      "   AIC:", format(2 * x$npar - 2 * x$loglik, digits = 6), "\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

#' @export
coef.dragon <- function(object, ...) object$coefficients
#' @export
logLik.dragon <- function(object, ...) {
  ll <- object$loglik; attr(ll, "df") <- object$npar; class(ll) <- "logLik"; ll
}
#' @export
summary.dragon <- function(object, ...) { print(object); invisible(object) }

#' Direction-versus-drift collinearity diagnostic for a DRAGON analysis
#'
#' Warns when a directional covariate and a drift (coalescence-rate) covariate are
#' too aligned to separately identify directional gene flow from a density /
#' effective-size gradient on a given dataset. Reports three signals: the covariate
#' correlation (\code{r_cov}), the disagreement between the drift-blind and
#' drift-joint directional estimates (\code{gap}, the practical confounding signal),
#' and a verdict. See \code{DRAGON_DESIGN.md} Section 7.3.
#'
#' @inheritParams dragon
#' @param drift Node-level drift covariate to test against the directional covariate.
#' @return A list with \code{r_cov}, \code{gdir_blind}, \code{gdir_joint},
#'   \code{gap}, and a \code{verdict} string.
#' @export
dragon_collinearity <- function(directional, data, S, nu, drift,
                                gap_warn = 0.15)
{
  elev <- as.numeric(directional); z <- as.numeric(drift)
  r_cov <- stats::cor(elev, z)
  fb <- dragon(directional, data, S, nu, coalescence = "uniform")
  fj <- dragon(directional, data, S, nu, coalescence = "drift", drift = drift)
  gb <- unname(fb$coefficients["gamma_dir"]); gj <- unname(fj$coefficients["gamma_dir"])
  gap <- abs(gb - gj) / max(abs(gb), 0.1)
  verdict <- if (abs(r_cov) < 0.7) {
    "OK: low covariate collinearity; fit jointly."
  } else if (gap > gap_warn) {
    paste("CONFOUNDED: covariates collinear and blind vs joint gamma_dir disagree.",
          "Report both as bounds, flag collinearity, seek independent N_e.")
  } else {
    "OK: covariates collinear but the coalescent separates them; joint fit trustworthy."
  }
  list(r_cov = r_cov, gdir_blind = gb, gdir_joint = gj, gap = gap, verdict = verdict)
}
