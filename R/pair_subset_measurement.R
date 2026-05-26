#' Pair-subset measurement models
#'
#' Build a measurement model that keeps all focal sites in the resistance
#' calculation but evaluates the likelihood using only selected pairs.
#'
#' @param measurement_model A regression-style measurement model. Currently
#'   \code{\link{leastsquares}} and \code{\link{mlpe}} are supported.
#' @param pairs A two-column matrix or data frame identifying the site pairs to
#'   use. Numeric pairs are interpreted as row/column positions in the response
#'   matrix. Character pairs are matched to \code{rownames(S)} at fit time.
#'
#' @details
#' This helper is useful when the sampling design or biological question calls
#' for retaining all sites in the landscape graph while fitting the measurement
#' model to a subset of pairwise genetic distances. The selected sites still
#' contribute right-hand sides to the Laplacian solve, but only the requested
#' pair rows contribute to the likelihood.
#'
#' The MLPE version preserves the pairwise dependence structure by constructing
#' the MLPE incidence matrix from the selected pairs only. Pair subsetting is
#' most natural for distance-regression measurement models; matrix-likelihood
#' models such as \code{\link{generalized_wishart}} and
#' \code{\link{wishart_covariance}} are intentionally not supported. Wishart
#' models evaluate a full distance or covariance matrix likelihood, so dropping
#' selected pair entries would change the model contract. For site-level
#' environmental effects in Wishart models, use
#' \code{\link{wishart_covariates}()} instead.
#'
#' @return A function of class \code{terradish_measurement_model}.
#'
#' @examples
#' data(melip)
#' selected_pairs <- rbind(
#'   c(1, 2),
#'   c(1, 3),
#'   c(2, 4),
#'   c(3, 5)
#' )
#' pair_mlpe <- pair_subset_measurement_model(mlpe, selected_pairs)
#' inherits(pair_mlpe, "terradish_measurement_model")
#'
#' @export
pair_subset_measurement_model <- function(measurement_model = mlpe, pairs)
{
  if (!inherits(measurement_model, c("terradish_measurement_model",
                                     "radish_measurement_model")))
    stop("`measurement_model` must be a terradish measurement model.", call. = FALSE)

  model_name <- .pair_subset_measurement_model_name(measurement_model)
  if (isTRUE(model_name %in% c("generalized_wishart", "wishart_covariance")))
    stop("`pair_subset_measurement_model()` cannot wrap `", model_name, "`. ",
         "Wishart models use a full-matrix likelihood, so selected pair rows ",
         "cannot be dropped without changing the model contract. Use ",
         "`wishart_covariates()` to add site-level covariance kernels for ",
         "`generalized_wishart` or `wishart_covariance` workflows.",
         call. = FALSE)

  if (!isTRUE(model_name %in% c("leastsquares", "mlpe")))
    stop("Pair-subset measurement models currently support only `leastsquares` and `mlpe`.",
         call. = FALSE)

  force(pairs)
  pair_cache <- new.env(parent = emptyenv())

  out <- function(E, S, phi, nu = NULL, gradient = TRUE, hessian = TRUE,
                  partial = TRUE, nonnegative = TRUE, validate = FALSE)
  {
    pair_info <- .pair_subset_info(pairs, E, S, pair_cache)

    if (identical(model_name, "leastsquares"))
      return(.leastsquares_pair_subset(E = E, S = S, pairs = pair_info$pairs,
                                       phi = if (missing(phi)) NULL else phi,
                                       nu = nu, gradient = gradient,
                                       hessian = hessian, partial = partial,
                                       nonnegative = nonnegative,
                                       validate = validate))

    .mlpe_pair_subset(E = E, S = S, pair_info = pair_info,
                      phi = if (missing(phi)) NULL else phi,
                      nu = nu, gradient = gradient,
                      hessian = hessian, partial = partial,
                      nonnegative = nonnegative,
                      validate = validate)
  }

  class(out) <- c("terradish_pair_subset_measurement_model",
                  "terradish_measurement_model",
                  "radish_measurement_model")
  attr(out, "base_model") <- model_name
  attr(out, "pairs") <- pairs
  out
}

.pair_subset_measurement_model_name <- function(model)
{
  known <- c("leastsquares", "mlpe", "generalized_wishart",
             "wishart_covariance")
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

.pair_subset_info <- function(pairs, E, S, cache)
{
  if (!(is.matrix(E) && is.matrix(S) && all(dim(E) == dim(S))))
    stop("`E` and `S` must be square matrices with matching dimensions.", call. = FALSE)

  n <- nrow(S)
  resolved <- .pair_subset_resolve_pairs(pairs, S)
  key <- paste(n, paste(c(t(resolved)), collapse = ","), sep = ":")
  cached <- get0(key, envir = cache, inherits = FALSE)
  if (!is.null(cached))
    return(cached)

  U <- sparseMatrix(i = rep(seq_len(nrow(resolved)), 2),
                    j = c(resolved[, 1], resolved[, 2]),
                    x = 1,
                    dims = c(nrow(resolved), n))
  eig <- eigen(as.matrix(t(U) %*% U))
  info <- list(pairs = resolved, U = U, eigen = eig, n_sites = n)
  assign(key, info, envir = cache)
  info
}

.pair_subset_resolve_pairs <- function(pairs, S)
{
  if (is.data.frame(pairs))
    pairs <- as.matrix(pairs)
  if (!is.matrix(pairs) || ncol(pairs) != 2L)
    stop("`pairs` must be a two-column matrix or data frame.", call. = FALSE)
  if (!nrow(pairs))
    stop("`pairs` must contain at least one pair.", call. = FALSE)

  n <- nrow(S)
  if (is.character(pairs))
  {
    site_names <- rownames(S)
    if (is.null(site_names))
      stop("Character `pairs` require row names on `S`.", call. = FALSE)
    pairs <- matrix(match(c(pairs), site_names), ncol = 2L)
    if (anyNA(pairs))
      stop("All character `pairs` must match row names on `S`.", call. = FALSE)
  }
  else
  {
    suppressWarnings(storage.mode(pairs) <- "integer")
  }

  if (anyNA(pairs) || any(pairs < 1L) || any(pairs > n))
    stop("`pairs` contain indices outside the response matrix.", call. = FALSE)
  if (any(pairs[, 1] == pairs[, 2]))
    stop("`pairs` cannot contain self-pairs.", call. = FALSE)

  # Store unordered pairs in lower-triangle order to match existing models.
  pairs <- cbind(pmax(pairs[, 1], pairs[, 2]),
                 pmin(pairs[, 1], pairs[, 2]))
  if (any(duplicated(paste(pairs[, 1], pairs[, 2], sep = ":"))))
    stop("`pairs` cannot contain duplicate unordered pairs.", call. = FALSE)

  pairs
}

.pair_subset_symm <- function(x)
{
  (x + t(x)) / 2
}

.pair_subset_distance <- function(E, pairs)
{
  E <- .pair_subset_symm(E)
  Ed <- diag(E)
  Ed[pairs[, 1]] + Ed[pairs[, 2]] - 2 * E[pairs]
}

.pair_subset_response <- function(S, pairs)
{
  S[pairs]
}

.pair_subset_backprop_distance <- function(n, pairs, q)
{
  q <- as.vector(q)
  if (length(q) != nrow(pairs))
    stop("Internal pair derivative dimension mismatch.", call. = FALSE)

  dR <- matrix(0, n, n)
  dR[pairs] <- 2 * q
  dR <- .pair_subset_symm(dR)
  ones <- matrix(1, n, 1)
  diag(n) * (dR %*% ones %*% t(ones)) - dR
}

.pair_subset_backprop_matrix <- function(n, pairs, q)
{
  q <- as.matrix(q)
  out <- matrix(0, n * n, ncol(q))
  for (k in seq_len(ncol(q)))
    out[, k] <- c(.pair_subset_backprop_distance(n, pairs, q[, k]))
  colnames(out) <- colnames(q)
  out
}

.pair_subset_response_matrix <- function(n, pairs, value)
{
  value <- as.vector(value)
  out <- matrix(0, n, n)
  out[pairs] <- value
  out[cbind(pairs[, 2], pairs[, 1])] <- value
  out
}

.pair_subset_response_partial <- function(n, pairs, value)
{
  value <- as.matrix(value)
  out <- matrix(0, nrow(value), n * n)
  idx1 <- pairs[, 1] + (pairs[, 2] - 1L) * n
  idx2 <- pairs[, 2] + (pairs[, 1] - 1L) * n
  out[, idx1] <- value
  out[, idx2] <- value
  out
}

.leastsquares_pair_subset <- function(E, S, pairs, phi = NULL, nu = NULL,
                                      gradient = TRUE, hessian = TRUE,
                                      partial = TRUE, nonnegative = TRUE,
                                      validate = FALSE)
{
  if (is.null(phi))
  {
    Rl <- .pair_subset_distance(E, pairs)
    Sl <- .pair_subset_response(S, pairs)
    fit <- gls(Sl ~ Rl, method = "ML")

    if (!nonnegative || coef(fit)[2] > 0)
    {
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = phi[2], "tau" = -2 * log(sigma(fit)))
    }
    else
    {
      fit <- gls(Sl ~ 1, method = "ML")
      phi <- coef(fit)
      names(phi) <- NULL
      phi <- c("alpha" = phi[1], "beta" = 0, "tau" = -2 * log(sigma(fit)))
    }

    return(list(phi = phi,
                lower = if (nonnegative) c(-Inf, 0, -Inf) else c(-Inf, -Inf, -Inf),
                upper = c(Inf, Inf, Inf)))
  }

  if (!(is.numeric(phi) && length(phi) == 3L))
    stop("invalid inputs", call. = FALSE)

  names(phi) <- c("alpha", "beta", "tau")
  n <- nrow(E)
  alpha <- phi["alpha"]
  beta <- phi["beta"]
  tau <- exp(phi["tau"])

  Rl <- .pair_subset_distance(E, pairs)
  Sl <- .pair_subset_response(S, pairs)
  unos <- matrix(1, length(Sl), 1)
  e <- Sl - alpha * unos - beta * Rl
  loglik <- -0.5 * tau * t(e) %*% e + 0.5 * nrow(e) * log(tau)

  Ed <- diag(.pair_subset_symm(E))
  Rfull <- Ed %*% t(matrix(1, n, 1)) + matrix(1, n, 1) %*% t(Ed) - 2 * .pair_subset_symm(E)
  fitted <- alpha + beta * Rfull

  if (gradient || hessian || partial)
  {
    dPhi <- matrix(0, length(phi), 1)
    ddPhi <- matrix(0, length(phi), length(phi))
    rownames(dPhi) <- colnames(ddPhi) <- rownames(ddPhi) <- names(phi)

    dPhi["alpha", ] <- t(unos) %*% e * tau
    dPhi["beta", ] <- t(Rl) %*% e * tau
    dPhi["tau", ] <- -0.5 * tau * t(e) %*% e + 0.5 * length(e)

    if (hessian || partial)
    {
      ddPhi["alpha", "alpha"] <- -t(unos) %*% unos * tau
      ddPhi["alpha", "beta"] <- -tau * t(unos) %*% Rl
      ddPhi["alpha", "tau"] <- tau * t(unos) %*% e
      ddPhi["beta", "beta"] <- -t(Rl) %*% Rl * tau
      ddPhi["beta", "tau"] <- tau * t(Rl) %*% e
      ddPhi["tau", "tau"] <- -0.5 * tau * t(e) %*% e
      ddPhi <- ddPhi + t(ddPhi)
      diag(ddPhi) <- diag(ddPhi) / 2

      if (partial)
      {
        q <- beta * tau * as.vector(e)
        dE <- .pair_subset_backprop_distance(n, pairs, q)
        dq_dphi <- cbind(
          "alpha" = -beta * tau * as.vector(unos),
          "beta" = tau * as.vector(e - beta * Rl),
          "tau" = beta * tau * as.vector(e)
        )
        ddEdPhi <- .pair_subset_backprop_matrix(n, pairs, dq_dphi)
        ddPhidS <- rbind("alpha" = as.vector(unos) * tau,
                         "beta" = as.vector(Rl) * tau,
                         "tau" = -tau * as.vector(e))
        ddPhidS <- .pair_subset_response_partial(n, pairs, ddPhidS)

        jacobian_E <- function(dE)
        {
          dR <- .pair_subset_distance(dE, pairs)
          - .pair_subset_backprop_distance(n, pairs, -beta^2 * tau * dR)
        }

        jacobian_S <- function(dE)
        {
          dR <- .pair_subset_distance(dE, pairs)
          .pair_subset_response_matrix(n, pairs, -beta * tau * dR)
        }
      }
    }
  }

  list(objective = -c(loglik),
       fitted = fitted,
       boundary = nonnegative && beta == 0,
       gradient = if (!gradient) NULL else -dPhi,
       hessian = if (!hessian) NULL else -ddPhi,
       gradient_E = if (!partial) NULL else -dE,
       partial_E = if (!partial) NULL else -ddEdPhi,
       partial_S = if (!partial) NULL else -ddPhidS,
       jacobian_E = if (!partial) NULL else jacobian_E,
       jacobian_S = if (!partial) NULL else jacobian_S,
       num_gradient = NULL,
       num_hessian = NULL,
       num_gradient_E = NULL,
       num_partial_E = NULL,
       num_partial_S = NULL,
       num_jacobian_E = NULL,
       num_jacobian_S = NULL)
}

.mlpe_pair_subset <- function(E, S, pair_info, phi = NULL, nu = NULL,
                              gradient = TRUE, hessian = TRUE,
                              partial = TRUE, nonnegative = TRUE,
                              validate = FALSE)
{
  pairs <- pair_info$pairs
  if (is.null(phi))
  {
    ls_start <- .leastsquares_pair_subset(E = E, S = S, pairs = pairs,
                                          nonnegative = nonnegative)$phi
    phi <- c(ls_start, "rho" = qlogis(0.2))
    return(list(phi = phi,
                lower = if (nonnegative) c(-Inf, 0, -Inf, -Inf)
                        else c(-Inf, -Inf, -Inf, -Inf),
                upper = c(Inf, Inf, Inf, Inf)))
  }

  if (!(is.numeric(phi) && length(phi) == 4L))
    stop("invalid inputs", call. = FALSE)

  names(phi) <- c("alpha", "beta", "tau", "rho")
  n <- nrow(E)
  alpha <- phi["alpha"]
  beta <- phi["beta"]
  tau <- exp(phi["tau"])
  rho <- 0.5 * plogis(phi["rho"])

  Rl <- .pair_subset_distance(E, pairs)
  Sl <- .pair_subset_response(S, pairs)
  U <- pair_info$U
  D <- pair_info$eigen$values
  P <- pair_info$eigen$vectors
  Dr <- D / (1 - 2 * rho) + 1 / rho

  SigmaInv <- function(x)
  {
    x <- as.matrix(x)
    Ax <- 1 / (1 - 2 * rho) * x
    y <- as.matrix(t(U) %*% Ax)
    y <- as.matrix(t(P) %*% y)
    y <- y / Dr
    y <- as.matrix(P %*% y)
    Ax - 1 / (1 - 2 * rho) * as.matrix(U %*% y)
  }

  SigmaLogDet <- sum(log(Dr)) + length(D) * log(rho) + length(Sl) * log(1 - 2 * rho)

  unos <- matrix(1, length(Sl), 1)
  e <- Sl - alpha * unos - beta * Rl
  Si_e <- SigmaInv(e)
  loglik <- -0.5 * tau * t(e) %*% Si_e + 0.5 * nrow(e) * log(tau) - 0.5 * SigmaLogDet

  Ed <- diag(.pair_subset_symm(E))
  Rfull <- Ed %*% t(matrix(1, n, 1)) + matrix(1, n, 1) %*% t(Ed) - 2 * .pair_subset_symm(E)
  fitted <- alpha + beta * Rfull

  if (gradient || hessian || partial)
  {
    dPhi <- matrix(0, length(phi), 1)
    ddPhi <- matrix(0, length(phi), length(phi))
    rownames(dPhi) <- colnames(ddPhi) <- rownames(ddPhi) <- names(phi)

    drho_Si_e <- t(U) %*% Si_e
    drho_Si_e <- as.matrix(2 * Si_e - U %*% drho_Si_e)
    drho_trans <- rho * (1 - 2 * rho)

    dPhi["alpha", ] <- t(unos) %*% Si_e * tau
    dPhi["beta", ] <- t(Rl) %*% Si_e * tau
    dPhi["tau", ] <- -0.5 * tau * t(e) %*% Si_e + 0.5 * length(e)
    dPhi["rho", ] <-
      (-0.5 * tau * t(Si_e) %*% drho_Si_e -
         0.5 * sum((2 * D / (1 - 2 * rho)^2 - 1 / rho^2) / Dr) -
         0.5 * length(D) / rho + length(Sl) / (1 - 2 * rho)) * drho_trans

    if (hessian || partial)
    {
      Si_unos <- SigmaInv(unos)
      Si_Rl <- SigmaInv(Rl)
      Si_drho_Si_e <- SigmaInv(drho_Si_e)

      ddPhi["alpha", "alpha"] <- -tau * t(unos) %*% Si_unos
      ddPhi["alpha", "beta"] <- -tau * t(unos) %*% Si_Rl
      ddPhi["alpha", "tau"] <- tau * t(unos) %*% Si_e
      ddPhi["alpha", "rho"] <- tau * t(Si_unos) %*% drho_Si_e * drho_trans
      ddPhi["beta", "beta"] <- -tau * t(Rl) %*% Si_Rl
      ddPhi["beta", "tau"] <- tau * t(Rl) %*% Si_e
      ddPhi["beta", "rho"] <- tau * t(Si_Rl) %*% drho_Si_e * drho_trans
      ddPhi["tau", "tau"] <- -0.5 * tau * t(e) %*% Si_e
      ddPhi["tau", "rho"] <- -0.5 * tau * t(Si_e) %*% drho_Si_e * drho_trans
      ddPhi["rho", "rho"] <-
        (-tau * t(drho_Si_e) %*% Si_drho_Si_e +
           -0.5 * sum((8 * D / (1 - 2 * rho)^3 + 2 / rho^3) / Dr) +
           0.5 * sum((2 * D / (1 - 2 * rho)^2 - 1 / rho^2)^2 / Dr^2) +
           0.5 * length(D) / rho^2 + 2 * length(Sl) / (1 - 2 * rho)^2) *
        drho_trans^2 + dPhi["rho", ] * (1 - 4 * rho)
      ddPhi <- ddPhi + t(ddPhi)
      diag(ddPhi) <- diag(ddPhi) / 2

      if (partial)
      {
        q <- beta * tau * as.vector(Si_e)
        dE <- .pair_subset_backprop_distance(n, pairs, q)
        dq_dphi <- cbind(
          "alpha" = -beta * tau * as.vector(Si_unos),
          "beta" = tau * as.vector(Si_e - beta * Si_Rl),
          "tau" = beta * tau * as.vector(Si_e),
          "rho" = beta * tau * as.vector(Si_drho_Si_e) * drho_trans
        )
        ddEdPhi <- .pair_subset_backprop_matrix(n, pairs, dq_dphi)
        ddPhidS <- rbind("alpha" = tau * as.vector(Si_unos),
                         "beta" = tau * as.vector(Si_Rl),
                         "tau" = -tau * as.vector(Si_e),
                         "rho" = -tau * as.vector(Si_drho_Si_e) * drho_trans)
        ddPhidS <- .pair_subset_response_partial(n, pairs, ddPhidS)

        jacobian_E <- function(dE)
        {
          dR <- .pair_subset_distance(dE, pairs)
          - .pair_subset_backprop_distance(n, pairs,
                                           -beta^2 * tau * as.vector(SigmaInv(dR)))
        }

        jacobian_S <- function(dE)
        {
          dR <- .pair_subset_distance(dE, pairs)
          .pair_subset_response_matrix(n, pairs,
                                       -beta * tau * as.vector(SigmaInv(dR)))
        }
      }
    }
  }

  list(objective = -c(loglik),
       fitted = fitted,
       boundary = nonnegative && beta == 0,
       gradient = if (!gradient) NULL else -dPhi,
       hessian = if (!hessian) NULL else -ddPhi,
       gradient_E = if (!partial) NULL else -dE,
       partial_E = if (!partial) NULL else -ddEdPhi,
       partial_S = if (!partial) NULL else -ddPhidS,
       jacobian_E = if (!partial) NULL else jacobian_E,
       jacobian_S = if (!partial) NULL else jacobian_S,
       num_gradient = NULL,
       num_hessian = NULL,
       num_gradient_E = NULL,
       num_partial_E = NULL,
       num_partial_S = NULL,
       num_jacobian_E = NULL,
       num_jacobian_S = NULL)
}
