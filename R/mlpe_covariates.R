#' Pairwise endpoint-difference covariates
#'
#' Constructs pairwise environmental-mismatch covariates from site-level
#' values.  The result can be used alongside isolation-by-resistance (IBR)
#' terms to model isolation by environment (IBE) within
#' \code{\link{mlpe_covariates}}.
#'
#' @param x Site-level covariates.  Supported inputs:
#'   \itemize{
#'     \item A numeric vector (one covariate per site).
#'     \item A numeric matrix or data frame (sites in rows, covariates in
#'       columns).
#'     \item A \code{terra::SpatRaster}: values are extracted at \code{coords}.
#'   }
#' @param coords Required when \code{x} is a raster.  Focal-point coordinates
#'   in the same projection as \code{x}; accepts the same inputs as
#'   \code{\link{conductance_surface}}.
#' @param transform How to collapse site-level values into a pairwise
#'   environmental-mismatch measure.
#'   \describe{
#'     \item{\code{"absdiff"}}{Absolute difference between endpoints:
#'       \eqn{|x_i - x_j|}.  One output column per covariate.  Symmetric;
#'       suitable when you expect environmental dissimilarity (regardless of
#'       direction) to drive differentiation.}
#'     \item{\code{"sqdiff"}}{Squared difference:
#'       \eqn{(x_i - x_j)^2}.  One output column per covariate.  Like
#'       \code{"absdiff"} but penalizes large differences more heavily.}
#'     \item{\code{"euclidean"}}{Euclidean distance across all covariates:
#'       \eqn{\sqrt{\sum_k (x_{ik} - x_{jk})^2}}.  Collapses multivariate
#'       data to a single dissimilarity measure.}
#'     \item{\code{"manhattan"}}{Manhattan (city-block) distance:
#'       \eqn{\sum_k |x_{ik} - x_{jk}|}.  Collapses to a single measure;
#'       less sensitive to outliers than Euclidean.}
#'   }
#' @param scale Logical.  If \code{TRUE}, standardize (z-score) site-level
#'   covariates before computing the pairwise transform.  Recommended when
#'   covariates are on different scales.
#'
#' @details
#' The function produces a \eqn{n(n-1)/2 \times p} matrix of pairwise
#' endpoint-difference covariates, where \eqn{n} is the number of focal sites
#' and \eqn{p} is the number of output columns (one per input covariate for
#' \code{"absdiff"} and \code{"sqdiff"}; one for \code{"euclidean"} and
#' \code{"manhattan"}).  Rows correspond to the lower triangle of the
#' pairwise distance matrix, scanned column by column.
#'
#' The returned object has class \code{"terradish_pairwise_covariates"}.  This
#' class carries the original site-level covariate matrix and the chosen
#' transform as attributes, so that \code{\link{terradish_cv}} can
#' automatically rebuild the correct pairwise covariates for each train/test
#' split.
#'
#' If \code{x} already has class \code{"terradish_pairwise_covariates"}, it is
#' returned unchanged (pass-through).
#'
#' @return A numeric matrix of class \code{"terradish_pairwise_covariates"}
#'   with one row per unordered site pair and named columns reflecting the
#'   chosen transform (e.g. \code{"absdiff_altitude"},
#'   \code{"euclidean"}).
#'
#' @seealso \code{\link{mlpe_covariates}}, \code{\link{terradish_cv}},
#'   \code{\link{terradish}}
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' z_abs <- pairwise_endpoint_covariates(melip.altitude, melip.coords,
#'                                       transform = "absdiff", scale = TRUE)
#' head(z_abs)
#'
#' site_env <- data.frame(
#'   altitude = scale(as.numeric(terra::extract(melip.altitude, melip.coords)[, 2])),
#'   xcoord = scale(terra::geom(melip.coords)[, "x"])
#' )
#' z_env <- pairwise_endpoint_covariates(site_env, transform = "euclidean")
#' head(z_env)
#'
#' @export
pairwise_endpoint_covariates <- function(x,
                                         coords = NULL,
                                         transform = c("absdiff", "sqdiff",
                                                       "euclidean", "manhattan"),
                                         scale = FALSE)
{
  if (inherits(x, c("terradish_pairwise_covariates",
                    "radish_pairwise_covariates")))
    return(x)

  transform <- match.arg(transform)
  site_covariates <- .pairwise_site_covariates(x, coords = coords, scale = scale)
  .make_pairwise_endpoint_covariates(site_covariates, transform = transform)
}

#' MLPE measurement model with pairwise endpoint-difference covariates
#'
#' Creates an MLPE measurement model whose mean structure includes both an IBR
#' term (resistance distance) and one or more fixed pairwise environmental
#' covariates (IBE terms).
#'
#' @param x Pairwise endpoint-difference covariates (the output of
#'   \code{\link{pairwise_endpoint_covariates}}), \emph{or} site-level
#'   environmental data in any format accepted by that helper (numeric vector,
#'   matrix, data frame, or \code{terra::SpatRaster}).
#' @param coords Optional focal-point coordinates, required when \code{x} is
#'   a raster.  Accepts the same formats as \code{\link{conductance_surface}}.
#' @param transform Pairwise transform applied when \code{x} is not already a
#'   \code{"terradish_pairwise_covariates"} object.  Passed through to
#'   \code{\link{pairwise_endpoint_covariates}}.
#' @param scale Logical.  Standardize site-level covariates before the
#'   pairwise transform?  Passed to \code{\link{pairwise_endpoint_covariates}}.
#' @param rho_start Starting value for the MLPE correlation parameter \eqn{\rho}
#'   on its natural scale, which must lie in \code{(0, 0.5)}.  The default
#'   \code{0.2} corresponds to a moderate shared-site correlation.  Adjust if
#'   the optimizer fails to improve from the starting value.
#'
#' @details
#' The fitted mean structure is:
#'
#' \deqn{S_{ij} = \alpha + \beta R_{ij} + Z_{ij}^\top \gamma + e_{ij}}
#'
#' where \eqn{R_{ij}} is the resistance distance (IBR) and \eqn{Z_{ij}} is the
#' row of pairwise endpoint-difference covariates for the pair \eqn{(i,j)}
#' (IBE).  The residuals follow the MLPE correlation structure of Clarke et al.
#' (2002).
#'
#' The nuisance parameters (intercept \eqn{\alpha}, IBR slope \eqn{\beta},
#' IBE coefficients \eqn{\gamma}, log-precision \code{tau}, and logit MLPE
#' correlation \code{rho}) are all profiled out during optimization; only the
#' conductance parameters in the formula are optimized by the outer loop.
#'
#' When \code{x} is site-level data or the output of
#' \code{pairwise_endpoint_covariates()}, the pairwise covariate matrix is
#' stored as an attribute of the returned function.  This allows
#' \code{\link{terradish_cv}} to rebuild the correct pairwise covariates for
#' each train/test split automatically.
#'
#' @return A function of class \code{"terradish_measurement_model"} suitable
#'   for use as the \code{measurement_model} argument of \code{\link{terradish}},
#'   \code{\link{terradish_grid}}, and \code{\link{terradish_cv}}.
#'
#' @seealso \code{\link{mlpe}}, \code{\link{pairwise_endpoint_covariates}},
#'   \code{\link{terradish_cv}}
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' covariates <- c(terra::scale(melip.altitude), terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' g_joint <- mlpe_covariates(melip.altitude, melip.coords,
#'                            transform = "absdiff", scale = TRUE)
#' fit_joint <- terradish(melip.Fst ~ altitude + forestcover, surface,
#'                     terradish::loglinear_conductance, g_joint)
#' summary(fit_joint)
#'
#' @export
mlpe_covariates <- function(x,
                            coords = NULL,
                            transform = c("absdiff", "sqdiff",
                                          "euclidean", "manhattan"),
                            scale = FALSE,
                            rho_start = 0.2)
{
  transform <- match.arg(transform)
  if (!is.numeric(rho_start) || length(rho_start) != 1L ||
      is.na(rho_start) || rho_start <= 0 || rho_start >= 0.5)
    stop("`rho_start` must be a single number in (0, 0.5)")

  pairwise_covariates <- if (inherits(x, c("terradish_pairwise_covariates",
                                           "radish_pairwise_covariates")))
    x
  else
    pairwise_endpoint_covariates(x, coords = coords, transform = transform,
                                 scale = scale)

  Z <- .as_pairwise_covariate_matrix(pairwise_covariates)

  g <- function(E, S, phi, nu = NULL, gradient = TRUE, hessian = TRUE,
                partial = TRUE, nonnegative = TRUE, validate = FALSE)
  {
    symm <- function(X) (X + t(X))/2

    if (!(is.matrix(E) && is.matrix(S) && all(dim(E) == dim(S))))
      stop("invalid inputs")

    ones <- matrix(1, nrow(E), 1)
    Ed   <- diag(E)
    R    <- Ed %*% t(ones) + ones %*% t(Ed) - 2 * symm(E)
    Rl   <- matrix(R[lower.tri(R)], ncol = 1)
    Sl   <- matrix(S[lower.tri(S)], ncol = 1)

    if (nrow(Z) != length(Sl))
      stop("pairwise endpoint covariates do not match the number of focal-point pairs")

    unos <- matrix(1, length(Sl), 1)
    X    <- cbind(unos, Rl, Z)

    coef_names <- c("alpha", "beta", colnames(Z))
    phi_names  <- c(coef_names, "tau", "rho")

    if (missing(phi))
    {
      coef_start <- .mlpe_covariate_start(X, Sl, nonnegative = nonnegative)
      resid <- Sl - X %*% matrix(coef_start, ncol = 1)
      sigma2 <- max(mean(resid^2), .Machine$double.eps)
      phi <- c(coef_start, "tau" = log(1 / sigma2),
               "rho" = qlogis(2 * rho_start))

      lower <- c(-Inf, if (nonnegative) 0 else -Inf,
                 rep(-Inf, ncol(Z)), -Inf, -Inf)
      upper <- rep(Inf, length(phi))
      names(phi) <- phi_names
      return(list(phi = phi, lower = lower, upper = upper))
    }
    else if (!(is.numeric(phi) && length(phi) == length(phi_names)))
      stop("invalid inputs")

    names(phi) <- phi_names

    coef_vec <- matrix(phi[coef_names], ncol = 1)
    beta  <- coef_vec[2, 1]
    tau   <- exp(phi["tau"])
    rho   <- 0.5 * plogis(phi["rho"])

    Ind <- which(lower.tri(R), arr.ind = TRUE)
    U   <- sparseMatrix(i = rep(seq_len(length(Sl)), 2), j = c(Ind), x = 1)

    eigUtU <- .get_mlpe_eigen(nrow(E))
    D      <- eigUtU$values
    P      <- eigUtU$vectors
    Dr     <- D/(1 - 2 * rho) + 1/rho

    SigmaInv <- function(x)
    {
      Ax <- 1/(1 - 2 * rho) * x
      x  <- t(U) %*% Ax
      x  <- t(P) %*% x
      x  <- x / Dr
      x  <- P %*% x
      x  <- Ax - 1/(1 - 2 * rho) * U %*% x
      as.matrix(x)
    }

    SigmaLogDet <- sum(log(Dr)) + length(D) * log(rho) +
      length(Sl) * log(1 - 2 * rho)

    e      <- Sl - X %*% coef_vec
    Si_e   <- SigmaInv(e)
    loglik <- -0.5 * tau * t(e) %*% Si_e + 0.5 * nrow(e) * log(tau) -
      0.5 * SigmaLogDet

    fitted <- matrix(0, nrow(S), ncol(S))
    fitted[lower.tri(fitted)] <- as.vector(X %*% coef_vec)
    fitted <- fitted + t(fitted)
    diag(fitted) <- coef_vec[1, 1]

    if (gradient || hessian || partial)
    {
      p        <- ncol(X)
      dPhi     <- matrix(0, length(phi), 1,
                         dimnames = list(phi_names, NULL))
      ddPhi    <- matrix(0, length(phi), length(phi),
                         dimnames = list(phi_names, phi_names))
      ddEdPhi  <- matrix(0, length(Rl), length(phi),
                         dimnames = list(NULL, phi_names))
      ddPhidS  <- matrix(0, length(phi), length(Sl),
                         dimnames = list(phi_names, NULL))

      drho_Si_e  <- t(U) %*% Si_e
      drho_Si_e  <- as.matrix(2 * Si_e - U %*% drho_Si_e)
      drho_trans <- rho * (1 - 2 * rho)

      Si_X <- lapply(seq_len(p), function(j)
        SigmaInv(matrix(X[, j], ncol = 1)))
      names(Si_X) <- coef_names

      dPhi[coef_names, 1] <- tau * crossprod(X, Si_e)[, 1]
      dPhi["tau", 1]      <- -0.5 * tau * t(e) %*% Si_e + 0.5 * length(e)
      dPhi["rho", 1]      <-
        (-0.5 * tau * t(Si_e) %*% drho_Si_e -
         0.5 * sum((2 * D/(1 - 2 * rho)^2 - 1/rho^2)/Dr) -
         0.5 * length(D)/rho + length(Sl)/(1 - 2 * rho)) * drho_trans

      if (hessian || partial)
      {
        Si_drho_Si_e <- SigmaInv(drho_Si_e)

        for (j in seq_len(p))
        {
          xj <- matrix(X[, j], ncol = 1)
          for (k in j:p)
          {
            val <- -tau * crossprod(xj, Si_X[[k]])[1, 1]
            ddPhi[coef_names[j], coef_names[k]] <- val
            ddPhi[coef_names[k], coef_names[j]] <- val
          }
          ddPhi[coef_names[j], "tau"] <- tau * crossprod(xj, Si_e)[1, 1]
          ddPhi["tau", coef_names[j]] <- ddPhi[coef_names[j], "tau"]
          ddPhi[coef_names[j], "rho"] <-
            tau * crossprod(Si_X[[j]], drho_Si_e)[1, 1] * drho_trans
          ddPhi["rho", coef_names[j]] <- ddPhi[coef_names[j], "rho"]
        }

        ddPhi["tau", "tau"] <- -0.5 * tau * t(e) %*% Si_e
        ddPhi["tau", "rho"] <- -0.5 * tau * t(Si_e) %*% drho_Si_e * drho_trans
        ddPhi["rho", "tau"] <- ddPhi["tau", "rho"]
        ddPhi["rho", "rho"] <-
          (-tau * t(drho_Si_e) %*% Si_drho_Si_e -
           0.5 * sum((8 * D/(1 - 2 * rho)^3 + 2/rho^3)/Dr) +
           0.5 * sum((2 * D/(1 - 2 * rho)^2 - 1/rho^2)^2/Dr^2) +
           0.5 * length(D)/rho^2 + 2 * length(Sl)/(1 - 2 * rho)^2) *
          drho_trans^2 + dPhi["rho", 1] * (1 - 4 * rho)

        if (partial)
        {
          dR <- matrix(0, nrow(R), ncol(R))
          dR[lower.tri(dR)] <- 2 * beta * tau * as.vector(Si_e)
          dR <- symm(dR)
          dE <- diag(nrow(R)) * (dR %*% ones %*% t(ones)) - dR

          ddEdPhi[, "alpha"] <- -2 * beta * tau * Si_X[["alpha"]][, 1]
          ddEdPhi[, "beta"]  <- 2 * tau * (Si_e[, 1] - beta * Si_X[["beta"]][, 1])
          if (ncol(Z) > 0)
          {
            for (nm in colnames(Z))
              ddEdPhi[, nm] <- -2 * beta * tau * Si_X[[nm]][, 1]
          }
          ddEdPhi[, "tau"] <- 2 * beta * tau * Si_e[, 1]
          ddEdPhi[, "rho"] <- 2 * beta * tau * Si_drho_Si_e[, 1] * drho_trans
          ddEdPhi <- apply(ddEdPhi, 2, function(x) {
            Xmat <- matrix(0, nrow(E), ncol(E))
            Xmat[lower.tri(Xmat)] <- x
            Xmat <- symm(Xmat)
            diag(nrow(E)) * (Xmat %*% ones %*% t(ones)) - Xmat
          })

          for (j in seq_len(p))
            ddPhidS[coef_names[j], ] <- tau * Si_X[[j]][, 1]
          ddPhidS["tau", ] <- -tau * Si_e[, 1]
          ddPhidS["rho", ] <- -tau * Si_drho_Si_e[, 1] * drho_trans

          jacobian_E <- function(dE)
          {
            ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
            ddEdE[lower.tri(ddEdE)] <-
              SigmaInv(matrix(ddEdE[lower.tri(ddEdE)], ncol = 1))[, 1]
            ddEdE[upper.tri(ddEdE)] <- 0
            ddEdE <- ddEdE + t(ddEdE)
            ddEdE <- -beta^2 * tau * ddEdE
            ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
            -ddEdE
          }

          jacobian_S <- function(dE)
          {
            ddEdE <- diag(dE) %*% t(ones) + ones %*% t(diag(dE)) - 2 * symm(dE)
            ddEdE[lower.tri(ddEdE)] <-
              SigmaInv(matrix(ddEdE[lower.tri(ddEdE)], ncol = 1))[, 1]
            ddEdE[upper.tri(ddEdE)] <- 0
            ddEdE <- ddEdE + t(ddEdE)
            ddEdE <- -beta * tau * ddEdE
            ddEdE <- diag(nrow(dE)) * (ddEdE %*% ones %*% t(ones)) - ddEdE
            diag(ddEdE) <- 0
            -ddEdE
          }
        }
      }
    }

    list(objective  = -c(loglik),
         fitted     = fitted,
         boundary   = nonnegative && beta == 0,
         gradient   = if (!gradient) NULL else -dPhi,
         hessian    = if (!hessian)  NULL else -ddPhi,
         gradient_E = if (!partial)  NULL else -dE,
         partial_E  = if (!partial)  NULL else -ddEdPhi,
         partial_S  = if (!partial)  NULL else -ddPhidS,
         jacobian_E = if (!partial)  NULL else jacobian_E,
         jacobian_S = if (!partial)  NULL else jacobian_S)
  }

  attr(g, "pairwise_covariates") <- pairwise_covariates
  attr(g, "subsetter") <- function(index)
    mlpe_covariates(.subset_pairwise_endpoint_covariates(pairwise_covariates, index),
                    rho_start = rho_start)
  class(g) <- unique(c("terradish_measurement_model",
                       "radish_measurement_model",
                       class(g)))
  g
}

.pairwise_site_covariates <- function(x, coords = NULL, scale = FALSE)
{
  if (inherits(x, "PackedSpatRaster"))
    x <- unwrap(x)

  site_covariates <- if (inherits(x, "SpatRaster"))
  {
    if (is.null(coords))
      stop("`coords` must be supplied when `x` is a raster")
    pts <- .coords_matrix(coords, x)
    vals <- extract(x, pts)
    if ("ID" %in% colnames(vals))
      vals <- vals[, colnames(vals) != "ID", drop = FALSE]
    as.matrix(vals)
  }
  else if (is.vector(x) && is.atomic(x))
  {
    matrix(as.numeric(x), ncol = 1)
  }
  else if (is.data.frame(x) || is.matrix(x))
  {
    as.matrix(x)
  }
  else
  {
    stop("`x` must be a numeric vector, matrix, data.frame, or SpatRaster")
  }

  if (!is.numeric(site_covariates))
    stop("site-level covariates must be numeric")
  if (anyNA(site_covariates))
    stop("missing values are not supported in endpoint covariates")

  if (is.null(colnames(site_covariates)))
    colnames(site_covariates) <- paste0("var", seq_len(ncol(site_covariates)))

  if (isTRUE(scale))
  site_covariates <- scale(site_covariates)

  as.matrix(site_covariates)
}

.make_pairwise_endpoint_covariates <- function(site_covariates, transform)
{
  site_covariates <- as.matrix(site_covariates)
  n <- nrow(site_covariates)
  if (n < 2)
    stop("need at least two focal points to construct pairwise covariates")

  pairwise <- switch(
    transform,
    absdiff = {
      out <- matrix(NA_real_, n * (n - 1)/2, ncol(site_covariates))
      for (j in seq_len(ncol(site_covariates)))
      {
        diff_j <- outer(site_covariates[, j], site_covariates[, j], "-")
        out[, j] <- abs(diff_j[lower.tri(diff_j)])
      }
      colnames(out) <- paste0("absdiff_", colnames(site_covariates))
      out
    },
    sqdiff = {
      out <- matrix(NA_real_, n * (n - 1)/2, ncol(site_covariates))
      for (j in seq_len(ncol(site_covariates)))
      {
        diff_j <- outer(site_covariates[, j], site_covariates[, j], "-")
        out[, j] <- diff_j[lower.tri(diff_j)]^2
      }
      colnames(out) <- paste0("sqdiff_", colnames(site_covariates))
      out
    },
    euclidean = {
      out <- matrix(as.vector(dist(site_covariates, method = "euclidean")),
                    ncol = 1)
      colnames(out) <- "euclidean"
      out
    },
    manhattan = {
      out <- matrix(as.vector(dist(site_covariates, method = "manhattan")),
                    ncol = 1)
      colnames(out) <- "manhattan"
      out
    }
  )

  structure(pairwise,
            site_covariates = site_covariates,
            transform = transform,
            class = unique(c("terradish_pairwise_covariates",
                             "radish_pairwise_covariates",
                             class(pairwise))))
}

.subset_pairwise_endpoint_covariates <- function(x, index)
{
  if (!inherits(x, c("terradish_pairwise_covariates",
                     "radish_pairwise_covariates")))
    stop("`x` must inherit from 'terradish_pairwise_covariates'")

  site_covariates <- attr(x, "site_covariates")
  transform <- attr(x, "transform")
  .make_pairwise_endpoint_covariates(site_covariates[index, , drop = FALSE],
                                     transform = transform)
}

.as_pairwise_covariate_matrix <- function(x)
{
  Z <- as.matrix(x)
  if (!is.numeric(Z))
    stop("pairwise endpoint covariates must be numeric")
  if (anyNA(Z))
    stop("missing values are not supported in pairwise endpoint covariates")
  if (is.null(colnames(Z)))
    colnames(Z) <- paste0("cov", seq_len(ncol(Z)))
  Z
}

.mlpe_covariate_start <- function(X, y, nonnegative = TRUE)
{
  coef_start <- tryCatch(as.numeric(qr.solve(X, y)),
                         error = function(e) rep(0, ncol(X)))

  if (nonnegative && coef_start[2] < 0)
  {
    X0 <- X[, -2, drop = FALSE]
    coef0 <- tryCatch(as.numeric(qr.solve(X0, y)),
                      error = function(e) rep(0, ncol(X0)))
    coef_start <- c(coef0[1], 0, coef0[-1])
  }

  setNames(coef_start, colnames(X))
}
