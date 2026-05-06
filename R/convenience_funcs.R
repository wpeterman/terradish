#' Rescale values to the unit interval
#'
#' Rescales numeric values so the minimum non-missing value maps to 0 and
#' the maximum maps to 1.  \code{terra::SpatRaster} inputs are rescaled
#' layer-by-layer using only the finite (non-\code{NA}) values.
#'
#' @param x Numeric vector, matrix, \code{dist} object, or
#'   \code{terra::SpatRaster}.  \code{PackedSpatRaster} objects (from
#'   \code{terra::wrap()}) are unwrapped automatically.
#'
#' @details
#' The transformation applied to each element is
#' \eqn{(x - \min) / (\max - \min)}, computed over non-\code{NA} values.
#' \code{NA} values are propagated unchanged.
#'
#' If all non-missing values are identical (constant layer or vector), the
#' output is set to 0 for all non-missing entries rather than producing
#' \code{NaN} from a zero-range denominator.
#'
#' For rasters, each layer is rescaled independently using its own minimum and
#' maximum.
#'
#' @return An object of the same class as \code{x} with values in
#'   \code{[0, 1]} (non-\code{NA} entries only).
#'
#' @seealso \code{\link{scale_covariates}} for standardized (z-score or min-max)
#'   scaling with metadata retained for back-transformation.
#'
#' @examples
#' scale_to_0_1(c(2, 4, 6))
#' scale_to_0_1(c(5, 5, NA))  # constant non-NA values become 0
#'
#' library(terra)
#' r <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 3, 5, NA))
#' scale_to_0_1(r)
#'
#' @export
scale_to_0_1 <- function(x)
{
  if (inherits(x, "PackedSpatRaster"))
    x <- unwrap(x)

  if (inherits(x, "SpatRaster"))
  {
    vals <- values(x, dataframe = FALSE)
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0)
    {
      vals[!is.na(vals)] <- 0
    }
    else
    {
      vals <- (vals - rng[1]) / diff(rng)
    }
    values(x) <- vals
    return(x)
  }

  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0)
  {
    out <- x
    out[!is.na(out)] <- 0
    return(out)
  }

  (x - rng[1]) / diff(rng)
}

#' Scale raster covariates and retain original-scale metadata
#'
#' Transforms each layer of a \code{SpatRaster} and stores the original-scale
#' offset and divisor as a \code{"terradish_scale"} attribute.  When the scaled
#' raster is used to build a \code{\link{conductance_surface}}, marginal-effect
#' plots (\code{plot(fit, type = "marginal", data = surface)}) can
#' automatically back-transform x-axes to the original covariate units.
#'
#' @param covariates A \code{SpatRaster} or \code{PackedSpatRaster} whose
#'   layers should be scaled.
#' @param method Scaling method applied to each layer independently.
#'   \code{"zscore"} (default) subtracts the layer mean and divides by the
#'   standard deviation, giving zero mean and unit variance.
#'   \code{"minmax"} maps the minimum to 0 and the maximum to 1.
#' @param center Logical.  Only used when \code{method = "zscore"}.  If
#'   \code{TRUE} (default), subtract the layer mean before dividing by the
#'   standard deviation.
#' @param scale Logical.  Only used when \code{method = "zscore"}.  If
#'   \code{TRUE} (default), divide by the layer standard deviation.  If
#'   \code{FALSE} with \code{center = TRUE}, layers are centered but not
#'   divided.
#'
#' @details
#' Scaling is strongly recommended before fitting conductance models.  In the
#' loglinear model, conductance is \eqn{\exp(\theta \cdot x)}.  When \code{x}
#' has a large range (e.g. altitude in metres), a small change in \eqn{\theta}
#' produces enormous conductance differences, making the objective surface
#' poorly conditioned and the optimizer slow to converge.
#'
#' Missing values (\code{NA}) are excluded from the computation of the mean,
#' standard deviation, minimum, and range, and are left as \code{NA} in the
#' output.
#'
#' If a layer is constant (all values identical), the scale is set to 1 to
#' avoid division by zero; the layer is then returned unchanged (or centered
#' to 0 if \code{center = TRUE}).
#'
#' The \code{"terradish_scale"} attribute is a named list with one element per
#' layer; each element is a named numeric vector with entries \code{center} and
#' \code{scale} (the subtracted offset and the divisor, respectively).
#' \code{crop_to_focal_buffer} and \code{conductance_surface} preserve this
#' attribute.
#'
#' @return A \code{SpatRaster} with transformed values and a
#'   \code{"terradish_scale"} attribute.
#'
#' @seealso \code{\link{scale_to_0_1}}, \code{\link{conductance_surface}}
#'
#' @examples
#' library(terra)
#' r <- c(terra::rast(nrows = 2, ncols = 2, vals = c(100, 200, 300, NA)),
#'        terra::rast(nrows = 2, ncols = 2, vals = c(0.2, 0.5, 0.4, NA)))
#' names(r) <- c("altitude", "forestcover")
#'
#' scaled_r <- scale_covariates(r)
#' attr(scaled_r, "terradish_scale")  # per-layer center and scale
#'
#' # min-max scaling to [0, 1]
#' scaled_mm <- scale_covariates(r, method = "minmax")
#'
#' @export
scale_covariates <- function(covariates,
                             method = c("zscore", "minmax"),
                             center = TRUE,
                             scale = TRUE)
{
  covariates <- .as_spatraster(covariates)
  method <- match.arg(method)
  stats <- lapply(seq_len(nlyr(covariates)), function(i) {
    vals <- values(covariates[[i]], dataframe = FALSE)[, 1]
    vals <- vals[is.finite(vals)]

    if (method == "minmax")
    {
      center_i <- min(vals)
      scale_i <- diff(range(vals))
      if (!is.finite(scale_i) || scale_i == 0)
        scale_i <- 1
    }
    else
    {
      center_i <- if (isTRUE(center)) mean(vals) else 0
      scale_i <- if (isTRUE(scale)) sd(vals) else 1
      if (!is.finite(scale_i) || scale_i == 0)
        scale_i <- 1
    }

    c(center = center_i, scale = scale_i)
  })
  names(stats) <- names(covariates)

  vals <- values(covariates, dataframe = FALSE)
  for (nm in names(stats))
  {
    vals[, nm] <- (vals[, nm] - stats[[nm]][["center"]]) /
      stats[[nm]][["scale"]]
  }
  values(covariates) <- vals
  attr(covariates, "terradish_scale") <- stats
  covariates
}

#' Extract the strict lower triangle of a matrix
#'
#' Returns the entries strictly below the main diagonal of a square matrix,
#' scanned column-by-column (standard R column-major order).
#'
#' @param x A square numeric matrix.
#'
#' @return A numeric vector of length \eqn{n(n-1)/2} containing the strict
#'   lower-triangular entries of \code{x}, where \eqn{n} is the dimension of
#'   \code{x}.
#'
#' @details
#' This is a thin wrapper around \code{x[lower.tri(x)]}.  It is provided as a
#' convenience for extracting vectorized pairwise measurements from symmetric
#' matrices such as distance or covariance matrices.
#'
#' @seealso \code{\link[base]{lower.tri}}
#'
#' @examples
#' m <- matrix(1:9, nrow = 3)
#' lower(m)  # entries (2,1), (3,1), (3,2)
#'
#' # Useful for extracting pairwise distances:
#' D <- matrix(c(0, 0.1, 0.3,
#'               0.1, 0, 0.2,
#'               0.3, 0.2, 0), nrow = 3)
#' lower(D)
#'
#' @export
lower <- function(x)
{
  stopifnot(is.matrix(x), nrow(x) == ncol(x))
  x[lower.tri(x)]
}

#' Pairwise genetic distances from PCA scores
#'
#' Runs a PCA on allele-frequency data from an \code{adegenet::genind} object
#' and returns Euclidean distances computed on the retained principal component
#' scores.
#'
#' @param gi A \code{genind} object from \pkg{adegenet}.
#' @param n_axes Number of principal component axes to retain.  Capped at the
#'   total number of available axes.
#' @param scale Logical.  If \code{TRUE} (default), the resulting distance
#'   matrix is rescaled to \code{[0, 1]} via \code{\link{scale_to_0_1}}.
#'
#' @details
#' Missing genotypes in \code{gi} are imputed by the allele-mean method
#' (equivalent to \code{adegenet::tab(gi, NA.method = "mean")}), then a
#' standard PCA via \code{\link[stats]{prcomp}} is applied to the allele
#' frequency table.  Euclidean distances are then computed on the first
#' \code{n_axes} scores.
#'
#' This is a convenience wrapper designed for users who prefer to work from
#' a \code{genind} object.  For biallelic SNP data the more direct
#' \code{\link{dist_from_biallelic}} is generally faster and does not require
#' the \pkg{adegenet} package.
#'
#' @return A symmetric numeric distance matrix with one row/column per
#'   individual in \code{gi}.
#'
#' @seealso \code{\link{dist_from_biallelic}}, \code{\link{cov_from_genetic_data}}
#'
#' @examples
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   data(nancycats, package = "adegenet")
#'   d <- pca_dist(nancycats, n_axes = 4)
#'   dim(d)
#' }
#'
#' @export
pca_dist <- function(gi, n_axes = 64, scale = TRUE)
{
  a_tab <- .adegenet_tab(gi, NA.method = "mean")
  pc <- prcomp(a_tab)
  n_axes <- min(as.integer(n_axes), ncol(pc$x))
  if (n_axes < 1)
    stop("`n_axes` must retain at least one principal component axis")

  pc_dist <- as.matrix(dist(pc$x[, seq_len(n_axes), drop = FALSE]))
  if (isTRUE(scale))
    pc_dist <- as.matrix(scale_to_0_1(pc_dist))

  pc_dist
}
