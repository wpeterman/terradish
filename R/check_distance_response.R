#' Check a generalized-Wishart distance response
#'
#' Tests whether a matrix supplied as the response to
#' \code{\link{generalized_wishart}} has squared-Euclidean geometry. The check
#' evaluates the positive semidefiniteness of the centered Gram matrix implied
#' by the observed response. It does not assess whether the chosen genetic
#' summary is biologically appropriate for a particular study.
#'
#' @param S Numeric square matrix. Entries are interpreted as the
#'   \strong{squared-distance response} \eqn{D} that would be supplied to
#'   \code{generalized_wishart}, not as ordinary distances \eqn{d}. If you have
#'   ordinary distances, check \code{S^2} and supply that squared matrix to the
#'   model.
#' @param tol Positive numeric scalar. Relative tolerance used for structural
#'   checks and for classifying negative eigenvalues. An eigenvalue is treated
#'   as numerical zero when it is no smaller than
#'   \eqn{-tol \lambda_{max}}, where \eqn{\lambda_{max}} is the largest
#'   eigenvalue of the centered Gram matrix.
#'
#' @details
#' For \eqn{n} sites, the function forms
#'
#' \deqn{B = -\frac{1}{2} H D H, \qquad
#'       H = I - \frac{1}{n} \mathbf{1}\mathbf{1}^{\mathsf T}.}
#'
#' The response is admissible when \eqn{B} is positive semidefinite within the
#' requested tolerance. Symmetry, nonnegative entries, and a zero diagonal are
#' necessary but are not sufficient.
#'
#' This function checks the matrix that you actually intend to analyze because
#' scaling, transformations, locus aggregation, and missing-data handling can
#' change its geometry. It never applies Cailliez, Lingoes, nearest-Euclidean,
#' or other corrections. Those procedures alter the observed response and may
#' alter its biological interpretation.
#'
#' For an independent check with \pkg{ade4}, note that
#' \code{ade4::is.euclid()} expects ordinary distances and squares them
#' internally. For an already squared response \code{S}, the corresponding call
#' is \code{ade4::is.euclid(as.dist(sqrt(S)))}, not
#' \code{ade4::is.euclid(as.dist(S))}.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{admissible}{Logical. \code{TRUE} when the centered Gram matrix is
#'   positive semidefinite within \code{tol}.}
#' \item{min_eigenvalue}{Smallest eigenvalue of the centered Gram matrix, on
#'   the squared-distance scale of \code{S}.}
#' \item{relative_min_eigenvalue}{The smallest eigenvalue divided by the
#'   largest eigenvalue. Values close to zero from below usually reflect
#'   numerical rounding; materially negative values indicate incompatibility.}
#' \item{n_negative}{Number of eigenvalues below the scaled negative
#'   tolerance.}
#' \item{tolerance}{The relative tolerance supplied in \code{tol}.}
#' \item{eigenvalue_threshold}{The absolute negative-eigenvalue magnitude used
#'   for the decision, equal to \code{tol * max_eigenvalue}.}
#' \item{max_eigenvalue}{Largest eigenvalue of the centered Gram matrix.}
#' \item{n_sites}{Number of rows and columns in \code{S}.}
#' }
#'
#' @seealso \code{\link{generalized_wishart}},
#'   \code{\link{wishart_covariance}}, \code{\link{mlpe}}
#'
#' @examples
#' # Squared Euclidean distances from known coordinates are admissible.
#' xy <- matrix(c(0, 0,
#'                1, 0,
#'                0, 1,
#'                1, 1), ncol = 2, byrow = TRUE)
#' D <- as.matrix(dist(xy))^2
#' check_distance_response(D)
#'
#' # Symmetry and a zero diagonal alone are not sufficient.
#' bad <- matrix(c(0, 1, 1,
#'                 1, 0, 9,
#'                 1, 9, 0), nrow = 3, byrow = TRUE)
#' check_distance_response(bad)$admissible
#'
#' @export
check_distance_response <- function(S, tol = 1e-7)
{
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("`tol` must be one positive, finite number.", call. = FALSE)
  if (!is.matrix(S))
    stop("`S` must be a numeric square matrix containing squared distances.",
         call. = FALSE)
  if (!is.numeric(S))
    stop("`S` must be numeric.", call. = FALSE)
  if (length(dim(S)) != 2L || nrow(S) != ncol(S))
    stop("`S` must be square.", call. = FALSE)
  if (nrow(S) < 2L)
    stop("`S` must contain at least two sites.", call. = FALSE)
  if (any(!is.finite(S)))
    stop("`S` must contain only finite values.", call. = FALSE)

  response_scale <- max(abs(S))
  structural_scale <- max(response_scale, .Machine$double.eps)
  structural_threshold <- tol * structural_scale

  max_asymmetry <- max(abs(S - t(S)))
  if (max_asymmetry > structural_threshold)
    stop("`S` must be symmetric within the requested tolerance; maximum ",
         "asymmetry is ", format(max_asymmetry, digits = 5), ".",
         call. = FALSE)

  max_diagonal <- max(abs(diag(S)))
  if (max_diagonal > structural_threshold)
    stop("`S` must have a zero diagonal within the requested tolerance; ",
         "maximum absolute diagonal entry is ",
         format(max_diagonal, digits = 5), ".", call. = FALSE)

  min_entry <- min(S)
  if (min_entry < -structural_threshold)
    stop("`S` must be nonnegative within the requested tolerance; minimum ",
         "entry is ", format(min_entry, digits = 5), ".", call. = FALSE)

  # Remove only numerical-scale structural deviations before the geometry test.
  D <- (S + t(S)) / 2
  diag(D) <- 0
  D[D < 0] <- 0

  n <- nrow(D)
  H <- diag(n) - matrix(1 / n, nrow = n, ncol = n)
  B <- -0.5 * H %*% D %*% H
  B <- (B + t(B)) / 2
  eigenvalues <- eigen(B, symmetric = TRUE, only.values = TRUE)$values

  max_eigenvalue <- max(eigenvalues)
  eigenvalue_threshold <- tol * max(max_eigenvalue, 0)
  min_eigenvalue <- min(eigenvalues)
  relative_min_eigenvalue <- if (max_eigenvalue > 0)
    min_eigenvalue / max_eigenvalue
  else
    0
  n_negative <- sum(eigenvalues < -eigenvalue_threshold)

  list(
    admissible = n_negative == 0L,
    min_eigenvalue = min_eigenvalue,
    relative_min_eigenvalue = relative_min_eigenvalue,
    n_negative = n_negative,
    tolerance = tol,
    eigenvalue_threshold = eigenvalue_threshold,
    max_eigenvalue = max_eigenvalue,
    n_sites = n
  )
}

.prepare_gw_response <- function(S, tol = 1e-7)
{
  checked <- attr(S, "terradish_gw_response_check", exact = TRUE)
  if (is.list(checked) && identical(checked$tolerance, tol) &&
      isTRUE(checked$admissible))
    return(S)

  checked <- check_distance_response(S, tol = tol)
  if (!isTRUE(checked$admissible))
  {
    stop(
      "The generalized-Wishart response is not squared Euclidean. Its centered ",
      "Gram matrix has ", checked$n_negative, " eigenvalue(s) below the ",
      "scaled tolerance; the minimum relative eigenvalue is ",
      format(checked$relative_min_eigenvalue, digits = 5), ". Use ",
      "`check_distance_response(S)` for diagnostics. Do not silently correct ",
      "the matrix; choose a scientifically appropriate admissible response, ",
      "or use `mlpe` for a distance response that does not meet this requirement.",
      call. = FALSE
    )
  }

  response_scale <- max(abs(S))
  structural_threshold <- tol * max(response_scale, .Machine$double.eps)
  S <- (S + t(S)) / 2
  diag(S) <- 0
  S[S < 0 & S >= -structural_threshold] <- 0
  attr(S, "terradish_gw_response_check") <- checked
  S
}

.validate_measurement_response <- function(measurement_model, S, tol = 1e-7)
{
  model_name <- attr(measurement_model, "base_model", exact = TRUE)
  if (is.null(model_name) &&
      identical(formals(measurement_model), formals(generalized_wishart)) &&
      identical(body(measurement_model), body(generalized_wishart)))
    model_name <- "generalized_wishart"
  if (identical(model_name, "generalized_wishart"))
    return(.prepare_gw_response(S, tol = tol))
  S
}
