.pairwise_site_count <- function(pairs) {
  n <- (1 + sqrt(1 + 8 * pairs)) / 2
  if (length(n) != 1L || !is.finite(n) || n < 2 || abs(n - round(n)) > 1e-8)
    stop("A pairwise matrix must have n(n-1)/2 rows. For site-level covariates, use a data frame or pairwise_endpoint_covariates().", call. = FALSE)
  as.integer(round(n))
}

#' Combine environmental and other pairwise covariates
#'
#' Combines differences calculated from site attributes with user-supplied
#' pairwise matrices, such as geographic distance, in one measurement-model
#' input. Values retain their supplied units and scaling.
#' @param ... Objects from \code{pairwise_endpoint_covariates()}, or named
#'   symmetric numeric matrices or \code{dist} objects. Matrices must have
#'   a zero diagonal. Every input must describe the same sites in the same order.
#' @return A classed numeric matrix with one row per unordered site pair and
#'   one column per covariate. Rows follow the lower-triangle order used by R's
#'   \code{dist()}. Both \code{mlpe_covariates()} and
#'   \code{wishart_covariates()} accept this object.
#' @details No centering, scaling, or geometric correction is applied. For a
#'   Wishart model, each supplied distance must additionally produce a positive
#'   semidefinite centered kernel; the measurement-model constructor checks it.
#'   Site labels, when supplied on square matrices, must agree across inputs.
#' @examples
#' climate <- pairwise_endpoint_covariates(data.frame(temp = c(4, 7, 9)))
#' geography <- dist(cbind(c(0, 1, 3), c(0, 2, 1)))
#' z <- pairwise_covariates(climate, geographic = geography)
#' colnames(z)
#' @seealso \code{\link{pairwise_endpoint_covariates}},
#'   \code{\link{mlpe_covariates}}, \code{\link{wishart_covariates}}
#' @export
#' @template pairwise-interpretation
pairwise_covariates <- function(...) {
  inputs <- list(...)
  if (!length(inputs)) stop("Supply at least one pairwise covariate.", call. = FALSE)
  input_names <- names(inputs)
  if (is.null(input_names)) input_names <- rep("", length(inputs))
  n_sites <- NULL
  site_names <- NULL
  columns <- lapply(seq_along(inputs), function(i) {
    x <- inputs[[i]]
    name <- input_names[i]
    if (inherits(x, c("terradish_pairwise_covariates", "radish_pairwise_covariates"))) {
      out <- .as_pairwise_covariate_matrix(x)
      n <- .pairwise_site_count(nrow(out))
    } else {
      if (!nzchar(name)) stop("Name each supplied square matrix or dist object.", call. = FALSE)
      if (!is.matrix(x) && !inherits(x, "dist"))
        stop("Inputs must be classed pairwise covariates, square matrices, or dist objects.", call. = FALSE)
      x <- as.matrix(x)
      if (!is.numeric(x) || nrow(x) < 2L || nrow(x) != ncol(x) ||
          any(!is.finite(x)) || !isTRUE(all.equal(x, t(x), check.attributes = FALSE)) ||
          any(abs(diag(x)) > 1e-10))
        stop("Pairwise matrices must be finite, square, symmetric, and have a zero diagonal.", call. = FALSE)
      n <- nrow(x)
      labels <- rownames(x)
      if (!is.null(labels) && !is.null(colnames(x)) && !identical(labels, colnames(x)))
        stop("Row and column site labels must agree.", call. = FALSE)
      if (!is.null(labels)) {
        if (!is.null(site_names) && !identical(labels, site_names))
          stop("Site labels or their order differ across inputs.", call. = FALSE)
        site_names <<- labels
      }
      out <- matrix(x[lower.tri(x)], ncol = 1L)
    }
    if (!is.null(n_sites) && n != n_sites)
      stop("All inputs must describe the same number of sites.", call. = FALSE)
    n_sites <<- n
    if (nzchar(name)) colnames(out) <- if (ncol(out) == 1L) name else
      paste(name, colnames(out), sep = ".")
    out
  })
  out <- do.call(cbind, columns)
  if (any(!is.finite(out)) || is.null(colnames(out)) ||
      any(!nzchar(colnames(out))) || anyDuplicated(colnames(out)))
    stop("Covariates must have finite values and unique, nonempty names.", call. = FALSE)
  structure(out, n_sites = n_sites, site_names = site_names,
            class = c("terradish_pairwise_covariates", "matrix", "array"))
}
