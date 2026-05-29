#' F_ST from biallelic allele counts
#'
#' Estimates pairwise F\eqn{_{ST}} from counts of the derived allele across
#' biallelic markers using the ratio-of-averages estimator of Bhatia et al.
#' (2013).
#'
#' @param Y Numeric matrix of derived-allele counts with rows as populations
#'   and columns as loci.
#' @param N Numeric matrix of sampled haploid chromosomes with the same
#'   dimensions as \code{Y}.  For diploid data without missing genotypes,
#'   \code{N = 2 * n_individuals} for each cell.
#'
#' @details
#' The estimator computes, for each locus, the squared allele-frequency
#' difference between populations, corrects for within-population
#' heterozygosity, and then forms a ratio of locus-averages (Bhatia et al.
#' 2013).  Diagonal entries are set to 0.
#'
#' The result is not guaranteed to lie in \eqn{[0, 1]} for every pair; small
#' negative values can arise from sampling noise in very similar populations.
#'
#' To use F\eqn{_{ST}} as the response in \code{\link{terradish}}, pass the
#' matrix directly to the formula left-hand side with \code{\link{mlpe}} or
#' \code{\link{generalized_wishart}}.
#'
#' @references
#' Bhatia G, Patterson N, Sankararaman S, Price AL. 2013. Estimating and
#' interpreting F\eqn{_{ST}}: the impact of rare variants. Genome Research
#' 23(9):1514-1521.
#'
#' @return A symmetric numeric matrix of pairwise F\eqn{_{ST}} with zero
#'   diagonal.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{dist_from_biallelic}}
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' fst_from_biallelic(Y, N)
#'
#' @export

fst_from_biallelic <- function(Y, N)
{
  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 1))
    stop("All sample sizes in `N` must exceed 1")

  Fr <- Y/N
  f2 <- apply(apply(Fr, 2, function(x) outer(x,x,"-")^2), 1, mean, na.rm=TRUE)
  cornum <- apply(Fr*(1-Fr)/(N-1), 1, mean, na.rm=TRUE)
  corden <- apply(Fr*(1-Fr)*N/(N-1), 1, mean, na.rm=TRUE)
  fst <- (f2 - c(outer(cornum, cornum, "+")))/(f2 - c(outer(cornum, cornum, "+")) + c(outer(corden, corden, "+")))
  fst <- matrix(fst, nrow(Y), nrow(Y))
  diag(fst) <- 0
  fst
}

#' Genetic covariance from biallelic allele counts
#'
#' Computes a sample covariance matrix from counts of the derived allele across
#' biallelic (SNP) markers.  Rows may be individuals, populations, or any other
#' sampled genetic unit.  The result is directly suitable as the response
#' \code{S} in \code{\link{wishart_covariance}}.
#'
#' @param Y Numeric matrix of derived-allele counts with sampled units in rows
#'   and loci in columns.  For diploid individuals genotyped as 0/1/2, this is
#'   the standard dosage matrix.
#' @param N Optional haploid sample-size specification.  Controls how many
#'   chromosomes were sampled at each locus for each unit:
#'   \itemize{
#'     \item \code{NULL}: use \code{ploidy} for all cells.
#'     \item A scalar: recycled to all cells.
#'     \item A vector of length \code{nrow(Y)}: row-specific sample sizes
#'       (same at every locus for that row).
#'     \item A vector of length \code{ncol(Y)}: locus-specific sample sizes
#'       (same across all rows).
#'     \item A matrix matching \code{dim(Y)}: fully general cell-level sizes.
#'   }
#' @param ploidy Haploid chromosome count per unit, used when \code{N = NULL}.
#'   The default \code{2} corresponds to diploid individuals with dosage coding
#'   0/1/2 (so each individual contributes two chromosomes at each locus).
#' @param monomorphic How to handle loci with pooled allele frequency 0 or 1
#'   (fixed across all rows).  \code{"drop"} (default) removes them and emits a
#'   warning; \code{"error"} stops with an error to preserve historical strict
#'   behavior.
#' @param tol Tolerance for identifying monomorphic loci.  A locus is treated
#'   as monomorphic when its pooled frequency is within \code{tol} of 0 or 1.
#'
#' @details
#' Let \eqn{p_\ell = \sum_i Y_{i\ell} / \sum_i N_{i\ell}} be the pooled
#' derived-allele frequency at locus \eqn{\ell}.  The normalized allele dosage
#' is:
#'
#' \deqn{Z_{i\ell} = \frac{Y_{i\ell} - N_{i\ell} p_\ell}
#'                        {\sqrt{N_{i\ell} p_\ell (1 - p_\ell)}}}
#'
#' The function returns \eqn{Z Z^\top / L}, where \eqn{L} is the number of
#' retained (polymorphic) loci.  The resulting matrix is the sample covariance
#' in normalized allele-frequency space, analogous to the genomic relationship
#' matrix of Yang et al. (2010).
#'
#' To convert the covariance to a pairwise distance matrix, use
#' \code{\link{dist_from_cov}}.  For non-biallelic or multiallelic data, use
#' \code{\link{cov_from_genetic_data}}.
#'
#' @references
#' Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR,
#' Madden PA, Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher PM.
#' 2010. Common SNPs explain a large proportion of the heritability for human
#' height. Nature Genetics 42(7):565-569.
#'
#' @return A symmetric positive-semidefinite numeric matrix of pairwise
#'   genetic covariances.  Row and column names match those of \code{Y} when
#'   present.
#'
#' @seealso \code{\link{cov_from_genetic_data}}, \code{\link{fst_from_biallelic}},
#'   \code{\link{dist_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{\link{wishart_covariance}}
#'
#' @examples
#' # Individual diploid dosage matrix: 3 individuals, 3 SNPs
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' cov_from_biallelic(Y)          # ploidy = 2 by default
#'
#' # Population-level counts with unequal sampling
#' Y2 <- matrix(c(10, 4, 1,
#'                5,  5, 3), nrow = 2, byrow = TRUE)
#' N2 <- matrix(c(20, 20, 10,
#'                10, 10, 6), nrow = 2, byrow = TRUE)
#' cov_from_biallelic(Y2, N = N2)
#'
#' @export

cov_from_biallelic <- function(Y,
                               N = NULL,
                               ploidy = 2,
                               monomorphic = c("drop", "error"),
                               tol = sqrt(.Machine$double.eps))
{
  monomorphic <- match.arg(monomorphic)

  if (is.data.frame(Y))
    Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L)
    stop("`Y` must be a numeric matrix or data frame with at least two rows and one column",
         call. = FALSE)
  storage.mode(Y) <- "double"

  N <- .biallelic_sample_size_matrix(N, Y, ploidy)

  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(!is.finite(Y)) || any(!is.finite(N)))
    stop("`Y` and `N` must contain only finite values", call. = FALSE)
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 0))
    stop("All sample sizes in `N` must be positive")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a finite nonnegative number", call. = FALSE)

  Fr <- colSums(Y) / colSums(N)
  variable <- is.finite(Fr) & Fr > tol & Fr < (1 - tol)
  if (!all(variable))
  {
    if (identical(monomorphic, "error"))
      stop("Each retained locus must have pooled allele frequency strictly between 0 and 1",
           call. = FALSE)
    if (!any(variable))
      stop("No variable biallelic loci remain after dropping monomorphic loci",
           call. = FALSE)
    dropped <- colnames(Y)[!variable]
    if (is.null(dropped))
      dropped <- which(!variable)
    warning("Dropping monomorphic biallelic loci: ",
            paste(dropped, collapse = ", "),
            call. = FALSE)
    Y <- Y[, variable, drop = FALSE]
    N <- N[, variable, drop = FALSE]
    Fr <- Fr[variable]
  }

  Y  <- (Y - N*Fr) / sqrt(N * Fr * (1-Fr))

  Y %*% t(Y) / ncol(Y)
}

.biallelic_sample_size_matrix <- function(N, Y, ploidy)
{
  if (is.null(N))
  {
    if (!is.numeric(ploidy) || length(ploidy) != 1L ||
        !is.finite(ploidy) || ploidy <= 0)
      stop("`ploidy` must be one finite positive number", call. = FALSE)
    N <- ploidy
  }

  if (is.data.frame(N))
    N <- as.matrix(N)
  if (is.matrix(N))
  {
    if (!all(dim(Y) == dim(N)))
      stop("Dimension mismatch", call. = FALSE)
    if (!is.numeric(N))
      stop("`N` must be numeric", call. = FALSE)
    storage.mode(N) <- "double"
    return(N)
  }
  if (!is.numeric(N))
    stop("`N` must be numeric", call. = FALSE)
  if (length(N) == 1L)
    return(matrix(N, nrow(Y), ncol(Y), dimnames = dimnames(Y)))
  if (length(N) == nrow(Y))
    return(matrix(N, nrow(Y), ncol(Y), dimnames = dimnames(Y)))
  if (length(N) == ncol(Y))
    return(matrix(N, nrow(Y), ncol(Y), byrow = TRUE, dimnames = dimnames(Y)))

  stop("`N` must be a scalar, a matrix matching `Y`, or a row- or locus-length vector",
       call. = FALSE)
}

#' Covariance from multivariate genetic data
#'
#' Calculates individual- or population-level genetic covariance from
#' multivariate genetic data. This is useful for microsatellite, multiallelic,
#' SNP dosage, PCA score, or other numeric encodings.
#'
#' @param x Genetic data. For \code{input = "features"}, a numeric matrix or data
#'   frame with individuals in rows and genetic features in columns. For
#'   \code{input = "allele_calls"}, a matrix or data frame of allele calls with
#'   one column per allele copy.
#' @param groups Optional factor, character, or integer vector assigning each row
#'   of \code{x} to a population. If \code{NULL}, each row is treated as an
#'   individual sampled unit and an individual-level covariance matrix is
#'   returned.
#' @param input Input format. \code{"features"} treats \code{x} as an already
#'   numeric feature matrix. \code{"allele_calls"} converts allele calls to
#'   per-locus allele dosage columns before calculating covariance.
#' @param loci Required for \code{input = "allele_calls"}; a vector with one
#'   entry per column of \code{x} identifying the locus for each allele-copy
#'   column. For a diploid microsatellite matrix with two columns per locus, the
#'   two columns for each locus should have the same \code{loci} value.
#' @param center Should feature columns be centered by their global mean before
#'   covariance calculation? Default \code{TRUE}.
#' @param scale Should feature columns be divided by their global standard
#'   deviation? Default \code{TRUE}. Constant features are dropped.
#' @param tol Tolerance used to identify constant features when
#'   \code{scale = TRUE}.
#' @param diagonal How to set the returned covariance diagonal. \code{"auto"}
#'   uses \code{"within"} for grouped population data with replication and
#'   \code{"gower"} for individual-level data. \code{"within"} replaces the
#'   Gower covariance diagonal with the within-population genetic variance,
#'   following the population-graph covariance construction. \code{"gower"}
#'   leaves the double-centered distance diagonal unchanged.
#' @param normalize Should covariance scale be normalized? \code{"none"} returns
#'   sums over retained features. \code{"features"} divides the covariance,
#'   squared distances, and within-population variances by the number of retained
#'   features.
#'
#' @details
#' The function first obtains a numeric feature matrix. Microsatellite data can
#' be supplied as allele calls by setting \code{input = "allele_calls"}; these are
#' converted to allele dosage columns, one column per locus-allele combination.
#' Missing allele calls are imputed to the modal observed allele within each
#' locus, with a message reporting the number of imputed calls.
#'
#' For Wishart models, the \code{nu} argument is the effective degrees of
#' freedom in the empirical covariance, not necessarily the biological locus
#' count.  With biallelic SNPs, each retained SNP contributes one independent
#' frequency column, so \code{nu} is usually the retained SNP count.  With
#' microsatellites, a locus with \eqn{K_l} observed alleles contributes
#' approximately \eqn{K_l - 1} independent allele frequencies because allele
#' frequencies sum to a constant within the locus.  A useful microsatellite
#' choice is therefore \eqn{nu = \sum_l (K_l - 1)}.  Markers in strong linkage
#' disequilibrium are not independent, so the effective \code{nu} is smaller
#' than the raw count; thinning to roughly unlinked loci gives a more honest
#' value.  Because \code{nu} controls the precision of every downstream Wishart
#' fit (it scales standard errors and model-selection statistics but leaves
#' point estimates unchanged), see \code{\link{wishart_covariance}} for the
#' full explanation before settling on a value.
#'
#' If \code{groups = NULL}, rows of the processed feature matrix are used
#' directly and squared Euclidean distances among individuals are transformed to
#' a covariance matrix by Gower double-centering. This produces a positive
#' semidefinite row-level covariance matrix up to numerical tolerance.
#'
#' If \code{groups} is supplied, population centroids are calculated in feature
#' space before the Gower transform. With \code{diagonal = "within"}, the
#' diagonal is replaced by the sum of within-population feature variances. This
#' matches the covariance construction used before population-graph
#' partial-correlation filtering in Dyer-style population graph workflows.
#'
#' Missing values are not currently supported because common microsatellite
#' imputation choices can change the resulting covariance.
#'
#' If the grouped Dyer-style result is used as \code{S} in
#' \code{\link{wishart_covariance}}, inspect the eigenvalues first. The
#' within-population diagonal is useful for population-graph workflows, but it
#' does not guarantee a positive-definite covariance matrix for every dataset.
#'
#' @return A symmetric numeric matrix with one row and column per individual
#'   (when \code{groups = NULL}) or population (when \code{groups} is
#'   supplied).  The matrix is positive semi-definite up to numerical
#'   tolerance, except when \code{diagonal = "within"} is used (see Details).
#'   The following attributes are attached:
#' \describe{
#'   \item{\code{centroids}}{Feature matrix of population centroids (or
#'     individual feature vectors when ungrouped).}
#'   \item{\code{centroid_distance2}}{Squared Euclidean distance matrix among
#'     centroids.}
#'   \item{\code{within_variance}}{Named vector of within-population variances
#'     (only when \code{groups} is supplied and some populations have \eqn{\geq
#'     2} individuals).}
#'   \item{\code{unit_size}}{Named integer vector of group sizes.}
#'   \item{\code{feature_center}, \code{feature_scale}}{Centering and scaling
#'     constants applied to each retained feature.}
#'   \item{\code{retained_features}}{Names of the features kept after constant
#'     features were dropped.}
#'   \item{\code{level}}{Either \code{"individual"} or \code{"population"}.}
#'   \item{\code{diagonal}}{The \code{diagonal} method actually used.}
#'   \item{\code{normalizer}}{Divisor applied to the covariance and distances
#'     (\code{1} for \code{normalize = "none"}, number of retained features for
#'     \code{normalize = "features"}).}
#' }
#'
#' @seealso \code{\link{wishart_covariance}}, \code{\link{cov_from_biallelic}},
#'   \code{\link{dist_from_cov}}, \code{\link{simulate_covariance_response}}
#'
#' @examples
#' # Numeric allele dosage / feature matrix
#' x <- matrix(c(0, 1,
#'               1, 1,
#'               2, 0,
#'               2, 1,
#'               0, 2,
#'               1, 2),
#'             ncol = 2, byrow = TRUE)
#' cov_from_genetic_data(x)
#'
#' groups <- rep(c("pop1", "pop2", "pop3"), each = 2)
#' cov_from_genetic_data(x, groups = groups)
#'
#' # Microsatellite-style allele-call columns: two allele copies per locus
#' alleles <- data.frame(
#'   loc1_a = c(100, 100, 102, 102, 104, 104),
#'   loc1_b = c(100, 102, 102, 104, 104, 100),
#'   loc2_a = c(200, 202, 200, 202, 204, 204),
#'   loc2_b = c(202, 202, 204, 204, 204, 200)
#' )
#' cov_from_genetic_data(
#'   alleles,
#'   groups = groups,
#'   input = "allele_calls",
#'   loci = c("loc1", "loc1", "loc2", "loc2")
#' )
#'
#' @export
cov_from_genetic_data <- function(x,
                                  groups = NULL,
                                  input = c("features", "allele_calls"),
                                  loci = NULL,
                                  center = TRUE,
                                  scale = TRUE,
                                  tol = sqrt(.Machine$double.eps),
                                  diagonal = c("auto", "within", "gower"),
                                  normalize = c("none", "features"))
{
  input <- match.arg(input)
  diagonal <- match.arg(diagonal)
  normalize <- match.arg(normalize)

  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x) || nrow(x) < 2L || ncol(x) < 1L)
    stop("`x` must be a matrix or data frame with at least two rows and one column",
         call. = FALSE)
  if (!is.null(groups))
  {
    if (length(groups) != nrow(x))
      stop("`groups` must have one entry per row of `x`", call. = FALSE)
    if (anyNA(groups))
      stop("missing values are not supported in `groups`", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center))
    stop("`center` must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
    stop("`scale` must be TRUE or FALSE", call. = FALSE)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a finite nonnegative number", call. = FALSE)

  features <- if (identical(input, "allele_calls"))
    .allele_call_feature_matrix(x, loci)
  else
    .numeric_genetic_feature_matrix(x)
  imputed_allele_calls <- attr(features, "imputed_allele_calls")
  imputed_modal_alleles <- attr(features, "imputed_modal_alleles")

  feature_center <- if (isTRUE(center))
    colMeans(features)
  else
    setNames(rep(0, ncol(features)), colnames(features))
  features <- sweep(features, 2L, feature_center, "-")

  feature_scale <- setNames(rep(1, ncol(features)), colnames(features))
  retained <- rep(TRUE, ncol(features))
  if (isTRUE(scale))
  {
    feature_scale <- apply(features, 2L, stats::sd)
    retained <- is.finite(feature_scale) & feature_scale > tol
    if (!any(retained))
      stop("No variable genetic features remain after scaling", call. = FALSE)
    if (any(!retained))
      warning("Dropping constant or near-constant genetic features: ",
              paste(colnames(features)[!retained], collapse = ", "),
              call. = FALSE)
    features <- sweep(features[, retained, drop = FALSE],
                      2L, feature_scale[retained], "/")
    feature_center <- feature_center[retained]
    feature_scale <- feature_scale[retained]
  }

  normalizer <- if (identical(normalize, "features")) ncol(features) else 1

  if (is.null(groups))
  {
    unit_features <- features
    unit_names <- rownames(features)
    if (is.null(unit_names) || anyNA(unit_names) || any(!nzchar(unit_names)))
      unit_names <- paste0("sample", seq_len(nrow(features)))
    else if (anyDuplicated(unit_names))
      unit_names <- make.unique(unit_names)
    rownames(unit_features) <- unit_names
    unit_size <- setNames(rep(1L, nrow(unit_features)), unit_names)
    within_variance <- NULL
    level <- "individual"
    if (identical(diagonal, "within"))
      stop("`diagonal = \"within\"` requires population groups; use `diagonal = \"gower\"` for individual-level covariance",
           call. = FALSE)
    diagonal <- "gower"
  }
  else
  {
    groups <- if (is.factor(groups))
      droplevels(groups)
    else
      factor(groups, levels = unique(groups))
    if (nlevels(groups) < 2L)
      stop("`groups` must contain at least two populations", call. = FALSE)
    group_size <- as.numeric(table(groups))
    names(group_size) <- levels(groups)

    unit_features <- rowsum(features, groups, reorder = TRUE)
    unit_features <- sweep(unit_features, 1L, group_size, "/")
    rownames(unit_features) <- levels(groups)
    unit_size <- group_size
    level <- if (all(group_size == 1L)) "individual" else "population"

    if (identical(diagonal, "auto"))
      diagonal <- if (identical(level, "individual")) "gower" else "within"
    if (identical(diagonal, "within"))
    {
      if (all(group_size < 2L))
        stop("`diagonal = \"within\"` requires at least one population with two or more samples",
             call. = FALSE)
      if (any(group_size < 2L))
        warning("Some populations contain fewer than two individuals; their within-population variance is set to zero.",
                call. = FALSE)
    }

    within_variance <- vapply(levels(groups), function(g)
    {
      xg <- features[groups == g, , drop = FALSE]
      if (nrow(xg) < 2L)
        return(0)
      sum(apply(xg, 2L, stats::var))
    }, numeric(1))
  }

  centroid_dist2 <- as.matrix(stats::dist(unit_features))^2 / normalizer
  row_mean <- rowMeans(centroid_dist2)
  col_mean <- colMeans(centroid_dist2)
  grand_mean <- mean(centroid_dist2)
  covariance <- -0.5 * (
    centroid_dist2 -
      outer(row_mean, rep(1, length(row_mean))) -
      outer(rep(1, length(col_mean)), col_mean) +
      grand_mean
  )

  if (!is.null(within_variance))
    within_variance <- within_variance / normalizer
  if (identical(diagonal, "within"))
    diag(covariance) <- within_variance

  covariance <- (covariance + t(covariance)) / 2
  rownames(covariance) <- colnames(covariance) <- rownames(unit_features)
  attr(covariance, "centroids") <- unit_features
  attr(covariance, "unit_features") <- unit_features
  if (!is.null(within_variance))
    attr(covariance, "within_variance") <- within_variance
  attr(covariance, "centroid_distance2") <- centroid_dist2
  attr(covariance, "unit_distance2") <- centroid_dist2
  attr(covariance, "unit_size") <- unit_size
  attr(covariance, "feature_center") <- feature_center
  attr(covariance, "feature_scale") <- feature_scale
  attr(covariance, "retained_features") <- colnames(features)
  attr(covariance, "input") <- input
  if (!is.null(imputed_allele_calls))
  {
    attr(covariance, "imputed_allele_calls") <- imputed_allele_calls
    attr(covariance, "imputed_modal_alleles") <- imputed_modal_alleles
  }
  attr(covariance, "level") <- level
  attr(covariance, "diagonal") <- diagonal
  attr(covariance, "normalize") <- normalize
  attr(covariance, "normalizer") <- normalizer

  covariance
}

.numeric_genetic_feature_matrix <- function(x)
{
  if (!is.numeric(x))
    stop("`x` must be numeric when `input = \"features\"`", call. = FALSE)
  storage.mode(x) <- "double"
  if (anyNA(x))
    stop("missing values are not currently supported in `x`", call. = FALSE)
  if (is.null(colnames(x)))
    colnames(x) <- paste0("feature", seq_len(ncol(x)))
  x
}

.allele_call_feature_matrix <- function(x, loci)
{
  if (is.null(loci))
    stop("`loci` is required when `input = \"allele_calls\"`", call. = FALSE)
  if (length(loci) != ncol(x))
    stop("`loci` must have one entry per allele-call column", call. = FALSE)
  if (anyNA(loci) || any(!nzchar(as.character(loci))))
    stop("`loci` values must be non-missing and non-empty", call. = FALSE)

  x_chr <- matrix(as.character(x), nrow = nrow(x), ncol = ncol(x))
  loci <- as.character(loci)
  locus_names <- unique(loci)
  out <- vector("list", length(locus_names))
  names(out) <- locus_names
  imputed <- integer(length(locus_names))
  names(imputed) <- locus_names
  imputed_allele <- setNames(rep(NA_character_, length(locus_names)),
                             locus_names)

  for (loc in locus_names)
  {
    loc_values <- x_chr[, loci == loc, drop = FALSE]
    missing <- is.na(loc_values)
    if (any(missing))
    {
      observed <- loc_values[!missing]
      if (!length(observed))
        stop("Cannot impute allele calls for locus `", loc,
             "` because all calls are missing", call. = FALSE)
      modal <- .modal_allele(observed)
      loc_values[missing] <- modal
      imputed[[loc]] <- sum(missing)
      imputed_allele[[loc]] <- modal
    }
    alleles <- sort(unique(c(loc_values)))
    loc_counts <- vapply(alleles, function(allele)
      rowSums(loc_values == allele),
      numeric(nrow(loc_values)))
    loc_counts <- as.matrix(loc_counts)
    colnames(loc_counts) <- paste0(loc, ":", alleles)
    out[[loc]] <- loc_counts
  }

  if (sum(imputed) > 0L)
    message("Imputed ", sum(imputed),
            " missing allele call(s) to the modal allele within each locus.")

  features <- do.call(cbind, out)
  attr(features, "imputed_allele_calls") <- imputed[imputed > 0L]
  attr(features, "imputed_modal_alleles") <- imputed_allele[imputed > 0L]
  features
}

.modal_allele <- function(x)
{
  counts <- table(x)
  modes <- names(counts)[counts == max(counts)]
  sort(modes)[1]
}

#' Squared-distance matrix from a covariance matrix
#'
#' Converts a covariance matrix to its associated squared pairwise distance
#' matrix.  This inverts the Gower double-centering used by
#' \code{\link{cov_from_genetic_data}}.
#'
#' @param Cov A square numeric covariance matrix.  Does not need to be
#'   full-rank (e.g. the output of \code{\link{cov_from_biallelic}} is
#'   positive semi-definite).
#'
#' @details
#' The squared distance between units \eqn{i} and \eqn{j} is:
#'
#' \deqn{D_{ij} = C_{ii} + C_{jj} - 2 C_{ij}}
#'
#' which is the law-of-cosines identity relating a covariance matrix to its
#' implied squared Euclidean distances.  The diagonal of the returned matrix
#' is zero.
#'
#' Use this function to convert a covariance matrix produced by
#' \code{\link{cov_from_biallelic}} or \code{\link{cov_from_genetic_data}} into
#' a distance matrix suitable for \code{\link{generalized_wishart}} or
#' \code{\link{mlpe}}.
#'
#' @return A symmetric numeric matrix of the same dimension as \code{Cov} with
#'   zero diagonal and non-negative off-diagonal entries.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{cov_from_genetic_data}},
#'   \code{\link{dist_from_biallelic}}
#'
#' @examples
#' Cov <- matrix(c(2.0, 1.2, 0.8,
#'                 1.2, 1.5, 0.7,
#'                 0.8, 0.7, 1.1), nrow = 3, byrow = TRUE)
#' dist_from_cov(Cov)
#'
#' # Round-trip: covariance -> distance -> (verify non-negative, zero diag)
#' Y <- matrix(c(2, 1, 0, 1, 1, 1, 0, 1, 2), nrow = 3, byrow = TRUE)
#' S_cov  <- cov_from_biallelic(Y)
#' S_dist <- dist_from_cov(S_cov)
#' diag(S_dist)   # all zeros
#'
#' @export

dist_from_cov <- function(Cov)
{
  ones <- matrix(1, nrow(Cov), 1)
  diag(Cov) %*% t(ones) + ones %*% t(diag(Cov)) - 2*Cov
}

#' Genetic distance matrix from biallelic allele counts
#'
#' Computes pairwise squared genetic distances from biallelic (SNP) allele
#' counts.  This is a convenience wrapper for
#' \code{\link{dist_from_cov}(\link{cov_from_biallelic}(Y, N))}.
#'
#' @param Y Numeric matrix of derived-allele counts (populations/individuals in
#'   rows, loci in columns).
#' @param N Numeric matrix of haploid sample sizes with the same dimensions as
#'   \code{Y}.  See \code{\link{cov_from_biallelic}} for flexible \code{N}
#'   specifications.
#'
#' @return A symmetric numeric matrix of pairwise squared distances with zero
#'   diagonal.  Suitable for use with \code{\link{mlpe}} or
#'   \code{\link{generalized_wishart}}.
#'
#' @seealso \code{\link{cov_from_biallelic}}, \code{\link{dist_from_cov}},
#'   \code{\link{fst_from_biallelic}}
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' dist_from_biallelic(Y, N)
#'
#' @export

dist_from_biallelic <- function(Y, N)
{
  dist_from_cov(cov_from_biallelic(Y, N))
}

#' Simulate covariance responses from a conductance surface
#'
#' Simulates one or more covariance response matrices from a known
#' conductance-resistance relationship, using the same Wishart covariance model
#' fitted by \code{\link{wishart_covariance}}.
#'
#' @param theta Conductance parameters. May be supplied as a numeric vector or
#'   a one-row matrix. If unnamed, values are matched in the order implied by
#'   \code{formula}.
#' @param formula Model formula specifying the conductance covariates. The left
#'   hand side, if supplied, is ignored.
#' @param data A \code{\link{conductance_surface}} object.
#' @param conductance_model Conductance-model factory, such as
#'   \code{\link{loglinear_conductance}}.
#' @param tau Nonnegative scaling applied to the conductance-implied covariance
#'   matrix.
#' @param sigma Nonnegative nugget variance added to the diagonal.
#' @param nu Effective Wishart degrees of freedom: the number of
#'   (approximately independent) genetic markers the simulated \code{S} should
#'   represent.  This controls the amount of sampling noise in the simulated
#'   covariance: larger \code{nu} gives tighter draws around \code{Sigma}.  For
#'   biallelic SNP-like simulations this is usually the number of retained
#'   SNPs; for microsatellite-like panels, use approximately
#'   \eqn{\sum_l (K_l - 1)} where \eqn{K_l} is the number of observed alleles at
#'   locus \eqn{l}.  See \code{\link{wishart_covariance}} for how \code{nu}
#'   propagates into estimation (it scales standard errors and model-selection
#'   statistics but not point estimates).
#' @param nsim Number of covariance matrices to simulate.
#' @param seed Optional random seed.
#' @param cores Number of cores passed to \code{\link{terradish_distance}}.
#'
#' @details
#' The helper first computes the conductance-implied covariance matrix
#' \code{E(theta)} from the supplied resistance surface. It then forms
#' \code{Sigma = tau * E(theta) + sigma * I} and simulates
#' \code{S ~ Wishart(nu, Sigma) / nu}. This makes the returned covariance
#' matrices directly compatible with \code{\link{wishart_covariance}}.
#'
#' @return A list containing:
#' \item{covariance}{A covariance matrix if \code{nsim = 1}, otherwise a
#'   three-dimensional array of simulated covariance matrices.}
#' \item{Sigma}{The covariance matrix used as the Wishart scale matrix after
#'   applying \code{tau} and \code{sigma}.}
#' \item{E}{The conductance-implied covariance matrix returned by
#'   \code{\link{terradish_distance}(covariance = TRUE)}.}
#' \item{theta}{The validated conductance parameter matrix used for simulation.}
#' \item{tau}{The supplied \code{tau}.}
#' \item{sigma}{The supplied \code{sigma}.}
#' \item{nu}{The supplied \code{nu}.}
#'
#' @examples
#' r1 <- terra::rast(nrows = 3, ncols = 3, vals = 1:9)
#' r2 <- terra::rast(nrows = 3, ncols = 3, vals = c(9, 1, 4, 3, 8, 2, 5, 7, 6))
#' covariates <- c(r1, r2)
#' names(covariates) <- c("x1", "x2")
#' pts <- terra::vect(matrix(c(0.5, 0.5,
#'                            1.5, 1.5,
#'                            2.5, 2.5), ncol = 2, byrow = TRUE),
#'                    type = "points")
#' surface <- conductance_surface(covariates, pts, directions = 4)
#' sim <- simulate_covariance_response(
#'   theta = c(x1 = 0.3, x2 = -0.2),
#'   formula = ~ x1 + x2,
#'   data = surface,
#'   tau = 0.8,
#'   sigma = 0.5,
#'   nu = 20,
#'   seed = 1
#' )
#' sim$covariance
#'
#' @export
simulate_covariance_response <- function(theta,
                                         formula,
                                         data,
                                         conductance_model = loglinear_conductance,
                                         tau = 1,
                                         sigma = 0,
                                         nu,
                                         nsim = 1,
                                         seed = NULL,
                                         cores = 1L)
{
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(is.numeric(tau), length(tau) == 1L, is.finite(tau), tau >= 0)
  stopifnot(is.numeric(sigma), length(sigma) == 1L, is.finite(sigma), sigma >= 0)
  stopifnot(is.numeric(nu), length(nu) == 1L, is.finite(nu), nu > 0)
  stopifnot(is.numeric(nsim), length(nsim) == 1L, nsim >= 1)
  stopifnot(length(cores) == 1L, is.numeric(cores), cores >= 1)

  terms_obj <- terms(formula)
  if (!length(attr(terms_obj, "term.labels")))
    stop("`formula` must include at least one conductance covariate.", call. = FALSE)
  model_formula <- reformulate(attr(terms_obj, "term.labels"))

  conductance_model_obj <- conductance_model(model_formula, data$x)
  default <- attr(conductance_model_obj, "default")

  if (is.null(dim(theta)))
    theta <- matrix(theta, nrow = 1L)
  stopifnot(is.matrix(theta))
  stopifnot(nrow(theta) == 1L)
  stopifnot(ncol(theta) == length(default))
  theta <- .validate_theta_grid(theta, names(default))

  if (!is.null(seed))
    set.seed(seed)

  E <- terradish_distance(
    theta = theta,
    formula = model_formula,
    data = data,
    conductance_model = conductance_model,
    conductance = TRUE,
    covariance = TRUE,
    cores = as.integer(cores)
  )$covariance[, , 1]

  Sigma <- tau * E + sigma * diag(nrow(E))
  covariance <- array(NA_real_, dim = c(nrow(E), ncol(E), as.integer(nsim)))
  for (i in seq_len(as.integer(nsim)))
    covariance[, , i] <- rWishart(1, df = as.integer(nu), Sigma = Sigma)[, , 1] / as.integer(nu)

  list(
    covariance = if (as.integer(nsim) == 1L) covariance[, , 1] else covariance,
    Sigma = Sigma,
    E = E,
    theta = theta,
    tau = tau,
    sigma = sigma,
    nu = as.integer(nu)
  )
}

# Internal simulation helper for benchmarking generalized Wishart fits.
wishart_simulate_distance <- function(seed, S, nu)
{
  set.seed(seed)
  dist_from_cov(solve(rWishart(1, nu, S)[,,1]))
}

wishart_simulate_experiment <- function(seed, N, P, K, nu, neval, timingOnly=FALSE)
{
  covariates <- list()
  for(k in 1:K)
  {
    sim <- matrix(.randomfields_rfsimulate(.randomfields_rmexp(scale = 10),
                                           expand.grid(x = 1:N, y = 1:N),
                                           seed = seed)@data[[1]], N, N)
    covariates[[paste0("var", k)]] <- rast(nrows = N, ncols = N,
                                           xmin = 0, xmax = N, ymin = 0, ymax = N)
    values(covariates[[paste0("var", k)]]) <- as.vector(sim)
    covariates[[paste0("var", k)]] <- covariates[[paste0("var", k)]] /
      max(values(covariates[[paste0("var", k)]])[,1], na.rm = TRUE)
  }
  covariates <- do.call(c, covariates)
  names(covariates) <- paste0("var", seq_len(K))

  set.seed(seed)
  demes <- sample(1:(N^2), P)
  coords <- xyFromCell(covariates[[1]], demes)

  surf <- conductance_surface(covariates, coords, directions = 4, saveStack = FALSE)

  beta <- matrix(rnorm(K), 1, K)
  E <- terradish_distance(loglinear_conductance, surf, beta, covariance = TRUE)$covariance[,,1]
  S <- wishart_simulate_distance(seed=seed, nu=nu, S=E)
  S <- S/(max(S)*10)

  # functions
  nograd <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = FALSE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  wgrad <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  whess <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = TRUE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }

  timings_whess <- c()
  timings_wgrad <- c()
  timings_nograd <- c()
  for(i in 1:neval)
  {
    timings_whess <- rbind(timings_whess, system.time(ll <- whess(rep(0,K))))
    timings_wgrad <- rbind(timings_wgrad, system.time(ll <- wgrad(rep(0,K))))
    timings_nograd <- rbind(timings_nograd, system.time(ll <- nograd(rep(0,K))))
  }

  #optimization
  if(!timingOnly){
  opt_newton <- list(fcall = NA, loglik = NA, fit = list(boundary = NA))
  fcall <- 0L
  bqfn <- function(par)
  {
    fcall <<- fcall + 1L
    nograd(par)$objective
  }
  opt_bobyqa <- .nloptr_bobyqa(rep(0, K), bqfn)
  opt_bobyqa$fcall <- fcall
} else {
  opt_bobyqa = list(fcall=NA, value=NA)
  opt_newton = list(fcall=NA, loglik=NA, fit=list(boundary=NA))
  }


  list(timings=list(nograd=timings_nograd, wgrad=timings_wgrad, whess=timings_whess), opt_newton=opt_newton, 
       #opt_lbfgs=opt_lbfgs, 
       opt_bobyqa=opt_bobyqa, seed=seed, beta=beta, K=K, N=N, P=P)
}

run_benchmarks <- function(K = c(1,2,4,8,16,32), N = c(100), P = c(30), reps=5, timingOnly=timingOnly)
{
  set.seed(1)
  seeds <- sample.int(100000, length(K)*length(N)*length(P)*reps)
  z <- 0
  out <- list()
  for(p in P)
  {
    out[[as.character(p)]] <- list()
    for(n in N)
    {
      out[[as.character(p)]][[as.character(n)]] <- list()
      for(k in K)
      {
        out[[as.character(p)]][[as.character(n)]][[as.character(k)]] <- list()
        for(rep in 1:reps)
        {
          z <- z + 1
          try({
            out[[as.character(p)]][[as.character(n)]][[as.character(k)]][[as.character(rep)]] <- 
              wishart_simulate_experiment(seeds[z], n, p, k, 100, 1, timingOnly=timingOnly)
          })
        }
      }
    }
  }
  out
}
#aahh <- run_benchmarks(c(1,2,4,8,16,32), c(200), c(30), 10, timingOnly=TRUE)

extract_benchmark_timing <- function(benchmarks)
{
  out <- c()
  for(p in names(benchmarks))
    for(n in names(benchmarks[[p]]))
      for(k in names(benchmarks[[p]][[n]]))
        for(r in names(benchmarks[[p]][[n]][[k]]))
        {
          tmp <- benchmarks[[p]][[n]][[k]][[r]]
          oo <- data.frame(P=as.numeric(p), N=as.numeric(n), K=as.numeric(k), rep=as.numeric(r),
                           timing_nograd = mean(tmp$timings$nograd[,1]),
                           timing_wgrad = mean(tmp$timings$wgrad[,1]),
                           timing_whess = mean(tmp$timings$whess[,1]),
                           newton_eval  = tmp$opt_newton$fcall,
                           newton_ll    = tmp$opt_newton$loglik,
                           newton_boundary = tmp$opt_newton$fit$boundary,
                           #lbfgs_eval  = tmp$opt_lbfgs$fcall,
                           bobyqa_eval  = tmp$opt_bobyqa$fcall,
                           bobyqa_ll    = tmp$opt_bobyqa$value
                           )
          out <- rbind(out, oo)
        }
  out
}
