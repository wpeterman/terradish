expected_gower_covariance <- function(centroids)
{
  centroid_dist2 <- as.matrix(stats::dist(centroids))^2
  row_mean <- rowMeans(centroid_dist2)
  col_mean <- colMeans(centroid_dist2)

  -0.5 * (
    centroid_dist2 -
      outer(row_mean, rep(1, length(row_mean))) -
      outer(rep(1, length(col_mean)), col_mean) +
      mean(centroid_dist2)
  )
}

test_that("cov_from_genetic_data constructs Dyer-style covariance from features", {
  x <- matrix(c(0, -1,
                0,  1,
                1, -1,
                1,  1,
                2, -1,
                2,  1),
              ncol = 2, byrow = TRUE)
  groups <- rep(c("pop_b", "pop_a", "pop_c"), each = 2)

  cov <- cov_from_genetic_data(x, groups, center = FALSE, scale = FALSE)

  centroids <- matrix(c(0, 0,
                        1, 0,
                        2, 0),
                      ncol = 2, byrow = TRUE,
                      dimnames = list(c("pop_b", "pop_a", "pop_c"), NULL))
  gower_cov <- expected_gower_covariance(centroids)
  within_variance <- c(pop_b = 2, pop_a = 2, pop_c = 2)
  expected <- gower_cov
  diag(expected) <- within_variance

  expect_equal(cov[, ], expected)
  expect_equal(rownames(cov), names(within_variance))
  expect_equal(attr(cov, "within_variance"), within_variance)
  expect_equal(attr(cov, "centroid_distance2"), as.matrix(stats::dist(centroids))^2)
})

test_that("cov_from_genetic_data constructs individual-level covariance", {
  x <- matrix(c(0, -1,
                0,  1,
                1, -1,
                1,  1),
              ncol = 2, byrow = TRUE,
              dimnames = list(c("ind1", "ind2", "ind3", "ind4"), NULL))

  cov <- cov_from_genetic_data(x, center = FALSE, scale = FALSE)

  expect_equal(cov[, ], expected_gower_covariance(x))
  expect_equal(rownames(cov), rownames(x))
  expect_equal(attr(cov, "level"), "individual")
  expect_equal(attr(cov, "diagonal"), "gower")
  expect_null(attr(cov, "within_variance"))
  expect_equal(attr(cov, "unit_size"), setNames(rep(1L, nrow(x)), rownames(x)))

  cov_ids <- cov_from_genetic_data(
    x,
    groups = rownames(x),
    center = FALSE,
    scale = FALSE
  )
  expect_equal(cov_ids[, ], cov[, ])
  expect_equal(attr(cov_ids, "level"), "individual")
  expect_equal(attr(cov_ids, "diagonal"), "gower")
  expect_error(
    cov_from_genetic_data(
      x,
      groups = rownames(x),
      center = FALSE,
      scale = FALSE,
      diagonal = "within"
    ),
    "at least one population"
  )
})

test_that("cov_from_genetic_data can keep the Gower diagonal", {
  x <- matrix(c(0, -1,
                0,  1,
                1, -1,
                1,  1,
                2, -1,
                2,  1),
              ncol = 2, byrow = TRUE)
  groups <- rep(c("pop1", "pop2", "pop3"), each = 2)

  cov <- cov_from_genetic_data(
    x,
    groups,
    center = FALSE,
    scale = FALSE,
    diagonal = "gower"
  )
  centroids <- matrix(c(0, 0,
                        1, 0,
                        2, 0),
                      ncol = 2, byrow = TRUE,
                      dimnames = list(c("pop1", "pop2", "pop3"), NULL))

  expect_equal(cov[, ], expected_gower_covariance(centroids))
  expect_equal(attr(cov, "diagonal"), "gower")
})

test_that("cov_from_genetic_data can normalize by retained features", {
  x <- matrix(c(0, -1,
                0,  1,
                1, -1,
                1,  1,
                2, -1,
                2,  1),
              ncol = 2, byrow = TRUE)
  groups <- rep(c("pop1", "pop2", "pop3"), each = 2)

  cov_raw <- cov_from_genetic_data(x, groups, center = FALSE, scale = FALSE)
  cov_norm <- cov_from_genetic_data(
    x,
    groups,
    center = FALSE,
    scale = FALSE,
    normalize = "features"
  )

  expect_equal(cov_norm[, ], cov_raw[, ] / ncol(x))
  expect_equal(
    attr(cov_norm, "within_variance"),
    attr(cov_raw, "within_variance") / ncol(x)
  )
  expect_equal(
    attr(cov_norm, "unit_distance2"),
    attr(cov_raw, "unit_distance2") / ncol(x)
  )
  expect_equal(attr(cov_norm, "normalizer"), ncol(x))
})

test_that("cov_from_genetic_data converts allele calls to locus-allele dosage", {
  calls <- data.frame(
    loc1_a = c(100, 100, 102, 102, 104, 104),
    loc1_b = c(100, 102, 102, 104, 104, 100),
    loc2_a = c(200, 202, 200, 202, 204, 204),
    loc2_b = c(202, 202, 204, 204, 204, 200)
  )
  groups <- rep(c("pop1", "pop2", "pop3"), each = 2)

  dosage <- cbind(
    "loc1:100" = c(2, 1, 0, 0, 0, 1),
    "loc1:102" = c(0, 1, 2, 1, 0, 0),
    "loc1:104" = c(0, 0, 0, 1, 2, 1),
    "loc2:200" = c(1, 0, 1, 0, 0, 1),
    "loc2:202" = c(1, 2, 0, 1, 0, 0),
    "loc2:204" = c(0, 0, 1, 1, 2, 1)
  )

  cov_calls <- cov_from_genetic_data(
    calls,
    groups,
    input = "allele_calls",
    loci = c("loc1", "loc1", "loc2", "loc2"),
    center = FALSE,
    scale = FALSE
  )
  cov_dosage <- cov_from_genetic_data(
    dosage,
    groups,
    center = FALSE,
    scale = FALSE
  )

  expect_equal(cov_calls[, ], cov_dosage[, ])
  expect_equal(attr(cov_calls, "retained_features"), colnames(dosage))
})

test_that("cov_from_genetic_data converts individual allele calls", {
  calls <- data.frame(
    loc1_a = c(100, 100, 102, 102),
    loc1_b = c(100, 102, 102, 104),
    loc2_a = c(200, 202, 200, 202),
    loc2_b = c(202, 202, 204, 204)
  )

  dosage <- cbind(
    "loc1:100" = c(2, 1, 0, 0),
    "loc1:102" = c(0, 1, 2, 1),
    "loc1:104" = c(0, 0, 0, 1),
    "loc2:200" = c(1, 0, 1, 0),
    "loc2:202" = c(1, 2, 0, 1),
    "loc2:204" = c(0, 0, 1, 1)
  )

  cov_calls <- cov_from_genetic_data(
    calls,
    input = "allele_calls",
    loci = c("loc1", "loc1", "loc2", "loc2"),
    center = FALSE,
    scale = FALSE
  )
  cov_dosage <- cov_from_genetic_data(
    dosage,
    center = FALSE,
    scale = FALSE
  )

  expect_equal(cov_calls[, ], cov_dosage[, ])
  expect_equal(attr(cov_calls, "level"), "individual")
})

test_that("cov_from_genetic_data imputes missing allele calls to locus modes", {
  calls <- data.frame(
    loc1_a = c(100, NA, 102, 102),
    loc1_b = c(100, 100, 102, 104),
    loc2_a = c(200, 202, 200, 202),
    loc2_b = c(202, NA, 204, 204)
  )
  imputed_calls <- calls
  imputed_calls$loc1_a[2] <- 100
  imputed_calls$loc2_b[2] <- 202

  expect_message(
    cov_missing <- cov_from_genetic_data(
      calls,
      input = "allele_calls",
      loci = c("loc1", "loc1", "loc2", "loc2"),
      center = FALSE,
      scale = FALSE
    ),
    "Imputed 2 missing allele call"
  )
  cov_imputed <- cov_from_genetic_data(
    imputed_calls,
    input = "allele_calls",
    loci = c("loc1", "loc1", "loc2", "loc2"),
    center = FALSE,
    scale = FALSE
  )

  expect_equal(cov_missing[, ], cov_imputed[, ])
  expect_equal(attr(cov_missing, "imputed_allele_calls"),
               c(loc1 = 1L, loc2 = 1L))
  expect_equal(attr(cov_missing, "imputed_modal_alleles"),
               c(loc1 = "100", loc2 = "202"))
})

test_that("cov_from_genetic_data validates inputs and drops constant features", {
  expect_error(
    cov_from_genetic_data(matrix(1:4, nrow = 2), "pop1"),
    "`groups` must have one entry per row"
  )
  expect_error(
    cov_from_genetic_data(data.frame(marker = c("A", "B")), c("pop1", "pop2")),
    "`x` must be numeric"
  )
  expect_error(
    cov_from_genetic_data(
      matrix(c(100, 102), ncol = 1),
      c("pop1", "pop2"),
      input = "allele_calls"
    ),
    "`loci` is required"
  )
  expect_error(
    cov_from_genetic_data(
      matrix(NA, nrow = 2, ncol = 1),
      input = "allele_calls",
      loci = "loc1"
    ),
    "all calls are missing"
  )
  expect_error(
    cov_from_genetic_data(matrix(1:4, nrow = 4), rep("pop1", 4)),
    "at least two populations"
  )
  expect_error(
    cov_from_genetic_data(matrix(1:4, nrow = 2), diagonal = "within"),
    "requires population groups"
  )

  expect_warning(
    cov <- cov_from_genetic_data(
      cbind(variable = c(1, 2, 3, 4), constant = 1),
      rep(c("pop1", "pop2"), each = 2)
    ),
    "Dropping constant"
  )
  expect_equal(attr(cov, "retained_features"), "variable")
})

test_that("cov_from_biallelic supports individual genotypes by default", {
  Y <- matrix(c(0, 1, 2,
                2, 1, 0,
                1, 1, 1),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("ind1", "ind2", "ind3"),
                              paste0("snp", 1:3)))

  cov_default <- cov_from_biallelic(Y)
  cov_explicit <- cov_from_biallelic(Y, matrix(2, nrow(Y), ncol(Y)))

  expect_equal(cov_default, cov_explicit)
  expect_equal(rownames(cov_default), rownames(Y))
  expect_equal(colnames(cov_default), rownames(Y))
})

test_that("cov_from_biallelic drops or errors on monomorphic loci", {
  Y <- matrix(c(0, 2,
                1, 2,
                2, 2),
              nrow = 3, byrow = TRUE,
              dimnames = list(NULL, c("variable", "fixed")))

  expect_warning(
    cov <- cov_from_biallelic(Y),
    "Dropping monomorphic"
  )
  expect_equal(cov, cov_from_biallelic(Y[, "variable", drop = FALSE]))
  expect_error(
    cov_from_biallelic(Y, monomorphic = "error"),
    "pooled allele frequency"
  )
  expect_error(
    cov_from_biallelic(matrix(2, nrow = 3, ncol = 2)),
    "No variable biallelic loci"
  )
})
