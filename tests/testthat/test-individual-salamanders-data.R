test_that("individual marker panels preserve matched sampling and spatial support", {
  env <- new.env(parent = emptyenv())
  utils::data("individual_salamanders", package = "terradish", envir = env)
  x <- env$individual_salamanders
  terms <- terra::rast(system.file("extdata", "individual-salamanders",
    "drying_terms_10m.tif", package = "terradish"))
  expect_equal(terra::res(terms), c(10, 10))
  expect_true(terra::same.crs(terms, "EPSG:26917"))
  expect_identical(names(terms), c("drying1", "drying2"))
  for (z in x$data) {
    expect_equal(dim(z$snp_dosage[["2500"]]), c(160L, 2500L))
    expect_equal(dim(z$microsat_alleles), c(160L, 48L))
    for (G in z$snp_dosage) {
      expect_identical(rownames(G), z$samples$sample_id)
      expect_true(all(G %in% 0:2))
      expect_true(all(pmin(colMeans(G) / 2, 1 - colMeans(G) / 2) >= 0.05))
      expect_identical(G, z$snp_dosage[["2500"]][, colnames(G), drop = FALSE])
    }
    expect_identical(rownames(z$microsat_alleles), z$samples$sample_id)
    expect_true(all(z$microsat_alleles %in% 5:45))
    expect_equal(as.integer(table(z$samples$fold)), rep(40L, 4))
    expect_equal(as.integer(table(z$samples$stratum)), rep(10L, 16))
    xy <- as.matrix(z$samples[c("easting", "northing")])
    expect_gte(min(stats::dist(xy)), 10)
    expect_equal(length(unique(terra::cellFromXY(terms, xy))), 160L)
  }
})

test_that("public covariance helpers reconstruct the documented genetic quantity", {
  env <- new.env(parent = emptyenv())
  utils::data("individual_salamanders", package = "terradish", envir = env)
  for (z in env$individual_salamanders$data) {
    G <- z$snp_dosage[["2500"]]
    C <- cov_from_genetic_data(G, scale = FALSE, normalize = "none") / (2 * ncol(G))
    # Independent direct construction checks the public helper's normalization.
    centered <- sweep(G, 2, colMeans(G), "-")
    expect_equal(as.vector(C), as.vector(tcrossprod(centered) / (2 * ncol(G))),
      tolerance = 1e-12)
    expect_equal(as.vector(C), as.vector(z$covariance$SNP2500), tolerance = 1e-12)
    A <- z$microsat_alleles
    L <- ncol(A) / 2
    C_ms <- cov_from_genetic_data(A, input = "allele_calls",
      loci = rep(seq_len(L), each = 2), scale = FALSE,
      normalize = "none") / (4 * L)
    expect_equal(as.vector(C_ms), as.vector(z$covariance$microsat24), tolerance = 1e-12)
    for (S in list(C, C_ms)) {
      D <- as.matrix(dist_from_cov(S))
      expect_gte(min(D), -1e-12)
      expect_lte(max(D), 2 + 1e-12)
      expect_lt(max(abs(rowSums(S))), 1e-10)
      expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
    }
    # The weighted biallelic API is coherent, but represents another response.
    N <- rep(2, nrow(G))
    C_weighted <- cov_from_biallelic(G, N)
    expect_equal(dist_from_biallelic(G, N), dist_from_cov(C_weighted))
    expect_gt(max(as.matrix(dist_from_cov(C_weighted))), 2)
  }
})

test_that("stored prediction comparisons retain the mixed microsatellite evidence", {
  files <- system.file("extdata", "individual-salamanders", package = "terradish")
  results <- readRDS(file.path(files, "assessment.rds"))
  predictions <- readRDS(file.path(files, "predictions.rds"))
  expect_length(predictions, 40L)
  for (comparison in results$original) {
    errors <- lapply(c("uniform", "IBR"), function(model) {
      records <- lapply(1:4, function(fold) predictions[[paste(
        comparison$arm, 5000, comparison$panel, model, fold, sep = "__")]])
      expect_true(all(vapply(records, nrow, integer(1)) == 780L))
      pooled <- do.call(rbind, records)
      expect_true(all(is.finite(as.matrix(pooled[c("observed", "predicted")]))))
      expect_true(all(pooled$predicted >= 0 & pooled$predicted <= 2))
      sqrt(mean((pooled$observed - pooled$predicted)^2))
    })
    expect_equal(errors[[1]], comparison$uniform_RMSE, tolerance = 1e-12)
    expect_equal(errors[[2]], comparison$IBR_RMSE, tolerance = 1e-12)
    expect_equal(1 - errors[[2]] / errors[[1]], comparison$RMSE_reduction, tolerance = 1e-12)
  }
  for (comparisons in list(results$original, results$independent_assay)) {
    expect_gt(comparisons$compensatory_500__SNP2500$RMSE_reduction, 0.05)
    expect_gt(comparisons$compensatory_500__SNP2500$total_score_gain, 0)
    expect_gte(comparisons$compensatory_500__SNP2500$positive_folds, 3)
    expect_lt(comparisons$compensatory_500__microsat24$RMSE_reduction, 0)
    expect_lt(comparisons$control_500__SNP2500$RMSE_reduction, 0)
  }
  expect_equal(dim(results$block_contrasts), c(3L, 2L))
  expect_true(all(results$block_contrasts > 1))
  assay <- readRDS(file.path(files, "independent-assay.rds"))
  env <- new.env(parent = emptyenv())
  utils::data("individual_salamanders", package = "terradish", envir = env)
  for (arm in names(assay)) {
    ids <- env$individual_salamanders$data[[arm]]$samples$sample_id
    expect_identical(rownames(assay[[arm]]$snp), ids)
    expect_identical(rownames(assay[[arm]]$microsat), ids)
  }
})
