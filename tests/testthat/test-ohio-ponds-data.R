test_that('pond genotype panels and pedigree metadata retain their alignment', {
  env <- new.env(parent = emptyenv())
  utils::data('ohio_ponds', package = 'terradish', envir = env)
  x <- env$ohio_ponds
  expect_equal(dim(x$snp), c(480L, 2500L))
  expect_equal(dim(x$microsat), c(480L, 48L))
  expect_equal(dim(x$raw_snp), c(780L, 2500L))
  expect_identical(rownames(x$snp), x$samples$individual)
  expect_identical(rownames(x$microsat), x$samples$individual)
  expect_identical(rownames(x$raw_snp), x$raw_samples$individual)
  expect_identical(rownames(x$raw_microsat), x$raw_samples$individual)
  expect_identical(colnames(x$snp), colnames(x$raw_snp))
  expect_identical(colnames(x$microsat), colnames(x$raw_microsat))
  retained <- match(x$samples$individual, x$raw_samples$individual)
  expect_false(anyNA(retained))
  expect_identical(x$snp, x$raw_snp[retained, , drop = FALSE])
  expect_identical(x$microsat, x$raw_microsat[retained, , drop = FALSE])
  expect_true(all(table(factor(x$samples$pond, levels = x$ponds$pond)) == 20))
  expect_true(all(table(x$raw_samples$pond) >= 30))
  expect_equal(length(unique(x$microsat_loci)), 24L)
  expect_true(all(table(x$microsat_loci) == 2L))
  expect_true(all(x$snp %in% 0:2))
  expect_true(all(diff(x$snp_positions) >= 10000))
  frequencies <- colMeans(x$snp) / 2
  expect_true(all(pmin(frequencies, 1 - frequencies) >= 0.05))
})

test_that('packed pond landscapes are portable and retain projected units', {
  env <- new.env(parent = emptyenv())
  utils::data('ohio_ponds', package = 'terradish', envir = env)
  x <- env$ohio_ponds
  raw <- terra::unwrap(x$covariates)
  context <- terra::unwrap(x$forest500)
  expect_equal(terra::res(raw), c(100, 100))
  expect_equal(c(nrow(raw), ncol(raw)), c(102, 102))
  expect_identical(names(raw), c('forest', 'slope', 'elevation'))
  expect_identical(names(context), 'forest500')
  expect_true(terra::same.crs(raw, 'EPSG:26917'))
  expect_true(terra::compareGeom(raw[[1]], context, stopOnError = FALSE))
  expect_false(anyNA(terra::values(raw)))
  expect_true(all(terra::values(raw[['forest']]) >= 0 &
                  terra::values(raw[['forest']]) <= 1))
  expect_identical(as.integer(terra::cellFromXY(raw, x$ponds[c('x', 'y')])),
                   as.integer(x$ponds$cell))
  expect_true(all(x$ponds$forest_buffer_100m >= 0.9))
  expect_true(all(x$ponds$forest_buffer_300m >= 0.6))
  expect_true(all(x$ponds$nlcd_30m %in% c(41, 42, 43)))
  path <- system.file('extdata', 'ohio-ponds', 'landcover100.rds', package = 'terradish')
  expect_true(nzchar(path))
  cover <- terra::unwrap(readRDS(path))
  expect_true(is.factor(cover))
  expect_true(terra::compareGeom(raw[[1]], cover, stopOnError = FALSE))
  expect_equal(levels(cover)[[1]]$class, c('forest', 'nonforest'))
})

test_that('the documented SNP frequency response has valid contrast geometry', {
  env <- new.env(parent = emptyenv())
  utils::data('ohio_ponds', package = 'terradish', envir = env)
  x <- env$ohio_ponds
  groups <- factor(x$samples$pond, levels = x$ponds$pond)
  Y <- rowsum(x$snp, groups, reorder = FALSE)
  N <- 2 * as.numeric(table(groups))
  S <- cov_from_biallelic(Y, N) / N[1]
  p <- colSums(Y) / sum(N)
  Z <- sweep(sweep(Y / N, 2, p), 2, sqrt(p * (1 - p)), '/')
  expected <- tcrossprod(Z) / ncol(Z)
  expect_equal(unname(S), unname(expected), ignore_attr = TRUE, tolerance = 1e-12)
  expect_lt(max(abs(rowSums(S))), 1e-12)
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
  D <- dist_from_cov(S)
  expect_silent(check_distance_response(D))
  H <- diag(nrow(S)) - 1 / nrow(S)
  expect_equal(-0.5 * H %*% D %*% H, S, ignore_attr = TRUE, tolerance = 1e-12)
  S_ms <- cov_from_genetic_data(x$microsat, groups, input = 'allele_calls',
    loci = x$microsat_loci, diagonal = 'gower', normalize = 'features')
  expect_equal(dim(S_ms), c(24, 24))
  expect_silent(check_distance_response(dist_from_cov(S_ms)))
})

test_that('provided pond blocks have no known immediate pedigree leakage', {
  env <- new.env(parent = emptyenv())
  utils::data('ohio_ponds', package = 'terradish', envir = env)
  x <- env$ohio_ponds
  samples <- x$samples
  fold <- x$ponds$fold[match(samples$pond, x$ponds$pond)]
  family <- lapply(seq_len(nrow(samples)), function(i)
    setdiff(unlist(samples[i, c('pedigree_id', 'parent_1', 'parent_2')]), -1))
  # Check first-degree relationships, not a claim of deeper unrelatedness.
  leakage <- vapply(seq_len(nrow(samples) - 1L), function(i) {
    j <- seq.int(i + 1L, nrow(samples))
    any(vapply(j[fold[j] != fold[i]], function(k)
      length(intersect(family[[i]], family[[k]])) > 0, logical(1)))
  }, logical(1))
  expect_false(any(leakage))
  expect_equal(sort(unique(fold)), 1:3)
})
