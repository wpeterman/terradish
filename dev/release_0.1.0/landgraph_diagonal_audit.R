# L3: repeat the supplied exp11_within design using the contrast prototype.
pkgload::load_all()
library(terra)
scratch <- new.env(parent = asNamespace("terradish"))
audit_root <- "C:/Users/peterman.73/OneDrive - The Ohio State University/Research/R_packages/terradish/review_work_20260924/package_audit/core"
sys.source(file.path(audit_root, "helper.R"), scratch)
sys.source("dev/check/phase2/R/wishart_covariance.R", scratch)
with(scratch, {
  landscape <- make_land(30, seed = 2, n_sites = 20)
  surface <- landscape$surf
  E <- terradish_distance(matrix(c(0.6, -0.4), 1), ~x1 + x2, surface,
                           covariance = TRUE)$covariance[, , 1]
  n <- nrow(E)
  Sigma <- 0.05 * (E / mean(diag(E))) + 0.01 * diag(n)
  records <- list()
  for (seed in 1:3) {
    set.seed(seed)
    loci <- 500L
    individuals <- 8L
    p0 <- runif(loci, 0.1, 0.9)
    Z <- t(chol(Sigma)) %*% matrix(rnorm(n * loci), n, loci)
    P <- pmin(pmax(matrix(p0, n, loci, byrow = TRUE) +
                    Z * matrix(sqrt(p0 * (1 - p0)), n, loci, byrow = TRUE), 0.001), 0.999)
    G <- do.call(rbind, lapply(seq_len(n), function(k)
      matrix(rbinom(individuals * loci, 2, rep(P[k, ], each = individuals)), individuals, loci)))
    groups <- rep(seq_len(n), each = individuals)
    for (diagonal in c("within", "gower")) {
      C <- suppressWarnings(cov_from_genetic_data(G, groups = groups, diagonal = diagonal))
      C <- matrix(C, n, n)
      fit <- suppressWarnings(terradish(C ~x1 + x2, surface, loglinear_conductance,
                                         wishart_covariance, nu = loci))
      theta <- coef(fit)
      se <- if (length(theta)) sqrt(diag(-solve(fit$mle$hessian))) else rep(NA_real_, 2)
      if (!length(theta)) theta <- rep(NA_real_, 2)
      record <- data.frame(seed, diagonal, term = c("x1", "x2"), truth = c(0.6, -0.4),
                             estimate = as.numeric(theta), se = as.numeric(se),
                             boundary = fit$fit$boundary)
      record$bias_se <- abs(record$estimate - record$truth) / record$se
      record$wrong_sign <- sign(record$estimate) != sign(record$truth)
      print(record)
      records[[length(records) + 1L]] <- record
      saveRDS(do.call(rbind, records), "dev/release_0.1.0/landgraph_diagonal_audit.rds")
    }
  }
})
