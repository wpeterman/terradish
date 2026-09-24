x <- readRDS("dev/release_0.1.0/phase4_predictive_effect.rds")
for (i in 1:8) for (family in c("wishart", "mlpe")) {
  z <- x[[i]]$cv[[family]]$results
  z <- z[nzchar(z$error), c("model", "fold", "error")]
  if (nrow(z)) { cat("replicate", i, family, "\n"); print(z) }
}
