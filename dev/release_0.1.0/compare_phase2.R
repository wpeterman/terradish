args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
base <- readRDS("dev/release_0.1.0/core_baseline_0.0.47.rds")$models
new <- readRDS(args[1])$models
for (name in names(base)) {
  a <- base[[name]]
  b <- new[[name]]
  stopifnot(!is.null(b), is.null(b$error))
  cat(name, "logLik:", a$logLik, "->", b$logLik, "\n")
  if (name %in% c("d_wishart_covariance", "g_gaussian_even")) next
  stopifnot(abs(a$logLik - b$logLik) / max(abs(a$logLik), .Machine$double.eps) <= 1e-6,
            max(abs(a$theta - b$theta)) <= 1e-5,
            max(abs(a$se - b$se)) <= 1e-5)
  if (!name %in% c("e_spline", "i_wishart_covariates"))
    stopifnot(max(abs(a$phi - b$phi)) <= 1e-5)
}
# F4 centers spline columns, so a nuisance scale absorbs the constant shift.
# F8 removes kernel normalization; compare the corresponding lambda units.
pkgload::load_all()
data(melip, package = "terradish")
rasters <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
names(rasters) <- c("altitude", "forestcover")
z <- terra::extract(scale_covariates(rasters), terra::unwrap(melip.coords), ID = FALSE)$altitude
kernel_scale <- mean((z - mean(z))^2)
a <- base$i_wishart_covariates$phi
b <- new$i_wishart_covariates$phi
expected <- as.numeric(a)
expected[2] <- expected[2] / kernel_scale
stopifnot(max(abs(as.numeric(b) - expected)) <= 1e-5)
cat("F8 kernel variance:", kernel_scale, "; old/new lambda:", a[2], b[2], "\n")
cat("All Phase 2 baseline contracts passed.\n")
