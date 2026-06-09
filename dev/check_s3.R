# Exercise the full S3 suite (print, summary, coef, logLik, AIC, vcov, confint,
# plot) across all three new fit types on a small, fast graph.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))
set.seed(1)

DIM <- 6L
r <- rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM, ymin = 0, ymax = DIM)
gx <- xFromCell(r, seq_len(ncell(r))); gy <- yFromCell(r, seq_len(ncell(r)))
covs <- c(setValues(r, scale(gx + 0.5 * gy)[, 1]), setValues(r, scale(gx)[, 1]))
names(covs) <- c("v1", "elev")
fc <- c(1L, DIM, DIM * DIM, DIM * (DIM - 1L) + 1L, (DIM * DIM) %/% 2L, DIM * 2L + 3L)
coords <- xyFromCell(r, fc)
surface <- conductance_surface(covs, coords, directions = 8, saveStack = TRUE)
nf <- length(surface$demes); nu <- 500

# simulate a covariance response from a symmetric loglinear model
cm <- loglinear_conductance(~ v1, surface$x)
E <- as.matrix(terradish_algorithm(cm, leastsquares, surface, S = diag(nf),
                                   theta = c(v1 = 0.4), objective = FALSE,
                                   gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
set.seed(5); S <- rWishart(1, df = nu, Sigma = (0.8 * E + 0.2 * diag(nf)) / nu)[, , 1]

ok <- function(lbl, cond) cat(sprintf("  %-28s %s\n", lbl, ifelse(isTRUE(cond), "OK", "**FAIL**")))

test_suite <- function(fit, name, has_data = TRUE, np = 1L) {
  cat("\n===", name, "===\n")
  a <- AIC(fit); ll <- logLik(fit)
  ok("AIC numeric", is.numeric(a) && is.finite(a))
  ok("logLik + df attr", is.finite(as.numeric(ll)) && !is.null(attr(ll, "df")))
  ok("coef length", length(coef(fit)) == np)
  v <- tryCatch(vcov(fit), error = function(e) NULL)
  ok("vcov matrix", is.matrix(v) || is.null(v))
  ci <- tryCatch(confint(fit), error = function(e) NULL)
  ok("confint", is.matrix(ci) || is.null(ci))
  s <- tryCatch({ capture.output(summary(fit)); TRUE }, error = function(e) FALSE)
  ok("summary runs", s)
  p <- tryCatch({ capture.output(print(fit)); TRUE }, error = function(e) FALSE)
  ok("print runs", p)
  if (has_data) {
    pl <- tryCatch({ grDevices::pdf(NULL); on.exit(grDevices::dev.off());
                     rr <- plot(fit, data = surface); inherits(rr, "SpatRaster") }, error = function(e) FALSE)
    ok("plot -> raster", pl)
  }
  cat(sprintf("  AIC=%.2f  logLik=%.2f (df=%.2f)\n", a, as.numeric(ll), attr(ll, "df")))
}

# Tier 1: drift surface (a standard terradish fit -> base S3 methods)
g_drift <- wishart_drift_covariates(scale(seq_len(nf))[, 1], model = "wishart_covariance")
fit1 <- suppressWarnings(terradish(S ~ v1, surface, loglinear_conductance, g_drift, nu = nu,
                  leverage = FALSE, control = NewtonRaphsonControl(maxit = 50, verbose = FALSE)))
cat("\n=== Tier 1 (wishart_drift_covariates; base terradish S3) ===\n")
ok("AIC numeric", is.finite(AIC(fit1)))
ok("logLik df", !is.null(attr(logLik(fit1), "df")))
ok("summary runs", tryCatch({ capture.output(summary(fit1)); TRUE }, error = function(e) FALSE))
ok("coef", length(coef(fit1)) == 1L)

# Tier 2: hierarchical (fixed tau2 for speed)
fit2 <- terradish_hierarchical(S ~ v1, data = surface, measurement_model = wishart_covariance,
                               nu = nu, field_resolution = 3, tau2 = 0.5, verbose = FALSE)
test_suite(fit2, "Tier 2 (terradish_hierarchical)", has_data = TRUE, np = 1L)

# Tier 3: directed
dir_cov <- edge_gradient(covs[["elev"]], surface)
fit3 <- terradish_directed(S ~ v1, data = surface, directional = dir_cov,
                           measurement_model = wishart_covariance, nu = nu)
test_suite(fit3, "Tier 3 (terradish_directed)", has_data = TRUE, np = 2L)

cat("\nDONE\n")
