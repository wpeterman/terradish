source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/cv/setup.R")
suppressWarnings({
coords_m <- crds(pts)
rf <- terradish_folds(coords_m, k = 4, method = "random", seed = 1)
ck <- tempfile("phase8-cv-checkpoint-", fileext = ".rds")
S <- simulate_covariance_response(theta = c(x1 = 0.8), formula = ~x1, data = surface, tau = 1, sigma = 0.1, nu = 25, seed = 5)$covariance
a <- terradish_cv_folds(surface, list(m = S ~ x1), folds = rf, model = wishart_covariance, nu = 25, checkpoint = ck)
S <- simulate_covariance_response(theta = c(x1 = -0.8), formula = ~x1, data = surface, tau = 1, sigma = 0.1, nu = 25, seed = 99)$covariance
message <- tryCatch({
  terradish_cv_folds(surface, list(m = S ~ x1), folds = rf,
    model = wishart_covariance, nu = 25, checkpoint = ck)
  "ERROR: changed response was accepted"
}, error = conditionMessage)
cat("Changed-response resume:", message, "\n")
stopifnot(grepl("not compatible", message, fixed = TRUE))
unlink(ck)
})
