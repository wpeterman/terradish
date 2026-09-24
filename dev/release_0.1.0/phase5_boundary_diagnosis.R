.libPaths(c("dev/check/phase2-library", .libPaths()))
pkgload::load_all("dev/check/phase5-package", quiet = TRUE)
source("dev/check/phase5-package/tests/testthat/helper-robustness.R")
g <- robustness_surface()
S <- diag(length(g$demes))
diagnostic <- function(...) withCallingHandlers(wishart_covariance(...),
 error = function(e) {
   frames <- sys.frames()
   for (i in seq_along(frames)) {
     f <- frames[[i]]
     if (exists("dgrad_dtau", f, inherits = FALSE)) {
       print(sys.calls()[[i]])
       for (nm in c("E", "A", "B", "ASA", "dgrad_dtau", "phi", "nu")) {
         if (exists(nm, f, inherits = FALSE)) {
           x <- get(nm, f)
           print(list(name = nm, dim = dim(x), value = if (length(x) < 10) x else range(x)))
         }
       }
     }
   }
 })

attributes(diagnostic) <- attributes(wishart_covariance)
terradish(S ~ 1, g, measurement_model = diagnostic)
