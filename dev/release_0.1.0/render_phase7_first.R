# Render these frozen vignette drafts against the phase 7 candidate.
root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkg <- file.path(root, "dev/check/phase7-package")
out <- file.path(root, "dev/check/phase7-render-first")
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
pkgload::load_all(pkg, quiet = TRUE, export_all = FALSE)
names <- c("getting-started", "large-landscapes", "model-comparison",
           "wishart-covariance", "spline-conductance", "gaussian-scale-optimization")
results <- list()
for (name in names) {
  cat("\nRENDER", name, "\n")
  started <- Sys.time()
  warnings <- character()
  result <- tryCatch(withCallingHandlers(
    rmarkdown::render(file.path(pkg, "vignettes", paste0(name, ".Rmd")),
      output_dir = out, envir = new.env(parent = globalenv()), quiet = FALSE),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) list(error = conditionMessage(e)))
  results[[name]] <- list(started = started, finished = Sys.time(),
                          result = result, warnings = warnings)
  saveRDS(results, file.path(root, "dev/release_0.1.0/phase7_render_first.rds"))
  print(results[[name]])
}
