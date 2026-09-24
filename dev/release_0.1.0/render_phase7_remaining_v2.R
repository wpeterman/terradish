root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkg <- file.path(root, "dev/check/phase7-package")
out <- file.path(root, "dev/check/phase7-render-remaining-v2")
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
pkgload::load_all(pkg, quiet = TRUE, export_all = FALSE)
results <- list()
for (name in c("ibe-ibr-workflow", "simulation-design")) {
  cat("\nRENDER", name, "\n")
  started <- Sys.time()
  env <- new.env(parent = globalenv())
  result <- tryCatch(rmarkdown::render(
    file.path(pkg, "vignettes", paste0(name, ".Rmd")), output_dir = out,
    envir = env, quiet = FALSE), error = function(e) list(error = conditionMessage(e)))
  convergence <- if (exists("fits", env, inherits = FALSE))
    vapply(env$fits, function(fit) fit$convergence$code, numeric(1)) else NULL
  results[[name]] <- list(started = started, finished = Sys.time(),
    result = result, convergence = convergence)
  saveRDS(results, file.path(root, "dev/release_0.1.0/phase7_render_remaining_v2.rds"))
  print(results[[name]])
}
