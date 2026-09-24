root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkg <- file.path(root, "dev/check/phase7-package")
pkgload::load_all(pkg, quiet = TRUE, export_all = FALSE)
results <- list()
for (file in list.files(file.path(pkg, "inst/examples"), full.names = TRUE)) {
  source(file, local = globalenv())
}
for (name in c("run_spline_recovery_example", "fit_melip_spline_example",
               "run_covariance_response_power_example")) {
  cat("\nEXAMPLE", name, "\n")
  started <- Sys.time()
  value <- tryCatch(get(name)(), error = function(e) list(error = conditionMessage(e)))
  result <- if (!is.null(value$error)) value else if (!is.null(value$fits)) {
    list(model_comparison = value$model_comparison,
      convergence = vapply(value$fits, function(fit) fit$convergence$code, numeric(1)))
  } else {
    list(scenario_summary = value$scenario_summary,
      parameter_summary = value$parameter_summary, planning = value$planning)
  }
  results[[name]] <- list(started = started, finished = Sys.time(), result = result)
  saveRDS(results, file.path(root, "dev/release_0.1.0/phase7_examples.rds"))
  print(results[[name]])
}
