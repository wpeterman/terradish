root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
Sys.setenv(NOT_CRAN = "true")
started <- Sys.time()
result <- tryCatch(covr::package_coverage("dev/check/phase8-coverage", quiet = TRUE),
  error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), "dev/release_0.1.0/phase8_coverage.rds")
if (inherits(result, "coverage")) {
  print(covr::percent_coverage(result))
  source("dev/release_0.1.0/report_phase8_coverage.R")
} else print(result)
