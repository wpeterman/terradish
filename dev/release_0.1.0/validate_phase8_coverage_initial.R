root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
Sys.setenv(NOT_CRAN = "true")
started <- Sys.time()
result <- tryCatch(covr::package_coverage("dev/check/phase8-package", quiet = TRUE),
  error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), "dev/release_0.1.0/phase8_coverage.rds")
if (inherits(result, "coverage")) {
  print(covr::percent_coverage(result))
  rows <- covr::tally_coverage(result, by = "file")
  write.csv(rows, "dev/release_0.1.0/phase8_coverage_by_file.csv", row.names = FALSE)
  print(rows)
} else print(result)
