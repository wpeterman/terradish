root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
Sys.setenv(NOT_CRAN = "true")
# The plan's Phases 2-5 target comprises R files. Trace all R expressions while
# retaining the normal native build; do not collect gcov data in this run.
options(covr.flags = character(), covr.gcov = "")
started <- Sys.time()
result <- tryCatch(covr::package_coverage("dev/check/phase8-r-coverage", quiet = TRUE),
  error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  scope = "R source coverage; normal native compilation; full test suite",
  session = sessionInfo()), "dev/release_0.1.0/phase8_r_coverage.rds")
if (inherits(result, "coverage")) source("dev/release_0.1.0/report_phase8_r_coverage.R") else print(result)
