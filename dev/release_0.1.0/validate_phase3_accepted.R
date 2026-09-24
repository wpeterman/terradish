# Do not modify this script or its candidate while either validation is running.
mode <- commandArgs(trailingOnly = TRUE)
stopifnot(length(mode) == 1L, mode %in% c("test", "check"))
root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
started <- Sys.time()
result <- tryCatch({
  if (mode == "test") {
    as.data.frame(devtools::test("dev/check/phase3-accepted-package", reporter = "summary",
                                 stop_on_failure = FALSE))
  } else {
    devtools::check("dev/check/phase3-accepted-package", args = "--no-manual",
      document = FALSE, error_on = "never", check_dir = "dev/check/phase3-accepted")
  }
}, error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
             session = sessionInfo()),
  file.path("dev/release_0.1.0", paste0("phase3_accepted_", mode, ".rds")))
if (is.data.frame(result)) {
  print(colSums(result[, c("passed", "failed", "error", "warning", "skipped")]))
} else {
  print(result)
}
