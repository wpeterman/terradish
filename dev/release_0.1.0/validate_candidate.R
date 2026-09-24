# The candidate is frozen before this runner starts. Never edit it mid-run.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c("test", "check"))
root <- getwd()
stage <- file.path(root, "dev/check/phase2-package")
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
started <- Sys.time()
result <- tryCatch({
  if (args == "test")
    as.data.frame(devtools::test(stage, reporter = "summary", stop_on_failure = FALSE))
  else
    devtools::check(stage, args = "--no-manual", document = FALSE,
                     error_on = "never", check_dir = file.path(root, "dev/check/phase2"))
}, error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
             session = sessionInfo()),
          file.path(root, "dev/release_0.1.0", paste0("phase2_", args, ".rds")))
if (is.data.frame(result)) {
  print(colSums(result[, c("passed", "failed", "error", "warning", "skipped")]))
} else {
  print(result)
}
