root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
# Leave NOT_CRAN unset while building, to exercise an ordinary source build.
Sys.unsetenv("NOT_CRAN")
started <- Sys.time()
result <- tryCatch(devtools::check("dev/check/phase8-package",
  args = c("--as-cran", "--no-tests"), document = FALSE, error_on = "never",
  check_dir = "dev/check/phase8-doc-outputs"),
  error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), "dev/release_0.1.0/phase8_doc_outputs_check.rds")
print(result)
