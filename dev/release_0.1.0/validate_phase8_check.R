root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
started <- Sys.time()
devtools::document("dev/check/phase8-package")
result <- tryCatch(devtools::check("dev/check/phase8-package",
  args = "--as-cran", document = FALSE, error_on = "never",
  check_dir = "dev/check/phase8"), error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), "dev/release_0.1.0/phase8_check.rds")
print(result)
