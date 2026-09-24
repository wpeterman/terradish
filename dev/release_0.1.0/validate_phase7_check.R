root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
started <- Sys.time()
result <- tryCatch(devtools::check("dev/check/phase7-package",
  args = c("--no-manual", "--no-tests"), document = FALSE,
  error_on = "never", check_dir = "dev/check/phase7"),
  error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), "dev/release_0.1.0/phase7_check.rds")
print(result)
