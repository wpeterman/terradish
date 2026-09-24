# Run this file with Source in RStudio, from the terradish project root.
# This is the literal interactive-session check requested by the release plan.
stopifnot(interactive(), identical(Sys.getenv("RSTUDIO"), "1"))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(root, "dev/check/phase8-final-package/DESCRIPTION")))
receipt <- file.path(root, "dev/release_0.1.0/phase8_rstudio_check.rds")
if (file.exists(receipt)) stop("The RStudio receipt already exists; preserve it.")
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
locale <- Sys.getenv()
bad <- names(locale)[grepl("^LC_", names(locale)) & locale == "C.UTF-8"]
if (length(bad)) Sys.unsetenv(bad)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
started <- Sys.time()
result <- devtools::check(file.path(root, "dev/check/phase8-final-package"),
  args = "--as-cran", document = FALSE, error_on = "never",
  check_dir = file.path(root, "dev/check/phase8-rstudio"))
saveRDS(list(started = started, finished = Sys.time(), result = result,
  session = sessionInfo()), receipt)
print(result)
