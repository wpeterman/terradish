# Resume validation on the current checkout without overwriting earlier evidence.
# Run from the package root after clearing inherited C.UTF-8 locale variables.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, grepl("^[A-Za-z0-9_-]+$", args))
label <- args[[1]]
root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
check_dir <- file.path("dev/check", label)
receipt <- file.path("dev/release_0.1.0", paste0(label, ".rds"))
if (dir.exists(check_dir) || file.exists(receipt))
  stop("Choose a new validation label; existing evidence will not be overwritten.")
started <- Sys.time()
head <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
source_files <- c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
  list.files("man", full.names = TRUE),
  list.files("tests/testthat", pattern = "\\.R$", full.names = TRUE))
hashes <- tools::md5sum(source_files)
result <- tryCatch(devtools::check(".", args = "--no-manual", document = FALSE,
  error_on = "never", check_dir = check_dir), error = function(e)
    list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), head = head,
  source_hashes = hashes, result = result, session = sessionInfo()), receipt)
print(result)
