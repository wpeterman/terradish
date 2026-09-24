# Build a separate, reviewable Phase 2 candidate while C1 checks finish.
root <- getwd()
stage <- file.path(root, "dev/check/phase2-package")
if (dir.exists(stage)) stop("Candidate already exists; inspect it before replacing anything.")
dir.create(stage, recursive = TRUE)
paths <- unique(system2("git", c("--no-optional-locks", "ls-files", "--cached",
                                 "--others", "--exclude-standard"), stdout = TRUE))
paths <- paths[file.exists(paths) & !grepl("^(dev/|\\.git/)", paths)]
for (path in paths) {
  dir.create(dirname(file.path(stage, path)), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(path, file.path(stage, path)))
}
for (path in list.files("dev/check/phase2/R", full.names = TRUE))
  stopifnot(file.copy(path, file.path(stage, "R", basename(path)), overwrite = TRUE))
for (path in list.files("dev/check/phase2/tests", pattern = "\\.R$", full.names = TRUE))
  stopifnot(file.copy(path, file.path(stage, "tests/testthat", basename(path)), overwrite = TRUE))
path <- file.path(stage, "DESCRIPTION")
text <- readLines(path)
text <- sub("landgraph (>= 0.0.1)", "landgraph (>= 0.0.3)", text, fixed = TRUE)
text <- sub("Version: 0.0.48", "Version: 0.0.50", text, fixed = TRUE)
writeLines(text, path)
# Retain explicit coverage of the old within-diagonal option while adopting
# the owner's new Gower default and covariance metadata.
stopifnot(file.copy("../landgraph/tests/testthat/test-genetic-covariance.R",
                     file.path(stage, "tests/testthat/test-genetic-covariance.R"), overwrite = TRUE))
path <- file.path(stage, "tests/testthat/test-wishart-covariates.R")
text <- readLines(path)
text <- gsub("kernel_altitude", "kernel_absdiff_altitude", text, fixed = TRUE)
text <- gsub("kernel_moisture", "kernel_absdiff_moisture", text, fixed = TRUE)
text <- gsub("lambda_altitude", "lambda_absdiff_altitude", text, fixed = TRUE)
writeLines(text, path)
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
devtools::document(stage)
cat("Prepared Phase 2 candidate:", stage, "\n")
