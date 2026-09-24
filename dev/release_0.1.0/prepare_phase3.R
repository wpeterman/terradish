root <- getwd()
stage <- file.path(root, "dev/check/phase3-package")
if (!dir.exists(stage)) {
  dir.create(stage)
  paths <- list.files("dev/check/phase2-package", all.files = TRUE, no.. = TRUE,
                       full.names = TRUE)
  stopifnot(all(file.copy(paths, stage, recursive = TRUE)))
}
for (path in list.files("dev/check/phase3/R", full.names = TRUE))
  stopifnot(file.copy(path, file.path(stage, "R", basename(path)), overwrite = TRUE))
for (path in list.files("dev/check/phase3/tests", pattern = "\\.R$", full.names = TRUE))
  stopifnot(file.copy(path, file.path(stage, "tests/testthat", basename(path)), overwrite = TRUE))
description <- readLines(file.path(stage, "DESCRIPTION"))
if (!any(grepl("^  tools,", description)))
  description <- sub("^  utils$", "  tools,\n  utils", description)
writeLines(description, file.path(stage, "DESCRIPTION"))
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
devtools::document(stage)
