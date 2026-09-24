old <- "dev/check/phase5-package/R"
new <- "dev/check/phase6-package/R"
files <- list.files(old, pattern = "[.]R$")
stopifnot(identical(files, list.files(new, pattern = "[.]R$")))
for (name in files) {
  a <- parse(file.path(old, name), keep.source = FALSE)
  b <- parse(file.path(new, name), keep.source = FALSE)
  if (!identical(a, b)) stop("Executable R code differs: ", name)
}
cat(length(files), "R files have identical executable expressions; Phase 6 changes documentation only.\n")
