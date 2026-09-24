root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkg <- file.path(root, "dev/check/phase7-package")
devtools::document(pkg, quiet = TRUE)
spelling <- spelling::spell_check_package(pkg, vignettes = TRUE)
print(spelling)
saveRDS(spelling, file.path(root, "dev/release_0.1.0/phase7_spelling.rds"))
# R changes since C6 should still be documentation-only (currently identical).
old <- file.path(root, "dev/check/phase6-package/R")
same <- vapply(list.files(old, full.names = TRUE), function(f) {
  identical(parse(f, keep.source = FALSE),
    parse(file.path(pkg, "R", basename(f)), keep.source = FALSE))
}, logical(1))
stopifnot(all(same))
cat(length(same), "R files retain their C6 executable expressions.\n")
