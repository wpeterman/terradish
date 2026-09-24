.libPaths(c("dev/check/phase2-library", .libPaths()))
pkg <- "dev/check/phase6-package"
devtools::document(pkg)
spelling <- spelling::spell_check_package(pkg, vignettes = FALSE)
print(spelling)
saveRDS(spelling, "dev/release_0.1.0/phase6_spelling.rds")
