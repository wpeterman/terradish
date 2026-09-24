.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
pkgload::load_all("dev/check/phase3-package")
scratch <- new.env(parent = asNamespace("terradish"))
sys.source("dev/check/phase4/R/spatial_cv.R", scratch)
result <- testthat::test_file("dev/check/phase4/tests/test-cv-contracts.R",
  env = scratch, reporter = "summary", stop_on_failure = FALSE)
saveRDS(as.data.frame(result), "dev/release_0.1.0/phase4_prototype_test.rds")
