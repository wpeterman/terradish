.libPaths(c("dev/check/phase2-library", .libPaths()))
devtools::document("dev/check/phase5-package")
result <- devtools::test("dev/check/phase5-package",
  filter = "core-robustness|new-surface-prediction|spline-reuse|gaussian-alignment|nu-fit-power|optimizer-stopping|cv-contracts|spatial-cv|terra-workflow",
  reporter = "summary", stop_on_failure = FALSE)
saveRDS(result, "dev/release_0.1.0/phase5_focused.rds")
print(as.data.frame(result)[, c("test", "passed", "failed", "error", "warning")])
