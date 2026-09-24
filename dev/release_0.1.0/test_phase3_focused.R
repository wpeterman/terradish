.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
result <- devtools::test("dev/check/phase3-package",
  filter = "core-inference|scale-inference|optimizer-stopping|warm-starts",
  reporter = "summary", stop_on_failure = FALSE)
saveRDS(as.data.frame(result), "dev/release_0.1.0/phase3_focused_test.rds")
