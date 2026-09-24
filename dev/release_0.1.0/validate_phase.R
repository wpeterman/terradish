# Usage: Rscript --vanilla dev/release_0.1.0/validate_phase.R base test|check
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("test", "check"))
phase <- args[1]
action <- args[2]
output <- file.path("dev/release_0.1.0", paste0(phase, "_", action, ".rds"))
started <- Sys.time()
result <- tryCatch({
  if (action == "test") {
    tests <- devtools::test(reporter = "summary", stop_on_failure = FALSE)
    as.data.frame(tests)
  } else {
    # R CMD build copies .git before applying .Rbuildignore on Windows.
    # Stage tracked working-tree files to avoid long internal Git ref paths.
    stage <- file.path(tempdir(), "terradish")
    dir.create(stage)
    tracked <- system2("git", c("--no-optional-locks", "ls-files",
      "--cached", "--others", "--exclude-standard"), stdout = TRUE)
    tracked <- unique(tracked[file.exists(tracked)])
    tracked <- tracked[!grepl("^(dev/|\\.git/)", tracked)]
    for (path in tracked) {
      dir.create(dirname(file.path(stage, path)), recursive = TRUE,
                 showWarnings = FALSE)
      stopifnot(file.copy(path, file.path(stage, path), overwrite = TRUE))
    }
    devtools::check(pkg = stage, args = "--no-manual", document = FALSE,
      error_on = "never", check_dir = file.path(getwd(), "dev/check", phase))
  }
}, error = function(e) list(error = conditionMessage(e)))
saveRDS(list(started = started, finished = Sys.time(), result = result,
             session = sessionInfo()), output)
print(result)
