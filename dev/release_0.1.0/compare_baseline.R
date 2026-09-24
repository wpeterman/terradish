# Compare a new run to the frozen 0.0.47 baseline. Optional remaining arguments
# select models expected to be unchanged in the current phase.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
base <- readRDS("dev/release_0.1.0/core_baseline_0.0.47.rds")$models
new <- readRDS(args[1])$models
models <- if (length(args) > 1L) args[-1] else names(base)
ok <- TRUE
for (name in models) {
  a <- base[[name]]
  b <- new[[name]]
  stopifnot(!is.null(a), !is.null(b), is.null(a$error), is.null(b$error))
  difference <- c(logLik_relative = abs(a$logLik - b$logLik) /
                   max(abs(a$logLik), .Machine$double.eps),
                  theta = max(abs(a$theta - b$theta)),
                  phi = max(abs(a$phi - b$phi)),
                  se = max(abs(a$se - b$se)))
  passed <- all(is.finite(difference)) &&
    difference[1] <= 1e-6 && all(difference[-1] <= 1e-5)
  cat(name, if (passed) "PASS" else "FAIL", "\n")
  print(difference)
  ok <- ok && passed
}
stopifnot(ok)
