.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
pkgload::load_all("dev/check/phase3-package")
code <- readLines("dev/release_0.1.0/core_baseline.R")
code <- code[code != "pkgload::load_all()"]
code <- sub('model = "generalized_wishart")',
            'model = "generalized_wishart", transform = "sqdiff")', code, fixed = TRUE)
code <- sub('boundary = fit$fit$boundary)',
            'boundary = fit$fit$boundary, convergence = fit$convergence)', code, fixed = TRUE)
eval(parse(text = code), envir = globalenv())
previous <- readRDS("dev/release_0.1.0/core_phase2.rds")$models
for (name in names(records$models)) {
  current <- records$models[[name]]
  old <- previous[[name]]
  stopifnot(current$convergence$code == 0,
    abs(current$logLik - old$logLik) / max(abs(old$logLik), .Machine$double.eps) <= 1e-6,
    max(abs(current$theta - old$theta)) <= 1e-5,
    max(abs(current$phi - old$phi)) <= 1e-5,
    max(abs(current$se - old$se)) <= 1e-5)
}
cat("All nine complete-package Phase 3 baselines converge and agree within tolerance.\n")
