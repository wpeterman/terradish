# Run the baseline against the frozen, complete Phase 2 package candidate.
.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
pkgload::load_all("dev/check/phase2-package")
code <- readLines("dev/release_0.1.0/core_baseline.R")
code <- code[code != "pkgload::load_all()"]
code <- sub('model = "generalized_wishart")',
            'model = "generalized_wishart", transform = "sqdiff")', code, fixed = TRUE)
eval(parse(text = code), envir = globalenv())
