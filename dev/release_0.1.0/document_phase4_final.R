.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
devtools::document("dev/check/phase4-package")
