.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
devtools::document()
source("dev/release_0.1.0/check_globals.R")
