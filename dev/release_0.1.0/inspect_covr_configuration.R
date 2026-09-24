cat(deparse(covr::package_coverage), sep = "\n")
cat("\nCoverage options\n")
print(options()[grep("covr", names(options()))])
