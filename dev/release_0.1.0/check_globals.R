# Compare unresolved globals to the original namespace, including existing
# codetools false positives from functions using nonstandard evaluation.
args <- commandArgs(trailingOnly = TRUE)
pkgload::load_all()
ns <- asNamespace("terradish")
missing <- list()
for (name in ls(ns, all.names = TRUE)) {
  object <- get(name, ns)
  if (is.function(object)) {
    globals <- codetools::findGlobals(object, merge = FALSE)$functions
    absent <- globals[!vapply(globals, exists, logical(1), envir = ns)]
    if (length(absent)) missing[[name]] <- sort(absent)
  }
}
reference <- "dev/release_0.1.0/base_globals.rds"
if (identical(args, "record")) {
  if (file.exists(reference)) stop("Original globals already recorded.")
  saveRDS(missing, reference)
  print(missing)
} else {
  original <- readRDS(reference)
  added <- lapply(names(missing), function(name)
    setdiff(missing[[name]], original[[name]]))
  names(added) <- names(missing)
  added <- added[lengths(added) > 0]
  print(added)
  stopifnot(length(added) == 0L)
}
