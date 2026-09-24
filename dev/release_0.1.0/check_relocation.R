# Phase 1a must preserve every namespace function's formals, body and class.
# Ignore source-reference locations because relocation intentionally changes them.
args <- commandArgs(trailingOnly = TRUE)
pkgload::load_all()
ns <- asNamespace("terradish")
snapshot <- list()
for (name in ls(ns, all.names = TRUE)) {
  object <- get(name, ns)
  if (is.function(object)) snapshot[[name]] <- list(
    formals = paste(deparse(formals(object)), collapse = "\n"),
    body = paste(deparse(body(object)), collapse = "\n"),
    class = class(object))
}
path <- "dev/release_0.1.0/base_functions.rds"
if (identical(args, "record")) {
  if (file.exists(path)) stop("Original functions already recorded.")
  saveRDS(snapshot, path)
  cat(length(snapshot), "functions recorded\n")
} else {
  before <- readRDS(path)
  stopifnot(identical(names(before), names(snapshot)))
  changed <- names(before)[!vapply(names(before), function(name)
    identical(before[[name]], snapshot[[name]]), logical(1))]
  print(changed)
  stopifnot(length(changed) == 0L)
}
