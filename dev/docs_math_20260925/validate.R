root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkg <- file.path(root, "dev/check/docs-math-20260925-package")
out <- file.path(root, "dev/check/docs-math-20260925-render")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
devtools::document(pkg, quiet = TRUE)
# Assert all implementation syntax is unchanged from the checked release.
old <- file.path(root, "dev/check/phase8-final-package")
for (f in list.files(file.path(pkg, "R"), pattern = "[.]R$")) {
 stopifnot(identical(parse(file.path(pkg, "R", f), keep.source = FALSE),
                     parse(file.path(old, "R", f), keep.source = FALSE)))
}
cat("All executable R syntax unchanged.\n")
file.copy(list.files(file.path(pkg, "man"), full.names = TRUE),
          file.path(root, "man"), overwrite = TRUE)
helpout <- file.path(out, "help")
dir.create(helpout, showWarnings = FALSE)
for (f in list.files(file.path(pkg, "man"), pattern = "[.]Rd$", full.names = TRUE)) {
 rd <- tools::parse_Rd(f)
 tools::checkRd(rd)
 tools::Rd2HTML(rd, out = file.path(helpout, sub("[.]Rd$", ".html", basename(f))),
                package = "terradish", dynamic = FALSE)
}
cat("All help pages parsed and rendered.\n")
results <- list()
files <- list.files(file.path(pkg, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE)
files <- files[order(basename(files) != "model-comparison.Rmd", basename(files))]
for (f in files) {
 name <- basename(f); cat("RENDER", name, "\n")
 started <- Sys.time()
 result <- tryCatch(rmarkdown::render(f, output_dir = out,
   envir = new.env(parent = globalenv()), quiet = TRUE), error = function(e) list(error = conditionMessage(e)))
 results[[name]] <- list(started = started, finished = Sys.time(), result = result)
 saveRDS(results, file.path(root, "dev/docs_math_20260925/render_results.rds"))
 print(results[[name]])
}
stopifnot(all(vapply(results, function(x) is.character(x$result), logical(1))))
cat("All eight vignettes rendered successfully.\n")
