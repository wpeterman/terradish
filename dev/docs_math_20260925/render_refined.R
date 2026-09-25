root <- getwd()
.libPaths(c(file.path(root, "dev/check/phase2-library"), .libPaths()))
pkgload::load_all(file.path(root, "dev/check/docs-math-20260925-package"), quiet=TRUE, export_all=FALSE)
out <- file.path(root, "dev/check/docs-math-20260925-refined")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
results <- list()
for (name in c("model-comparison", "gaussian-scale-optimization")) {
 f <- paste0(name, ".Rmd")
 file.copy(file.path(root, "vignettes", f), file.path(out,f), overwrite=FALSE)
 results[[name]] <- rmarkdown::render(file.path(out,f), output_dir=out,
   envir=new.env(parent=globalenv()),quiet=TRUE)
 cat("Rendered final",name,"\n")
 saveRDS(results,file.path(root,"dev/docs_math_20260925/refined_results.rds"))
}
