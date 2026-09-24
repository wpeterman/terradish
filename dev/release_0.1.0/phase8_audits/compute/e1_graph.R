suppressMessages({library(terradish); library(terra); library(Matrix)})
set.seed(1)
# ---- E1a: scale_to_0_1 per-layer claim ----
r <- c(rast(nrows=2,ncols=2,vals=c(1,2,3,4)), rast(nrows=2,ncols=2,vals=c(100,200,300,400)))
cat("E1a scale_to_0_1 on 2-layer raster (docs: per-layer):\n")
print(values(scale_to_0_1(r)))


stopifnot(isTRUE(all.equal(values(scale_to_0_1(r))[, 1], values(scale_to_0_1(r))[, 2])))
