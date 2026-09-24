suppressMessages({library(terra); library(multiScaleR)})
set.seed(1)
mk <- function(res) { r <- rast(nrows=30, ncols=30, xmin=0, xmax=30*res, ymin=0, ymax=30*res); values(r) <- 0; r[15,15] <- 1; names(r) <- "a"; r }
for (res in c(1, 100)) {
  r <- mk(res)
  s <- kernel_scale.raster(r, sigma = 2, kernel = "gaussian", verbose = FALSE)
  m <- as.matrix(s, wide=TRUE)
  cat("res =", res, " sigma=2: center value", round(m[15,15],4), " value 2 cells away", signif(m[15,17],4), "\n")
  s2 <- kernel_scale.raster(r, sigma = 2*res, kernel = "gaussian", verbose = FALSE)
  m2 <- as.matrix(s2, wide=TRUE)
  cat("res =", res, " sigma=2*res: center value", round(m2[15,15],4), " value 2 cells away", signif(m2[15,17],4), "\n")
}
data(melip, package="terradish")
fc <- unwrap(melip.forestcover); print(res(fc)); print(dim(fc)); print(ext(fc))
fc3 <- aggregate(fc, 3, na.rm=TRUE); print(res(fc3)); print(dim(fc3))
