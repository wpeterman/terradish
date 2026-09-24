robustness_surface <- function() {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
  terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
  names(r) <- "x"
  coords <- terra::xyFromCell(r, c(1, 8, 20, 37, 57, 64))
  conductance_surface(r, coords, saveStack = TRUE)
}

