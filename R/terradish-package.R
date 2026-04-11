#' terradish package
#'
#' Gradient-based optimization of conductance and resistance surfaces.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom Matrix Cholesky Diagonal forceSymmetric rowSums solve sparseMatrix t update
#' @importFrom methods new setRefClass
#' @importFrom multiScaleR kernel_scale.raster
#' @importFrom nlme gls
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply stopCluster
#' @importFrom stats AIC D anova as.dist as.formula coef cov2cor delete.response dist
#' @importFrom stats dnorm fft fitted lm logLik model.matrix optimize pchisq plogis
#' @importFrom stats pnorm printCoefmat prcomp qlogis qnorm reformulate residuals rnorm
#' @importFrom stats rWishart sd setNames sigma simulate terms
#' @importFrom terra adjacent cellFromXY geom global is.factor levels ncell nlyr patches rast
#' @importFrom terra ext extract is.lonlat res rowColFromCell unwrap values values<- xyFromCell
#' @importFrom utils globalVariables modifyList write.csv
#' @importFrom grDevices terrain.colors
#' @useDynLib terradish, .registration = TRUE
"_PACKAGE"

globalVariables(c(
  "distance", "distance_lower", "distance_upper", "est", "label",
  "label_y", "observed", "upper", "weight", "x", "y"
))
