#' terradish package
#'
#' Gradient-based optimization of conductance and resistance surfaces.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom Matrix Cholesky Diagonal forceSymmetric rowSums solve sparseMatrix t update
#' @importFrom methods new setRefClass
#' @importFrom nlme gls
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply stopCluster
#' @importFrom stats AIC D anova as.dist as.formula coef cov2cor delete.response dist fitted
#' @importFrom stats logLik model.matrix pchisq plogis pnorm printCoefmat prcomp qlogis
#' @importFrom stats qnorm reformulate residuals rnorm rWishart sigma simulate terms
#' @importFrom terra adjacent cellFromXY geom global is.factor levels ncell nlyr patches rast
#' @importFrom terra unwrap values values<- xyFromCell
#' @importFrom utils globalVariables modifyList write.csv
#' @useDynLib terradish, .registration = TRUE
"_PACKAGE"

utils::globalVariables(c("est", "observed", "upper", "x", "y"))
