#' Measurement models
#'
#' Functions that fit a matrix observed genetic distances to a submatrix of the inverse graph Laplacian
#'
#' @details
#' For the ordinary two-component \code{\link{wishart_covariance}} and
#' \code{\link{generalized_wishart}} models with their standard bounds,
#' nuisance fitting first profiles the overall covariance variance analytically
#' and searches the relative graph/nugget contribution on site contrasts.
#' This deterministic initialization avoids false convergence when a warm
#' log nugget has an almost zero derivative. The existing Newton optimizer
#' then refines the result and computes derivatives at the same parameters.
#' Fixed predictive parameters and custom bounds retain their supplied values
#' and constraints. The finite profile search does not guarantee a global
#' optimum for every possible covariance matrix. A nugget approaching zero
#' remains a limiting estimate and can imply weak nuisance identification.
#'
#' @name terradish_measurement_model
#' @seealso \code{\link{leastsquares}}, \code{\link{mlpe}}, \code{\link{generalized_wishart}}
terradish_measurement_model <- NULL

#' Legacy radish measurement-model class alias
#'
#' The legacy \code{"radish_measurement_model"} class name is retained for
#' backward compatibility with older fitted objects and helper code.
#'
#' @name radish_measurement_model
#' @keywords internal
NULL
