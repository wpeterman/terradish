# Simulation and benchmarking utilities. The genetic covariance/distance helpers
# (cov_from_biallelic, cov_from_genetic_data, fst_from_biallelic, dist_from_cov,
# dist_from_biallelic) moved to the landgraph package and are re-exported via
# R/reexports.R; the functions below still call them through terradish's namespace.

#' Simulate covariance responses from a conductance surface
#'
#' Simulates one or more covariance response matrices from a known
#' conductance-resistance relationship, using the same Wishart covariance model
#' fitted by \code{\link{wishart_covariance}}.
#'
#' @param theta Conductance parameters. May be supplied as a numeric vector or
#'   a one-row matrix. If unnamed, values are matched in the order implied by
#'   \code{formula}.
#' @param formula Model formula specifying the conductance covariates. The left
#'   hand side, if supplied, is ignored.
#' @param data A \code{\link{conductance_surface}} object.
#' @param conductance_model Conductance-model factory, such as
#'   \code{\link{loglinear_conductance}}.
#' @param tau Nonnegative scaling applied to the conductance-implied covariance
#'   matrix.
#' @param sigma Nonnegative nugget variance added to the diagonal.
#' @param nu Effective Wishart degrees of freedom: controls the amount of
#'   sampling noise in the simulated covariance: larger \code{nu} gives tighter
#'   draws around \code{Sigma}.  Choose a value that matches what you would use
#'   when fitting the model to real data.  For biallelic SNP-like simulations
#'   use the number of retained polymorphic SNPs. For microsatellite-like
#'   panels, use the number of loci as the primary value;
#'   \eqn{\sum_l (K_l - 1)}, where \eqn{K_l} is the number of observed alleles
#'   at locus \eqn{l}, is a larger sensitivity value. The simulation requires
#'   \code{nu} to be at least the covariance-matrix dimension. See
#'   \code{\link{wishart_covariance}} for how \code{nu} propagates into
#'   estimation (it scales standard errors and model-selection statistics but
#'   not point estimates).
#' @param nsim Number of covariance matrices to simulate.
#' @param seed Optional random seed.
#' @param cores Number of cores passed to \code{\link{terradish_distance}}.
#'
#' @details
#' The helper first computes the conductance-implied covariance matrix
#' \code{E(theta)} from the supplied resistance surface. It then forms
#' \code{Sigma = tau * E(theta) + sigma * I} and simulates
#' \code{S ~ Wishart(nu, Sigma) / nu}. This makes the returned covariance
#' matrices directly compatible with \code{\link{wishart_covariance}}.
#'
#' @return A list containing:
#' \item{covariance}{A covariance matrix if \code{nsim = 1}, otherwise a
#'   three-dimensional array of simulated covariance matrices.}
#' \item{Sigma}{The covariance matrix used as the Wishart scale matrix after
#'   applying \code{tau} and \code{sigma}.}
#' \item{E}{The conductance-implied covariance matrix returned by
#'   \code{\link{terradish_distance}(covariance = TRUE)}.}
#' \item{theta}{The validated conductance parameter matrix used for simulation.}
#' \item{tau}{The supplied \code{tau}.}
#' \item{sigma}{The supplied \code{sigma}.}
#' \item{nu}{The supplied \code{nu}.}
#'
#' @examples
#' r1 <- terra::rast(nrows = 3, ncols = 3, vals = 1:9)
#' r2 <- terra::rast(nrows = 3, ncols = 3, vals = c(9, 1, 4, 3, 8, 2, 5, 7, 6))
#' covariates <- c(r1, r2)
#' names(covariates) <- c("x1", "x2")
#' pts <- terra::vect(matrix(c(0.5, 0.5,
#'                            1.5, 1.5,
#'                            2.5, 2.5), ncol = 2, byrow = TRUE),
#'                    type = "points")
#' surface <- conductance_surface(covariates, pts, directions = 4)
#' sim <- simulate_covariance_response(
#'   theta = c(x1 = 0.3, x2 = -0.2),
#'   formula = ~ x1 + x2,
#'   data = surface,
#'   tau = 0.8,
#'   sigma = 0.5,
#'   nu = 20,
#'   seed = 1
#' )
#' sim$covariance
#'
#' @export
simulate_covariance_response <- function(theta,
                                         formula,
                                         data,
                                         conductance_model = loglinear_conductance,
                                         tau = 1,
                                         sigma = 0,
                                         nu,
                                         nsim = 1,
                                         seed = NULL,
                                         cores = 1L)
{
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(is.numeric(tau), length(tau) == 1L, is.finite(tau), tau >= 0)
  stopifnot(is.numeric(sigma), length(sigma) == 1L, is.finite(sigma), sigma >= 0)
  stopifnot(is.numeric(nu), length(nu) == 1L, is.finite(nu), nu > 0)
  stopifnot(is.numeric(nsim), length(nsim) == 1L, nsim >= 1)
  stopifnot(length(cores) == 1L, is.numeric(cores), cores >= 1)

  terms_obj <- terms(formula)
  if (!length(attr(terms_obj, "term.labels")))
    stop("`formula` must include at least one conductance covariate.", call. = FALSE)
  model_formula <- reformulate(attr(terms_obj, "term.labels"))

  conductance_model_obj <- conductance_model(model_formula, data$x)
  default <- attr(conductance_model_obj, "default")

  if (is.null(dim(theta)))
    theta <- matrix(theta, nrow = 1L)
  stopifnot(is.matrix(theta))
  stopifnot(nrow(theta) == 1L)
  stopifnot(ncol(theta) == length(default))
  theta <- .validate_theta_grid(theta, names(default))

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (old_seed_exists)
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed))
    set.seed(seed)

  E <- terradish_distance(
    theta = theta,
    formula = model_formula,
    data = data,
    conductance_model = conductance_model,
    conductance = TRUE,
    covariance = TRUE,
    cores = as.integer(cores)
  )$covariance[, , 1]

  Sigma <- tau * E + sigma * diag(nrow(E))
  if (nu < nrow(Sigma))
    stop("`nu` must be at least the covariance-matrix dimension for Wishart simulation.",
         call. = FALSE)
  covariance <- array(NA_real_, dim = c(nrow(E), ncol(E), as.integer(nsim)))
  for (i in seq_len(as.integer(nsim)))
    covariance[, , i] <- rWishart(1, df = nu, Sigma = Sigma)[, , 1] / nu

  list(
    covariance = if (as.integer(nsim) == 1L) covariance[, , 1] else covariance,
    Sigma = Sigma,
    E = E,
    theta = theta,
    tau = tau,
    sigma = sigma,
    nu = nu
  )
}
