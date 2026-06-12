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
#'   use the number of retained polymorphic SNPs.  For microsatellite-like
#'   panels, use the number of loci as the conservative default; the true
#'   effective degrees of freedom is likely between the locus count \eqn{L} and
#'   \eqn{\sum_l (K_l - 1)}, where \eqn{K_l} is the number of observed alleles
#'   at locus \eqn{l}.  See \code{\link{wishart_covariance}} for how \code{nu}
#'   propagates into estimation (it scales standard errors and model-selection
#'   statistics but not point estimates).
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
  covariance <- array(NA_real_, dim = c(nrow(E), ncol(E), as.integer(nsim)))
  for (i in seq_len(as.integer(nsim)))
    covariance[, , i] <- rWishart(1, df = as.integer(nu), Sigma = Sigma)[, , 1] / as.integer(nu)

  list(
    covariance = if (as.integer(nsim) == 1L) covariance[, , 1] else covariance,
    Sigma = Sigma,
    E = E,
    theta = theta,
    tau = tau,
    sigma = sigma,
    nu = as.integer(nu)
  )
}

# Internal simulation helper for benchmarking generalized Wishart fits.
wishart_simulate_distance <- function(seed, S, nu)
{
  set.seed(seed)
  dist_from_cov(solve(rWishart(1, nu, S)[,,1]))
}

wishart_simulate_experiment <- function(seed, N, P, K, nu, neval, timingOnly=FALSE)
{
  covariates <- list()
  for(k in 1:K)
  {
    sim <- matrix(.randomfields_rfsimulate(.randomfields_rmexp(scale = 10),
                                           expand.grid(x = 1:N, y = 1:N),
                                           seed = seed)@data[[1]], N, N)
    covariates[[paste0("var", k)]] <- rast(nrows = N, ncols = N,
                                           xmin = 0, xmax = N, ymin = 0, ymax = N)
    values(covariates[[paste0("var", k)]]) <- as.vector(sim)
    covariates[[paste0("var", k)]] <- covariates[[paste0("var", k)]] /
      max(values(covariates[[paste0("var", k)]])[,1], na.rm = TRUE)
  }
  covariates <- do.call(c, covariates)
  names(covariates) <- paste0("var", seq_len(K))

  set.seed(seed)
  demes <- sample(1:(N^2), P)
  coords <- xyFromCell(covariates[[1]], demes)

  surf <- conductance_surface(covariates, coords, directions = 4, saveStack = FALSE)

  beta <- matrix(rnorm(K), 1, K)
  E <- terradish_distance(loglinear_conductance, surf, beta, covariance = TRUE)$covariance[,,1]
  S <- wishart_simulate_distance(seed=seed, nu=nu, S=E)
  S <- S/(max(S)*10)

  # functions
  nograd <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = FALSE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  wgrad <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = FALSE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }
  whess <- function(par)
  {
    terradish_algorithm(f = loglinear_conductance, g = leastsquares, s = surf, S = S, theta = c(par), 
                     gradient = TRUE,
                     hessian = TRUE,
                     partial = FALSE, 
                     nonnegative = TRUE)
  }

  timings_whess <- c()
  timings_wgrad <- c()
  timings_nograd <- c()
  for(i in 1:neval)
  {
    timings_whess <- rbind(timings_whess, system.time(ll <- whess(rep(0,K))))
    timings_wgrad <- rbind(timings_wgrad, system.time(ll <- wgrad(rep(0,K))))
    timings_nograd <- rbind(timings_nograd, system.time(ll <- nograd(rep(0,K))))
  }

  #optimization
  if(!timingOnly){
  opt_newton <- list(fcall = NA, loglik = NA, fit = list(boundary = NA))
  fcall <- 0L
  bqfn <- function(par)
  {
    fcall <<- fcall + 1L
    nograd(par)$objective
  }
  opt_bobyqa <- .nloptr_bobyqa(rep(0, K), bqfn)
  opt_bobyqa$fcall <- fcall
} else {
  opt_bobyqa = list(fcall=NA, value=NA)
  opt_newton = list(fcall=NA, loglik=NA, fit=list(boundary=NA))
  }


  list(timings=list(nograd=timings_nograd, wgrad=timings_wgrad, whess=timings_whess), opt_newton=opt_newton, 
       #opt_lbfgs=opt_lbfgs, 
       opt_bobyqa=opt_bobyqa, seed=seed, beta=beta, K=K, N=N, P=P)
}

run_benchmarks <- function(K = c(1,2,4,8,16,32), N = c(100), P = c(30), reps=5, timingOnly=timingOnly)
{
  set.seed(1)
  seeds <- sample.int(100000, length(K)*length(N)*length(P)*reps)
  z <- 0
  out <- list()
  for(p in P)
  {
    out[[as.character(p)]] <- list()
    for(n in N)
    {
      out[[as.character(p)]][[as.character(n)]] <- list()
      for(k in K)
      {
        out[[as.character(p)]][[as.character(n)]][[as.character(k)]] <- list()
        for(rep in 1:reps)
        {
          z <- z + 1
          try({
            out[[as.character(p)]][[as.character(n)]][[as.character(k)]][[as.character(rep)]] <- 
              wishart_simulate_experiment(seeds[z], n, p, k, 100, 1, timingOnly=timingOnly)
          })
        }
      }
    }
  }
  out
}
#aahh <- run_benchmarks(c(1,2,4,8,16,32), c(200), c(30), 10, timingOnly=TRUE)

extract_benchmark_timing <- function(benchmarks)
{
  out <- c()
  for(p in names(benchmarks))
    for(n in names(benchmarks[[p]]))
      for(k in names(benchmarks[[p]][[n]]))
        for(r in names(benchmarks[[p]][[n]][[k]]))
        {
          tmp <- benchmarks[[p]][[n]][[k]][[r]]
          oo <- data.frame(P=as.numeric(p), N=as.numeric(n), K=as.numeric(k), rep=as.numeric(r),
                           timing_nograd = mean(tmp$timings$nograd[,1]),
                           timing_wgrad = mean(tmp$timings$wgrad[,1]),
                           timing_whess = mean(tmp$timings$whess[,1]),
                           newton_eval  = tmp$opt_newton$fcall,
                           newton_ll    = tmp$opt_newton$loglik,
                           newton_boundary = tmp$opt_newton$fit$boundary,
                           #lbfgs_eval  = tmp$opt_lbfgs$fcall,
                           bobyqa_eval  = tmp$opt_bobyqa$fcall,
                           bobyqa_ll    = tmp$opt_bobyqa$value
                           )
          out <- rbind(out, oo)
        }
  out
}
