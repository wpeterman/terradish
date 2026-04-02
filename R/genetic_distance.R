#' Fst from allele counts
#'
#' Calculates Fst using the estimator of Bhatia et al.
#' from counts of the derived allele across biallelic markers
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @return A matrix containing pairwise Fst
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' fst_from_biallelic(Y, N)
#'
#' @export

fst_from_biallelic <- function(Y, N)
{
  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 1))
    stop("All sample sizes in `N` must exceed 1")

  Fr <- Y/N
  f2 <- apply(apply(Fr, 2, function(x) outer(x,x,"-")^2), 1, mean, na.rm=TRUE)
  cornum <- apply(Fr*(1-Fr)/(N-1), 1, mean, na.rm=TRUE)
  corden <- apply(Fr*(1-Fr)*N/(N-1), 1, mean, na.rm=TRUE)
  fst <- (f2 - c(outer(cornum, cornum, "+")))/(f2 - c(outer(cornum, cornum, "+")) + c(outer(corden, corden, "+")))
  fst <- matrix(fst, nrow(Y), nrow(Y))
  diag(fst) <- 0
  fst
}

#' Covariance from allele counts
#'
#' Calculates covariance from counts of the derived allele across biallelic markers,
#' using normalized allele frequencies
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @details
#' Let \code{p[l] = sum(Y[, l]) / sum(N[, l])} be the pooled derived-allele
#' frequency at locus \code{l}. The function returns
#' \code{Z \%*\% t(Z) / L}, where
#' \code{Z[i, l] = (Y[i, l] - N[i, l] * p[l]) / sqrt(N[i, l] * p[l] * (1 - p[l]))}
#' and \code{L} is the number of loci.
#'
#' @return A matrix containing pairwise covariance
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' cov_from_biallelic(Y, N)
#'
#' @export

cov_from_biallelic <- function(Y, N)
{
  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (anyNA(Y) || anyNA(N))
    stop("missing values are not currently supported")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")
  if (any(N <= 0))
    stop("All sample sizes in `N` must be positive")

  Fr <- colSums(Y) / colSums(N)
  if (any(Fr <= 0 | Fr >= 1))
    stop("Each locus must have pooled allele frequency strictly between 0 and 1")
  Y  <- (Y - N*Fr) / sqrt(N * Fr * (1-Fr))

  Y %*% t(Y) / ncol(Y)
}

#' Distance matrix from covariance matrix
#'
#' Returns the squared-distance matrix associated with a given covariance matrix
#'
#' @param Cov A covariance matrix (does not need to be full-rank)
#'
#' @details
#' The returned squared-distance matrix is
#' \code{D[i, j] = Cov[i, i] + Cov[j, j] - 2 * Cov[i, j]}.
#'
#' @return A distance matrix of the same dimension as the input
#'
#' @examples
#' Cov <- matrix(c(2.0, 1.2, 0.8,
#'                 1.2, 1.5, 0.7,
#'                 0.8, 0.7, 1.1), nrow = 3, byrow = TRUE)
#' dist_from_cov(Cov)
#'
#' @export

dist_from_cov <- function(Cov)
{
  ones <- matrix(1, nrow(Cov), 1)
  diag(Cov) %*% t(ones) + ones %*% t(diag(Cov)) - 2*Cov
}

#' Distance from allele counts
#'
#' Calculates genetic distance from counts of the derived allele across biallelic markers,
#' using normalized allele frequencies
#'
#' @param Y A matrix containing allele counts, of dimension (number of demes) x (number of loci)
#' @param N A matrix containing the number of sampled haplotypes, of dimension (number of demes) x (number of loci)
#'
#' @details
#' This is a convenience wrapper for
#' \code{dist_from_cov(cov_from_biallelic(Y, N))}.
#'
#' @return A matrix containing pairwise distance
#'
#' @examples
#' Y <- matrix(c(2, 1, 0,
#'               1, 1, 1,
#'               0, 1, 2), nrow = 3, byrow = TRUE)
#' N <- matrix(2, nrow = 3, ncol = 3)
#' dist_from_biallelic(Y, N)
#'
#' @export

dist_from_biallelic <- function(Y, N)
{
  dist_from_cov(cov_from_biallelic(Y, N))
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
