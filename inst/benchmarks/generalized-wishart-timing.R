# Benchmark harness for the generalized-Wishart fitting path.
#
# Moved out of R/ because nothing in the package calls it: it is a developer
# timing script, it sets the global RNG seed without restoring it, and it
# depends on RandomFields, which is archived on CRAN. The terradish_algorithm()
# and terradish_distance() calls below use an older argument order and need
# updating before the harness will run again.
#
# Usage: source() this file with terradish, terra, nloptr, and a Gaussian
# random field generator attached.

.randomfields_rmexp <- function(...)
  terradish:::.suggested_export("RandomFields", "RMexp",
                                context = "benchmark simulation helpers")(...)

.randomfields_rfsimulate <- function(...)
  terradish:::.suggested_export("RandomFields", "RFsimulate",
                                context = "benchmark simulation helpers")(...)

.nloptr_bobyqa <- function(...)
  terradish:::.suggested_export("nloptr", "bobyqa",
                                context = "benchmark optimization helpers")(...)

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
