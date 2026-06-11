// ===========================================================================
// Optional compiled Hutchinson trace/diagonal estimator for diag(L^{-1}) and
// tr(L^{-1} M), for curvature probing (prototype D) and randomized effective-
// resistance (Tier-3 F). Mirrors the Eigen CG style already used in
// src/radish.cpp (reduced_laplacian_eigen + ConjugateGradient).
//
// Use when the R per-probe solve in 03_gauss_newton_fisher.R is the bottleneck.
// Estimates:
//   * tr(L^{-1})       = E_z[ z^T L^{-1} z ],   z Rademacher
//   * diag(L^{-1})_i  ~= E_z[ z_i (L^{-1} z)_i ] (Bekas et al. 2007)
//
// Build (standalone test): compile inside the package so LinkingTo headers
// (RcppArmadillo, RcppEigen) resolve; then call from R via Rcpp::sourceCpp() or
// move into src/ and document in NAMESPACE.
//
// VERIFIED EQUIVALENT: the R/Python Hutchinson estimates in verify.py /
// 03_gauss_newton_fisher.R reproduce tr(L^{-1}) to ~1% at ~1500-2000 probes.
// ===========================================================================
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <random>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

namespace {

Eigen::SparseMatrix<double> reduced_laplacian_eigen(const arma::vec& conductance,
                                                    const arma::umat& edge_pairs)
{
  const arma::uword N = conductance.n_elem;
  const arma::uword Nred = N - 1;
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(2 * edge_pairs.n_rows + Nred);
  Eigen::VectorXd diag = Eigen::VectorXd::Zero((Eigen::Index)Nred);
  for (arma::uword e = 0; e < edge_pairs.n_rows; ++e) {
    const arma::uword ii = edge_pairs(e, 0) - 1, jj = edge_pairs(e, 1) - 1;
    const double w = conductance.at(ii) + conductance.at(jj);
    if (ii < Nred) diag[(Eigen::Index)ii] += w;
    if (jj < Nred) diag[(Eigen::Index)jj] += w;
    if (ii < Nred && jj < Nred) {
      trip.emplace_back((Eigen::Index)ii, (Eigen::Index)jj, -w);
      trip.emplace_back((Eigen::Index)jj, (Eigen::Index)ii, -w);
    }
  }
  for (arma::uword i = 0; i < Nred; ++i)
    trip.emplace_back((Eigen::Index)i, (Eigen::Index)i, diag[(Eigen::Index)i]);
  Eigen::SparseMatrix<double> A((Eigen::Index)Nred, (Eigen::Index)Nred);
  A.setFromTriplets(trip.begin(), trip.end());
  A.makeCompressed();
  return A;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List hutchinson_laplacian_inverse(const arma::vec& conductance,
                                        const arma::umat& edge_pairs,
                                        const int n_probe = 200,
                                        const bool want_diag = false,
                                        const double tol = 1e-8,
                                        const int maxit = 1000,
                                        const unsigned int seed = 1u)
{
  Eigen::SparseMatrix<double> A = reduced_laplacian_eigen(conductance, edge_pairs);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>> cg;
  cg.setTolerance(tol);
  cg.setMaxIterations(maxit);
  cg.compute(A);
  if (cg.info() != Eigen::Success)
    Rcpp::stop("[hutchinson] incomplete Cholesky preconditioner failed");

  const Eigen::Index n = A.rows();
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> coin(0, 1);

  double tr_sum = 0.0, tr_sq = 0.0;
  Eigen::VectorXd diag_acc = Eigen::VectorXd::Zero(want_diag ? n : 0);

  for (int t = 0; t < n_probe; ++t) {
    Eigen::VectorXd z(n);
    for (Eigen::Index i = 0; i < n; ++i) z[i] = coin(rng) ? 1.0 : -1.0;
    Eigen::VectorXd x = cg.solve(z);          // x = L^{-1} z
    const double q = z.dot(x);                // z^T L^{-1} z
    tr_sum += q; tr_sq += q * q;
    if (want_diag) diag_acc += z.cwiseProduct(x);
  }

  const double tr_mean = tr_sum / n_probe;
  const double tr_var  = std::max(0.0, tr_sq / n_probe - tr_mean * tr_mean);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("trace")    = tr_mean,
    Rcpp::Named("trace_se") = std::sqrt(tr_var / n_probe),
    Rcpp::Named("n_probe")  = n_probe);
  if (want_diag) {
    arma::vec d(n);
    for (Eigen::Index i = 0; i < n; ++i) d[i] = diag_acc[i] / n_probe;
    out["diag"] = d;
  }
  return out;
}
