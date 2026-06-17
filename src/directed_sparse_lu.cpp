#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/SparseLU>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

namespace {

typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> SpMat;
typedef Eigen::Triplet<double> Triplet;

class DirectedSparseLUFactor {
public:
  int n;
  int absorber;
  std::vector<int> idx;
  std::vector<int> rid;
  Eigen::VectorXd hfull;
  Eigen::SparseLU<SpMat> solver;

  DirectedSparseLUFactor(int n_, int absorber_)
    : n(n_), absorber(absorber_), rid(n_, -1), hfull(Eigen::VectorXd::Zero(n_))
  {
    idx.reserve(n_ - 1);
    for (int v = 0; v < n_; ++v) {
      if (v == absorber_) continue;
      rid[v] = static_cast<int>(idx.size());
      idx.push_back(v);
    }
  }
};

void validate_directed_inputs(const arma::umat& edges,
                              const arma::vec& rate,
                              const arma::uvec& focal,
                              const int n)
{
  if (n < 2) Rcpp::stop("`n` must be at least 2.");
  if (edges.n_cols != 2) Rcpp::stop("`edges` must have two columns.");
  if (edges.n_rows != rate.n_elem)
    Rcpp::stop("`rate` must have one entry per directed edge.");
  for (arma::uword k = 0; k < edges.n_rows; ++k) {
    if (edges(k, 0) < 1 || edges(k, 0) > static_cast<unsigned int>(n) ||
        edges(k, 1) < 1 || edges(k, 1) > static_cast<unsigned int>(n)) {
      Rcpp::stop("`edges` contains a vertex outside 1:n.");
    }
    if (!R_finite(rate[k]) || rate[k] <= 0.0)
      Rcpp::stop("`rate` must contain finite positive values.");
  }
  for (arma::uword k = 0; k < focal.n_elem; ++k) {
    if (focal[k] < 1 || focal[k] > static_cast<unsigned int>(n))
      Rcpp::stop("`focal` contains a vertex outside 1:n.");
  }
}

DirectedSparseLUFactor* factor_one_absorber(const arma::umat& edges,
                                            const arma::vec& rate,
                                            const int n,
                                            const int absorber)
{
  DirectedSparseLUFactor* fac = new DirectedSparseLUFactor(n, absorber);

  std::vector<double> row_sum(n, 0.0);
  for (arma::uword k = 0; k < edges.n_rows; ++k) {
    const int a = static_cast<int>(edges(k, 0)) - 1;
    row_sum[a] += rate[k];
  }

  std::vector<Triplet> triplets;
  triplets.reserve(static_cast<std::size_t>(edges.n_rows + n));
  for (int v = 0; v < n; ++v) {
    if (v == absorber) continue;
    triplets.push_back(Triplet(fac->rid[v], fac->rid[v], -row_sum[v]));
  }

  for (arma::uword k = 0; k < edges.n_rows; ++k) {
    const int a = static_cast<int>(edges(k, 0)) - 1;
    const int b = static_cast<int>(edges(k, 1)) - 1;
    if (a == absorber || b == absorber) continue;
    triplets.push_back(Triplet(fac->rid[a], fac->rid[b], rate[k]));
  }

  SpMat Q(n - 1, n - 1);
  Q.setFromTriplets(triplets.begin(), triplets.end());
  Q.makeCompressed();

  fac->solver.compute(Q);
  if (fac->solver.info() != Eigen::Success) {
    delete fac;
    Rcpp::stop("Directed SparseLU factorization failed.");
  }

  Eigen::VectorXd rhs = Eigen::VectorXd::Constant(n - 1, -1.0);
  Eigen::VectorXd hred = fac->solver.solve(rhs);
  if (fac->solver.info() != Eigen::Success) {
    delete fac;
    Rcpp::stop("Directed SparseLU forward solve failed.");
  }

  for (int r = 0; r < n - 1; ++r) fac->hfull[fac->idx[r]] = hred[r];
  return fac;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List directed_sparse_lu_forward_cpp(const arma::umat& edges,
                                          const arma::vec& rate,
                                          const arma::uvec& focal,
                                          const int n)
{
  validate_directed_inputs(edges, rate, focal, n);
  const int nf = static_cast<int>(focal.n_elem);
  arma::mat Hf(nf, nf, arma::fill::zeros);
  Rcpp::List hcache(nf);

  for (int fj = 0; fj < nf; ++fj) {
    const int absorber = static_cast<int>(focal[fj]) - 1;
    DirectedSparseLUFactor* fac = factor_one_absorber(edges, rate, n, absorber);

    for (int fi = 0; fi < nf; ++fi) {
      const int v = static_cast<int>(focal[fi]) - 1;
      Hf(fi, fj) = fac->hfull[v];
    }

    Rcpp::XPtr<DirectedSparseLUFactor> ptr(fac, true);
    hcache[fj] = ptr;
  }

  return Rcpp::List::create(
    Rcpp::Named("Hf") = Hf,
    Rcpp::Named("hcache") = hcache
  );
}

// [[Rcpp::export]]
arma::vec directed_sparse_lu_adjoint_cpp(Rcpp::List hcache,
                                         const arma::umat& edges,
                                         const arma::uvec& focal,
                                         const arma::mat& dL_dH,
                                         const int n)
{
  const int nf = static_cast<int>(focal.n_elem);
  if (hcache.size() != nf) Rcpp::stop("`hcache` length must match `focal`.");
  if (dL_dH.n_rows != focal.n_elem || dL_dH.n_cols != focal.n_elem)
    Rcpp::stop("`dL_dH` must be a square focal-by-focal matrix.");

  arma::vec dL_drate(edges.n_rows, arma::fill::zeros);

  for (int fj = 0; fj < nf; ++fj) {
    Rcpp::XPtr<DirectedSparseLUFactor> fac(hcache[fj]);
    if (fac->n != n) Rcpp::stop("SparseLU factor cache has incompatible `n`.");

    Eigen::VectorXd bred = Eigen::VectorXd::Zero(n - 1);
    for (int fi = 0; fi < nf; ++fi) {
      const int v = static_cast<int>(focal[fi]) - 1;
      if (v == fac->absorber) continue;
      bred[fac->rid[v]] = dL_dH(fi, fj);
    }

    Eigen::VectorXd ared = fac->solver.transpose().solve(bred);
    if (fac->solver.info() != Eigen::Success)
      Rcpp::stop("Directed SparseLU transpose solve failed.");

    Eigen::VectorXd afull = Eigen::VectorXd::Zero(n);
    for (int r = 0; r < n - 1; ++r) afull[fac->idx[r]] = ared[r];

    for (arma::uword k = 0; k < edges.n_rows; ++k) {
      const int a = static_cast<int>(edges(k, 0)) - 1;
      const int b = static_cast<int>(edges(k, 1)) - 1;
      if (a == fac->absorber) continue;
      dL_drate[k] += afull[a] * (fac->hfull[a] - fac->hfull[b]);
    }
  }

  return dL_drate;
}
