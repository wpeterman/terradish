#include <RcppArmadillo.h>
#include <RcppEigen.h>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/cg.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

namespace {

using AmgBackend = amgcl::backend::builtin<double>;
using AmgPreconditioner = amgcl::amg<
  AmgBackend,
  amgcl::coarsening::smoothed_aggregation,
  amgcl::relaxation::spai0
>;
using AmgIterativeSolver = amgcl::solver::cg<AmgBackend>;
using AmgSolver = amgcl::make_solver<AmgPreconditioner, AmgIterativeSolver>;
using RowSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

void validate_edge_pairs(const arma::vec& conductance, const arma::umat& edge_pairs)
{
  if (edge_pairs.n_cols != 2)
    Rcpp::stop("`edge_pairs` must have two columns");

  const arma::uword n_vertices = conductance.n_elem;
  if (n_vertices < 2)
    Rcpp::stop("need at least two vertices");

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword i = edge_pairs(edge, 0);
    const arma::uword j = edge_pairs(edge, 1);
    if (i < 1 || j < 1 || i > n_vertices || j > n_vertices || i == j)
      Rcpp::stop("vertex index out of bounds");
  }
}

RowSparseMatrix reduced_laplacian_row_major(
    const arma::vec& conductance,
    const arma::umat& edge_pairs)
{
  const arma::uword n_vertices = conductance.n_elem;
  const arma::uword n_reduced = n_vertices - 1;

  std::vector<Eigen::Triplet<double> > triplets;
  triplets.reserve(2 * edge_pairs.n_rows + n_reduced);
  std::vector<double> diagonal(static_cast<std::size_t>(n_reduced), 0.0);

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword ii = edge_pairs(edge, 0) - 1;
    const arma::uword jj = edge_pairs(edge, 1) - 1;
    const double w = conductance.at(ii) + conductance.at(jj);

    if (ii < n_reduced)
      diagonal[static_cast<std::size_t>(ii)] += w;
    if (jj < n_reduced)
      diagonal[static_cast<std::size_t>(jj)] += w;

    if (ii < n_reduced && jj < n_reduced)
    {
      triplets.emplace_back(static_cast<int>(ii), static_cast<int>(jj), -w);
      triplets.emplace_back(static_cast<int>(jj), static_cast<int>(ii), -w);
    }
  }

  for (arma::uword i = 0; i < n_reduced; ++i)
    triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), diagonal[static_cast<std::size_t>(i)]);

  RowSparseMatrix matrix(static_cast<int>(n_reduced), static_cast<int>(n_reduced));
  matrix.setFromTriplets(triplets.begin(), triplets.end());
  matrix.makeCompressed();
  return matrix;
}

void export_crs_arrays(const RowSparseMatrix& matrix,
                       std::vector<ptrdiff_t>& ptr,
                       std::vector<ptrdiff_t>& col,
                       std::vector<double>& val)
{
  const std::size_t n_reduced = static_cast<std::size_t>(matrix.rows());
  ptr.resize(n_reduced + 1);
  for (std::size_t i = 0; i < n_reduced + 1; ++i)
    ptr[i] = static_cast<ptrdiff_t>(matrix.outerIndexPtr()[static_cast<Eigen::Index>(i)]);

  const std::size_t nnz = static_cast<std::size_t>(matrix.nonZeros());
  col.resize(nnz);
  val.resize(nnz);
  for (std::size_t k = 0; k < nnz; ++k)
  {
    col[k] = static_cast<ptrdiff_t>(matrix.innerIndexPtr()[static_cast<Eigen::Index>(k)]);
    val[k] = matrix.valuePtr()[static_cast<Eigen::Index>(k)];
  }
}

class ReducedLaplacianAmgSolver {
public:
  ReducedLaplacianAmgSolver(
      const arma::vec& conductance,
      const arma::umat& edge_pairs,
      const double tol,
      const int maxit,
      const int coarse_enough,
      const int npre,
      const int npost,
      const double sa_relax,
      const double aggr_eps_strong,
      const bool estimate_spectral_radius,
      const int power_iters)
    : n_reduced_(0),
      nnz_(0),
      setup_time_(0.0),
      solver_()
  {
    validate_edge_pairs(conductance, edge_pairs);

    const auto setup_start = std::chrono::steady_clock::now();
    RowSparseMatrix matrix = reduced_laplacian_row_major(conductance, edge_pairs);
    n_reduced_ = static_cast<std::size_t>(matrix.rows());
    nnz_ = static_cast<std::size_t>(matrix.nonZeros());
    export_crs_arrays(matrix, ptr_, col_, val_);
    const auto matrix_view = std::tie(n_reduced_, ptr_, col_, val_);

    solver_.reset(new AmgSolver(
      matrix_view,
      build_params(
        tol,
        maxit,
        coarse_enough,
        npre,
        npost,
        sa_relax,
        aggr_eps_strong,
        estimate_spectral_radius,
        power_iters
      )
    ));

    const auto setup_end = std::chrono::steady_clock::now();
    setup_time_ = std::chrono::duration_cast<std::chrono::duration<double> >(setup_end - setup_start).count();
  }

  Rcpp::List rebuild(const arma::vec& conductance, const arma::umat& edge_pairs)
  {
    validate_edge_pairs(conductance, edge_pairs);

    const auto setup_start = std::chrono::steady_clock::now();
    RowSparseMatrix matrix = reduced_laplacian_row_major(conductance, edge_pairs);

    if (static_cast<std::size_t>(matrix.rows()) != n_reduced_)
      Rcpp::stop("[amg_reduced_laplacian_rebuild] matrix dimension changed");

    nnz_ = static_cast<std::size_t>(matrix.nonZeros());
    export_crs_arrays(matrix, ptr_, col_, val_);
    const auto matrix_view = std::tie(n_reduced_, ptr_, col_, val_);
    solver_->precond().rebuild(matrix_view);

    const auto setup_end = std::chrono::steady_clock::now();
    setup_time_ = std::chrono::duration_cast<std::chrono::duration<double> >(setup_end - setup_start).count();

    return Rcpp::List::create(
      Rcpp::Named("setup_time") = setup_time_,
      Rcpp::Named("n_reduced") = static_cast<double>(n_reduced_),
      Rcpp::Named("nnz") = static_cast<double>(nnz_)
    );
  }

  Rcpp::List solve(const arma::mat& rhs,
                   Rcpp::Nullable<arma::mat> x0,
                   Rcpp::Nullable<double> tol = R_NilValue,
                   Rcpp::Nullable<int> maxit = R_NilValue)
  {
    if (rhs.n_rows != n_reduced_)
      Rcpp::stop("[amg_reduced_laplacian_solve] rhs has incompatible number of rows");

    if (tol.isNotNull())
      solver_->solver().prm.tol = std::max(0.0, Rcpp::as<double>(tol));
    if (maxit.isNotNull())
      solver_->solver().prm.maxiter = static_cast<std::size_t>(std::max(1, Rcpp::as<int>(maxit)));

    arma::mat initial_guess;
    if (x0.isNotNull())
    {
      initial_guess = Rcpp::as<arma::mat>(x0);
      if (initial_guess.n_rows != rhs.n_rows || initial_guess.n_cols != rhs.n_cols)
        Rcpp::stop("[amg_reduced_laplacian_solve] x0 must match rhs dimensions");
    }

    arma::mat solution(rhs.n_rows, rhs.n_cols, arma::fill::zeros);
    Rcpp::LogicalVector converged(rhs.n_cols, false);
    Rcpp::IntegerVector iterations(rhs.n_cols, 0);
    arma::vec residual_norm(rhs.n_cols, arma::fill::zeros);

    const auto solve_start = std::chrono::steady_clock::now();
    for (arma::uword col = 0; col < rhs.n_cols; ++col)
    {
      std::vector<double> b(rhs.colptr(col), rhs.colptr(col) + rhs.n_rows);
      std::vector<double> x(rhs.n_rows, 0.0);

      if (x0.isNotNull())
        std::copy(initial_guess.colptr(col), initial_guess.colptr(col) + rhs.n_rows, x.begin());

      std::size_t iter = 0;
      double resid = 0.0;
      std::tie(iter, resid) = (*solver_)(b, x);

      iterations[col] = static_cast<int>(iter);
      residual_norm.at(col) = resid;
      converged[col] = std::isfinite(resid) && resid <= solver_->solver().prm.tol;
      std::copy(x.begin(), x.end(), solution.colptr(col));
    }
    const auto solve_end = std::chrono::steady_clock::now();
    const double solve_time =
      std::chrono::duration_cast<std::chrono::duration<double> >(solve_end - solve_start).count();

    return Rcpp::List::create(
      Rcpp::Named("solution") = solution,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("residual_norm") = residual_norm,
      Rcpp::Named("setup_time") = setup_time_,
      Rcpp::Named("solve_time") = solve_time,
      Rcpp::Named("n_reduced") = static_cast<double>(n_reduced_),
      Rcpp::Named("nnz") = static_cast<double>(nnz_)
    );
  }

private:
  static AmgSolver::params build_params(const double tol,
                                        const int maxit,
                                        const int coarse_enough,
                                        const int npre,
                                        const int npost,
                                        const double sa_relax,
                                        const double aggr_eps_strong,
                                        const bool estimate_spectral_radius,
                                        const int power_iters)
  {
    AmgSolver::params prm;
    prm.solver.tol = tol;
    prm.solver.maxiter = static_cast<std::size_t>(std::max(1, maxit));
    prm.precond.coarse_enough = static_cast<unsigned>(std::max(1, coarse_enough));
    prm.precond.npre = static_cast<unsigned>(std::max(0, npre));
    prm.precond.npost = static_cast<unsigned>(std::max(0, npost));
    prm.precond.coarsening.relax = static_cast<float>(sa_relax);
    prm.precond.coarsening.aggr.eps_strong = static_cast<float>(aggr_eps_strong);
    prm.precond.coarsening.estimate_spectral_radius = estimate_spectral_radius;
    prm.precond.coarsening.power_iters = std::max(0, power_iters);
    prm.precond.allow_rebuild = true;
    return prm;
  }

  std::size_t n_reduced_;
  std::size_t nnz_;
  std::vector<ptrdiff_t> ptr_;
  std::vector<ptrdiff_t> col_;
  std::vector<double> val_;
  double setup_time_;
  std::unique_ptr<AmgSolver> solver_;
};

using AmgSolverPtr = Rcpp::XPtr<ReducedLaplacianAmgSolver>;

} // namespace

// [[Rcpp::export]]
SEXP amg_reduced_laplacian_create(const arma::vec& conductance,
                                  const arma::umat& edge_pairs,
                                  const double tol = 1e-8,
                                  const int maxit = 1000,
                                  const int coarse_enough = 500,
                                  const int npre = 1,
                                  const int npost = 1,
                                  const double sa_relax = 1.0,
                                  const double aggr_eps_strong = 0.08,
                                  const bool estimate_spectral_radius = true,
                                  const int power_iters = 4)
{
  AmgSolverPtr ptr(new ReducedLaplacianAmgSolver(
      conductance,
      edge_pairs,
      tol,
      maxit,
      coarse_enough,
      npre,
      npost,
      sa_relax,
      aggr_eps_strong,
      estimate_spectral_radius,
      power_iters
  ), true);
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List amg_reduced_laplacian_solve(SEXP solver_ptr,
                                       const arma::mat& rhs,
                                       Rcpp::Nullable<arma::mat> x0 = R_NilValue,
                                       Rcpp::Nullable<double> tol = R_NilValue,
                                       Rcpp::Nullable<int> maxit = R_NilValue)
{
  AmgSolverPtr ptr(solver_ptr);
  return ptr->solve(rhs, x0, tol, maxit);
}

// [[Rcpp::export]]
Rcpp::List amg_reduced_laplacian_rebuild(SEXP solver_ptr,
                                         const arma::vec& conductance,
                                         const arma::umat& edge_pairs)
{
  AmgSolverPtr ptr(solver_ptr);
  return ptr->rebuild(conductance, edge_pairs);
}
