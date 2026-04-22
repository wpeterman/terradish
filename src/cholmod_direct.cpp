#include <RcppArmadillo.h>
// Matrix has shipped its CHOLMOD registration wrappers under two include
// layouts. Support either path and fail with a clear message otherwise.
#if defined(__has_include)
#  if __has_include(<Matrix_stubs.c>)
#    include <Matrix_stubs.c>
#  elif __has_include(<Matrix/stubs.c>)
#    include <Matrix/stubs.c>
#  else
#    error "Compatible Matrix CHOLMOD interface headers not found. Install Matrix >= 1.6-2."
#  endif
#else
#  include <Matrix_stubs.c>
#endif
#include <chrono>
#include <limits>
#include <string>

// [[Rcpp::depends(RcppArmadillo, Matrix)]]

namespace {

class CholmodCommon
{
public:
  CholmodCommon() : started_(false)
  {
    if (!M_cholmod_start(&common_))
      Rcpp::stop("failed to initialize CHOLMOD");
    started_ = true;
  }

  ~CholmodCommon()
  {
    if (started_)
      M_cholmod_finish(&common_);
  }

  CHM_CM get()
  {
    return &common_;
  }

private:
  cholmod_common common_;
  bool started_;
};

void check_rhs_dimensions(const arma::mat& rhs)
{
  if (rhs.n_rows < 1 || rhs.n_cols < 1)
    Rcpp::stop("`rhs` must be a non-empty numeric matrix");
  if (rhs.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max()) ||
      rhs.n_cols > static_cast<arma::uword>(std::numeric_limits<int>::max()))
    Rcpp::stop("`rhs` dimensions exceed the Matrix CHOLMOD interface limits");
}

CHM_SP sparse_matrix_view(SEXP matrix, cholmod_sparse& sparse_view)
{
  if (!Rf_isS4(matrix))
    Rcpp::stop("`matrix` must be a Matrix sparse object");

  CHM_SP A = M_sexp_as_cholmod_sparse(&sparse_view, matrix, FALSE, FALSE);
  if (A == NULL)
    Rcpp::stop("failed to read Matrix sparse object");
  if (A->nrow != A->ncol)
    Rcpp::stop("CHOLMOD direct solver requires a square matrix");

  return A;
}

void configure_cholmod_common(CHM_CM common,
                              const std::string& factorization,
                              const bool perm)
{
  if (factorization == "simplicial_ldl")
  {
    common->supernodal = CHOLMOD_SIMPLICIAL;
    common->final_super = FALSE;
    common->final_ll = FALSE;
  }
  else if (factorization == "simplicial_ll")
  {
    common->supernodal = CHOLMOD_SIMPLICIAL;
    common->final_super = FALSE;
    common->final_ll = TRUE;
  }
  else if (factorization == "supernodal_ll")
  {
    common->supernodal = CHOLMOD_SUPERNODAL;
    common->final_super = TRUE;
    common->final_ll = TRUE;
  }
  else
    Rcpp::stop("Unknown direct factorization mode: " + factorization);

  if (!perm)
  {
    common->nmethods = 1;
    common->method[0].ordering = CHOLMOD_NATURAL;
    common->postorder = FALSE;
  }
}

arma::mat cholmod_solve_dense(CHM_FR factor,
                              const arma::mat& rhs,
                              CHM_CM common,
                              double& solve_time)
{
  check_rhs_dimensions(rhs);

  // CHOLMOD's dense wrapper takes mutable memory, so use a local copy even
  // though the solve itself does not conceptually alter the right-hand side.
  arma::mat rhs_copy(rhs);
  cholmod_dense rhs_view;
  CHM_DN B = M_numeric_as_cholmod_dense(
    &rhs_view,
    rhs_copy.memptr(),
    static_cast<int>(rhs_copy.n_rows),
    static_cast<int>(rhs_copy.n_cols)
  );
  if (B == NULL)
    Rcpp::stop("failed to wrap right-hand side for CHOLMOD");

  const auto start = std::chrono::steady_clock::now();
  CHM_DN X = M_cholmod_solve(CHOLMOD_A, factor, B, common);
  const auto finish = std::chrono::steady_clock::now();
  if (X == NULL)
    Rcpp::stop("CHOLMOD solve failed");

  arma::mat solution(
    static_cast<double*>(X->x),
    rhs_copy.n_rows,
    rhs_copy.n_cols,
    true
  );
  M_cholmod_free_dense(&X, common);

  solve_time = std::chrono::duration<double>(finish - start).count();
  return solution;
}

class CholmodDirectHandle
{
public:
  CholmodDirectHandle(SEXP matrix,
                      const std::string& factorization,
                      const bool perm) :
    common_(),
    factor_(NULL),
    n_(0)
  {
    configure_cholmod_common(common_.get(), factorization, perm);
    analyze_and_factorize(matrix);
  }

  ~CholmodDirectHandle()
  {
    if (factor_ != NULL)
      M_cholmod_free_factor(&factor_, common_.get());
  }

  Rcpp::List update(SEXP matrix)
  {
    cholmod_sparse sparse_view;
    CHM_SP A = sparse_matrix_view(matrix, sparse_view);
    if (A->nrow != n_)
      Rcpp::stop("matrix dimension changed; cannot reuse CHOLMOD factorization");

    const auto start = std::chrono::steady_clock::now();
    if (!M_cholmod_factorize(A, factor_, common_.get()))
      Rcpp::stop("CHOLMOD numeric refactorization failed");
    const auto finish = std::chrono::steady_clock::now();

    return Rcpp::List::create(
      Rcpp::Named("update_time") =
        std::chrono::duration<double>(finish - start).count()
    );
  }

  Rcpp::List solve(const arma::mat& rhs)
  {
    if (rhs.n_rows != static_cast<arma::uword>(n_))
      Rcpp::stop("`rhs` row count does not match the CHOLMOD factorization");

    double solve_time = 0.0;
    arma::mat solution = cholmod_solve_dense(
      factor_,
      rhs,
      common_.get(),
      solve_time
    );

    return Rcpp::List::create(
      Rcpp::Named("solution") = solution,
      Rcpp::Named("solve_time") = solve_time
    );
  }

private:
  void analyze_and_factorize(SEXP matrix)
  {
    cholmod_sparse sparse_view;
    CHM_SP A = sparse_matrix_view(matrix, sparse_view);
    n_ = A->nrow;

    factor_ = M_cholmod_analyze(A, common_.get());
    if (factor_ == NULL)
      Rcpp::stop("CHOLMOD symbolic analysis failed");

    if (!M_cholmod_factorize(A, factor_, common_.get()))
    {
      M_cholmod_free_factor(&factor_, common_.get());
      Rcpp::stop("CHOLMOD numeric factorization failed");
    }
  }

  CholmodCommon common_;
  CHM_FR factor_;
  size_t n_;
};

} // namespace

// [[Rcpp::export]]
Rcpp::List cholmod_factor_solve(SEXP factor, const arma::mat& rhs)
{
  if (!Rf_isS4(factor))
    Rcpp::stop("`factor` must be a Matrix CHMfactor object");

  CholmodCommon common;

  cholmod_factor factor_view;
  CHM_FR L = M_sexp_as_cholmod_factor(&factor_view, factor);
  if (L == NULL)
    Rcpp::stop("failed to read Matrix CHMfactor object");

  double solve_time = 0.0;
  arma::mat solution = cholmod_solve_dense(L, rhs, common.get(), solve_time);

  return Rcpp::List::create(
    Rcpp::Named("solution") = solution,
    Rcpp::Named("solve_time") = solve_time
  );
}

// [[Rcpp::export]]
SEXP cholmod_direct_create(SEXP matrix,
                           const std::string& factorization = "simplicial_ldl",
                           const bool perm = true)
{
  Rcpp::XPtr<CholmodDirectHandle> ptr(
    new CholmodDirectHandle(matrix, factorization, perm),
    true
  );
  return ptr;
}

// [[Rcpp::export]]
Rcpp::List cholmod_direct_update(SEXP solver_ptr, SEXP matrix)
{
  Rcpp::XPtr<CholmodDirectHandle> ptr(solver_ptr);
  return ptr->update(matrix);
}

// [[Rcpp::export]]
Rcpp::List cholmod_direct_solve(SEXP solver_ptr, const arma::mat& rhs)
{
  Rcpp::XPtr<CholmodDirectHandle> ptr(solver_ptr);
  return ptr->solve(rhs);
}
