#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

namespace {

void validate_edge_pairs(const arma::vec& conductance, const arma::umat& edge_pairs)
{
  if (edge_pairs.n_cols != 2)
    Rcpp::stop("`edge_pairs` must have two columns");

  const arma::uword N = conductance.n_elem;
  if (N < 2)
    Rcpp::stop("need at least two vertices");

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword i = edge_pairs(edge, 0);
    const arma::uword j = edge_pairs(edge, 1);
    if (i < 1 || j < 1 || i > N || j > N || i == j)
      Rcpp::stop("vertex index out of bounds");
  }
}

arma::vec reduced_laplacian_diagonal(const arma::vec& conductance, const arma::umat& edge_pairs)
{
  const arma::uword N = conductance.n_elem;
  const arma::uword Nred = N - 1;
  arma::vec diagonal(Nred, arma::fill::zeros);

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword ii = edge_pairs(edge, 0) - 1;
    const arma::uword jj = edge_pairs(edge, 1) - 1;
    const double w = conductance.at(ii) + conductance.at(jj);

    if (ii < Nred)
      diagonal.at(ii) += w;
    if (jj < Nred)
      diagonal.at(jj) += w;
  }

  return diagonal;
}

arma::vec reduced_laplacian_matvec(const arma::vec& conductance,
                                   const arma::umat& edge_pairs,
                                   const arma::vec& x)
{
  const arma::uword N = conductance.n_elem;
  const arma::uword Nred = N - 1;
  if (x.n_elem != Nred)
    Rcpp::stop("dimension mismatch in reduced_laplacian_matvec");

  arma::vec y(Nred, arma::fill::zeros);
  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword ii = edge_pairs(edge, 0) - 1;
    const arma::uword jj = edge_pairs(edge, 1) - 1;
    const double w = conductance.at(ii) + conductance.at(jj);

    if (ii < Nred)
      y.at(ii) += w * x.at(ii);
    if (jj < Nred)
      y.at(jj) += w * x.at(jj);

    if (ii < Nred && jj < Nred)
    {
      y.at(ii) -= w * x.at(jj);
      y.at(jj) -= w * x.at(ii);
    }
  }

  return y;
}

Eigen::SparseMatrix<double> reduced_laplacian_eigen(const arma::vec& conductance,
                                                    const arma::umat& edge_pairs)
{
  const arma::uword N = conductance.n_elem;
  const arma::uword Nred = N - 1;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(2 * edge_pairs.n_rows + Nred);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Nred));

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword ii = edge_pairs(edge, 0) - 1;
    const arma::uword jj = edge_pairs(edge, 1) - 1;
    const double w = conductance.at(ii) + conductance.at(jj);

    if (ii < Nred)
      diagonal[static_cast<Eigen::Index>(ii)] += w;
    if (jj < Nred)
      diagonal[static_cast<Eigen::Index>(jj)] += w;

    if (ii < Nred && jj < Nred)
    {
      triplets.emplace_back(static_cast<Eigen::Index>(ii), static_cast<Eigen::Index>(jj), -w);
      triplets.emplace_back(static_cast<Eigen::Index>(jj), static_cast<Eigen::Index>(ii), -w);
    }
  }

  for (arma::uword i = 0; i < Nred; ++i)
    triplets.emplace_back(static_cast<Eigen::Index>(i),
                          static_cast<Eigen::Index>(i),
                          diagonal[static_cast<Eigen::Index>(i)]);

  Eigen::SparseMatrix<double> A(static_cast<Eigen::Index>(Nred),
                                static_cast<Eigen::Index>(Nred));
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}

} // namespace

// [[Rcpp::export]]
arma::sp_mat assemble_reduced_laplacian(const arma::vec& conductance, const arma::umat& edge_pairs)
{
  validate_edge_pairs(conductance, edge_pairs);

  const arma::uword N = conductance.n_elem;
  const arma::uword Nred = N - 1;
  arma::vec diagonal = reduced_laplacian_diagonal(conductance, edge_pairs);

  std::vector<arma::uword> rows;
  std::vector<arma::uword> cols;
  std::vector<double> vals;
  rows.reserve(2 * edge_pairs.n_rows + Nred);
  cols.reserve(2 * edge_pairs.n_rows + Nred);
  vals.reserve(2 * edge_pairs.n_rows + Nred);

  for (arma::uword edge = 0; edge < edge_pairs.n_rows; ++edge)
  {
    const arma::uword i = edge_pairs(edge, 0);
    const arma::uword j = edge_pairs(edge, 1);

    const arma::uword ii = i - 1;
    const arma::uword jj = j - 1;
    const double w = conductance.at(ii) + conductance.at(jj);

    if (ii < Nred && jj < Nred)
    {
      rows.push_back(ii);
      cols.push_back(jj);
      vals.push_back(-w);
      rows.push_back(jj);
      cols.push_back(ii);
      vals.push_back(-w);
    }
  }

  for (arma::uword i = 0; i < Nred; ++i)
  {
    rows.push_back(i);
    cols.push_back(i);
    vals.push_back(diagonal.at(i));
  }

  arma::umat locations(2, rows.size());
  arma::vec values(vals.size());
  for (arma::uword k = 0; k < rows.size(); ++k)
  {
    locations(0, k) = rows[k];
    locations(1, k) = cols[k];
    values.at(k) = vals[k];
  }

  return arma::sp_mat(locations, values, Nred, Nred, true, true);
}

// [[Rcpp::export]]
Rcpp::List pcg_reduced_laplacian_ic(const arma::mat& rhs,
                                    const arma::vec& conductance,
                                    const arma::umat& edge_pairs,
                                    Rcpp::Nullable<arma::mat> x0 = R_NilValue,
                                    const double tol = 1e-8,
                                    const int maxit = 1000)
{
  validate_edge_pairs(conductance, edge_pairs);

  const arma::uword Nred = conductance.n_elem - 1;
  if (rhs.n_rows != Nred)
    Rcpp::stop("[pcg_reduced_laplacian_ic] rhs has incompatible number of rows");

  arma::mat initial_guess;
  if (x0.isNotNull())
  {
    initial_guess = Rcpp::as<arma::mat>(x0);
    if (initial_guess.n_rows != rhs.n_rows || initial_guess.n_cols != rhs.n_cols)
      Rcpp::stop("[pcg_reduced_laplacian_ic] x0 must match rhs dimensions");
  }

  Eigen::SparseMatrix<double> A = reduced_laplacian_eigen(conductance, edge_pairs);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>> solver;
  solver.setTolerance(tol);
  solver.setMaxIterations(maxit);
  solver.compute(A);
  if (solver.info() != Eigen::Success)
    Rcpp::stop("[pcg_reduced_laplacian_ic] incomplete Cholesky factorization failed");

  arma::mat solution(rhs.n_rows, rhs.n_cols, arma::fill::zeros);
  Rcpp::LogicalVector converged(rhs.n_cols, false);
  Rcpp::IntegerVector iterations(rhs.n_cols, maxit);
  arma::vec residual_norm(rhs.n_cols, arma::fill::zeros);

  for (arma::uword col = 0; col < rhs.n_cols; ++col)
  {
    Eigen::Map<const Eigen::VectorXd> b(rhs.colptr(col), static_cast<Eigen::Index>(rhs.n_rows));
    Eigen::VectorXd x;
    if (x0.isNotNull())
    {
      Eigen::Map<const Eigen::VectorXd> guess(initial_guess.colptr(col),
                                              static_cast<Eigen::Index>(rhs.n_rows));
      x = solver.solveWithGuess(b, guess);
    }
    else
      x = solver.solve(b);

    const double bnorm = std::max(1.0, b.norm());
    converged[col] = solver.info() == Eigen::Success;
    iterations[col] = solver.iterations();
    residual_norm.at(col) = solver.error() * bnorm;
    std::copy(x.data(), x.data() + x.size(), solution.colptr(col));
  }

  return Rcpp::List::create(
    Rcpp::Named("solution") = solution,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("residual_norm") = residual_norm
  );
}

// [[Rcpp::export]]
Rcpp::List pcg_reduced_laplacian(const arma::mat& rhs,
                                 const arma::vec& conductance,
                                 const arma::umat& edge_pairs,
                                 Rcpp::Nullable<arma::mat> x0 = R_NilValue,
                                 const double tol = 1e-8,
                                 const int maxit = 1000)
{
  validate_edge_pairs(conductance, edge_pairs);

  const arma::uword Nred = conductance.n_elem - 1;
  if (rhs.n_rows != Nred)
    Rcpp::stop("[pcg_reduced_laplacian] rhs has incompatible number of rows");

  arma::mat initial_guess;
  if (x0.isNotNull())
  {
    initial_guess = Rcpp::as<arma::mat>(x0);
    if (initial_guess.n_rows != rhs.n_rows || initial_guess.n_cols != rhs.n_cols)
      Rcpp::stop("[pcg_reduced_laplacian] x0 must match rhs dimensions");
  }

  arma::vec diagonal = reduced_laplacian_diagonal(conductance, edge_pairs);
  if (arma::any(diagonal <= 0))
    Rcpp::stop("[pcg_reduced_laplacian] non-positive diagonal encountered");

  arma::mat solution(rhs.n_rows, rhs.n_cols, arma::fill::zeros);
  Rcpp::LogicalVector converged(rhs.n_cols, false);
  Rcpp::IntegerVector iterations(rhs.n_cols, maxit);
  arma::vec residual_norm(rhs.n_cols, arma::fill::zeros);

  for (arma::uword col = 0; col < rhs.n_cols; ++col)
  {
    const arma::vec b = rhs.col(col);
    arma::vec x = x0.isNotNull() ? initial_guess.col(col) : arma::vec(rhs.n_rows, arma::fill::zeros);
    arma::vec r = b - reduced_laplacian_matvec(conductance, edge_pairs, x);
    const double bnorm = std::max(arma::norm(b, 2), 1.0);
    double rnorm = arma::norm(r, 2);

    if (rnorm <= tol * bnorm)
    {
      converged[col] = true;
      iterations[col] = 0;
      residual_norm.at(col) = rnorm;
      solution.col(col) = x;
      continue;
    }

    arma::vec z = r / diagonal;
    arma::vec p = z;
    double rz_old = arma::dot(r, z);

    for (int iter = 0; iter < maxit; ++iter)
    {
      arma::vec Ap = reduced_laplacian_matvec(conductance, edge_pairs, p);
      const double denom = arma::dot(p, Ap);
      if (denom <= 0)
        Rcpp::stop("[pcg_reduced_laplacian] non-positive search direction");

      const double alpha = rz_old / denom;
      x += alpha * p;
      r -= alpha * Ap;
      rnorm = arma::norm(r, 2);

      if (rnorm <= tol * bnorm)
      {
        converged[col] = true;
        iterations[col] = iter + 1;
        residual_norm.at(col) = rnorm;
        solution.col(col) = x;
        break;
      }

      z = r / diagonal;
      const double rz_new = arma::dot(r, z);
      const double beta = rz_new / rz_old;
      p = z + beta * p;
      rz_old = rz_new;

      if (iter == maxit - 1)
      {
        residual_norm.at(col) = rnorm;
        solution.col(col) = x;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("solution") = solution,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("residual_norm") = residual_norm
  );
}

// [[Rcpp::export]]
arma::vec backpropagate_laplacian_to_conductance (const arma::mat& tGl, const arma::mat& tGr, const arma::umat& tadj)
{
  if (tadj.n_rows != 2)
    Rcpp::stop("[backpropagate_gradient] tadj.n_rows != 2");

  if (tGl.n_rows != tGr.n_rows || tGl.n_cols != tGr.n_cols)
    Rcpp::stop("[backpropagate_gradient] dim(tGl) != dim(tGr)");

  const unsigned N = tGl.n_cols;

  arma::vec djj (N);
  for (unsigned vertex = 0; vertex < N; ++vertex)
    djj.at(vertex) = -arma::dot(tGl.col(vertex), tGr.col(vertex));

  arma::vec dl_dC (N+1, arma::fill::zeros);
  unsigned i, j;
  double dij;

  for (unsigned edge = 0; edge < tadj.n_cols; ++edge)
  {
    j = tadj.at(0, edge);
    i = tadj.at(1, edge);
    if (i <= j) 
      continue;
    else if (i > N)
      Rcpp::stop("[backpropagate_gradient] vertex out of bounds");
    else if (i == N)
    {
      dl_dC.at(i) += djj.at(j);
      dl_dC.at(j) += djj.at(j);
    }
    else
    {
      dij = djj.at(i) + djj.at(j) + 
        arma::dot(tGl.col(i), tGr.col(j)) +
        arma::dot(tGl.col(j), tGr.col(i));
      dl_dC.at(i) += dij;
      dl_dC.at(j) += dij;
    }
  }

  return dl_dC;
}

// [[Rcpp::export]]
arma::sp_mat backpropagate_conductance_to_laplacian (const arma::vec& dgrad__ddl_dC, const arma::umat& tadj)
{
  if (tadj.n_rows != 2)
    Rcpp::stop("[backpropagate_conductance] tadj.n_rows != 2");

  const unsigned N = dgrad__ddl_dC.n_elem - 1;

  unsigned i, j;
  double dij;

  arma::vec diagonal (N, arma::fill::zeros),
            offdiagonal (tadj.n_cols, arma::fill::zeros);

  for (unsigned edge = 0; edge < tadj.n_cols; ++edge)
  {
    j = tadj.at(0, edge);
    i = tadj.at(1, edge);

    if (i <= j)
      continue;
    else if (i > N)
      Rcpp::stop("[backpropagate_conductance] vertex out of bounds");

    dij = dgrad__ddl_dC.at(i) + dgrad__ddl_dC.at(j);
    diagonal.at(j) += dij;
    if (i != N)
    {
      offdiagonal.at(edge) = -dij;
      diagonal.at(i) += dij;
    }
  }

  arma::sp_mat dgrad__ddl_dQn (tadj, offdiagonal, N, N);
  dgrad__ddl_dQn.diag() = diagonal;

  return arma::trimatu(dgrad__ddl_dQn);
}

// [[Rcpp::export]]
arma::mat laplacian_derivative_matrix_product(const arma::vec& dgrad__ddl_dC,
                                              const arma::umat& tadj,
                                              const arma::mat& G)
{
  if (tadj.n_rows != 2)
    Rcpp::stop("[laplacian_derivative_matrix_product] tadj.n_rows != 2");

  const unsigned N = dgrad__ddl_dC.n_elem - 1;
  if (G.n_rows != N)
    Rcpp::stop("[laplacian_derivative_matrix_product] G has incompatible number of rows");

  arma::mat dgrad__ddl_dQnG(N, G.n_cols, arma::fill::zeros);

  for (unsigned edge = 0; edge < tadj.n_cols; ++edge)
  {
    const unsigned j = tadj.at(0, edge);
    const unsigned i = tadj.at(1, edge);

    if (i <= j)
      continue;
    else if (i > N)
      Rcpp::stop("[laplacian_derivative_matrix_product] vertex out of bounds");

    const double dij = dgrad__ddl_dC.at(i) + dgrad__ddl_dC.at(j);
    if (i == N)
    {
      dgrad__ddl_dQnG.row(j) += dij * G.row(j);
    }
    else
    {
      dgrad__ddl_dQnG.row(j) += dij * (G.row(j) - G.row(i));
      dgrad__ddl_dQnG.row(i) += dij * (G.row(i) - G.row(j));
    }
  }

  return dgrad__ddl_dQnG;
}

// [[Rcpp::export]]
arma::mat graph_rhs_matrix_product(const arma::uvec& demes,
                                   const int n_vertices,
                                   const arma::mat& X)
{
  if (n_vertices < 2)
    Rcpp::stop("[graph_rhs_matrix_product] need at least two vertices");
  if (X.n_rows != demes.n_elem)
    Rcpp::stop("[graph_rhs_matrix_product] X has incompatible number of rows");

  const arma::uword N = static_cast<arma::uword>(n_vertices);
  const arma::uword Nred = N - 1;

  arma::rowvec baseline = -arma::sum(X, 0) / static_cast<double>(N);
  arma::mat out(Nred, X.n_cols, arma::fill::zeros);
  out.each_row() = baseline;

  for (arma::uword j = 0; j < demes.n_elem; ++j)
  {
    const arma::uword deme = demes.at(j);
    if (deme < 1 || deme > N)
      Rcpp::stop("[graph_rhs_matrix_product] deme index out of bounds");
    if (deme < N)
      out.row(deme - 1) += X.row(j);
  }

  return out;
}

// [[Rcpp::export]]
arma::mat graph_rhs_crossprod(const arma::uvec& demes,
                              const int n_vertices,
                              const arma::mat& X)
{
  if (n_vertices < 2)
    Rcpp::stop("[graph_rhs_crossprod] need at least two vertices");

  const arma::uword N = static_cast<arma::uword>(n_vertices);
  const arma::uword Nred = N - 1;
  if (X.n_rows != Nred)
    Rcpp::stop("[graph_rhs_crossprod] X has incompatible number of rows");

  arma::rowvec baseline = -arma::sum(X, 0) / static_cast<double>(N);
  arma::mat out(demes.n_elem, X.n_cols, arma::fill::zeros);
  out.each_row() = baseline;

  for (arma::uword j = 0; j < demes.n_elem; ++j)
  {
    const arma::uword deme = demes.at(j);
    if (deme < 1 || deme > N)
      Rcpp::stop("[graph_rhs_crossprod] deme index out of bounds");
    if (deme < N)
      out.row(j) += X.row(deme - 1);
  }

  return out;
}

namespace {

arma::vec block_cg_diagonal(const arma::vec& conductance, const arma::umat& edge_pairs)
{
  const arma::uword Nred = conductance.n_elem - 1;
  arma::vec d(Nred, arma::fill::zeros);
  for (arma::uword e = 0; e < edge_pairs.n_rows; ++e)
  {
    const arma::uword i = edge_pairs(e, 0) - 1, j = edge_pairs(e, 1) - 1;
    const double w = conductance.at(i) + conductance.at(j);
    if (i < Nred) d.at(i) += w;
    if (j < Nred) d.at(j) += w;
  }
  return d;
}

arma::mat block_cg_matmat(const arma::vec& conductance, const arma::umat& edge_pairs,
                          const arma::mat& X)
{
  const arma::uword Nred = conductance.n_elem - 1;
  arma::mat Y(Nred, X.n_cols, arma::fill::zeros);
  for (arma::uword e = 0; e < edge_pairs.n_rows; ++e)
  {
    const arma::uword i = edge_pairs(e, 0) - 1, j = edge_pairs(e, 1) - 1;
    const double w = conductance.at(i) + conductance.at(j);
    if (i < Nred) Y.row(i) += w * X.row(i);
    if (j < Nred) Y.row(j) += w * X.row(j);
    if (i < Nred && j < Nred) { Y.row(i) -= w * X.row(j); Y.row(j) -= w * X.row(i); }
  }
  return Y;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List block_cg_reduced_laplacian(const arma::mat& rhs,
                                      const arma::vec& conductance,
                                      const arma::umat& edge_pairs,
                                      Rcpp::Nullable<arma::mat> x0 = R_NilValue,
                                      const double tol = 1e-8,
                                      const int maxit = 1000)
{
  validate_edge_pairs(conductance, edge_pairs);

  const arma::uword Nred = conductance.n_elem - 1, s = rhs.n_cols;
  if (rhs.n_rows != Nred)
    Rcpp::stop("[block_cg_reduced_laplacian] rhs has incompatible number of rows");

  arma::vec d = block_cg_diagonal(conductance, edge_pairs);
  if (arma::any(d <= 0))
    Rcpp::stop("[block_cg_reduced_laplacian] non-positive diagonal encountered");

  arma::mat X = x0.isNotNull() ? Rcpp::as<arma::mat>(x0)
                               : arma::mat(Nred, s, arma::fill::zeros);
  if (X.n_rows != Nred || X.n_cols != s)
    Rcpp::stop("[block_cg_reduced_laplacian] x0 must match rhs dimensions");

  arma::mat R = rhs - block_cg_matmat(conductance, edge_pairs, X);
  arma::mat Z = R.each_col() / d;
  arma::mat P = Z;
  arma::mat RtZ = R.t() * Z;
  arma::vec bnorm = arma::sqrt(arma::sum(arma::square(rhs), 0)).t();
  bnorm.transform([](double v){ return std::max(v, 1.0); });

  int iter = 0;
  arma::vec rnorm = arma::sqrt(arma::sum(arma::square(R), 0)).t();
  for (iter = 0; iter < maxit; ++iter)
  {
    arma::mat Q = block_cg_matmat(conductance, edge_pairs, P);
    // pseudo-inverse for the small s x s systems: as RHS columns converge the
    // residual block loses rank (classic block-CG breakdown); pinv is stable.
    arma::mat alpha = arma::pinv(arma::mat(P.t() * Q)) * RtZ;
    X += P * alpha;
    R -= Q * alpha;
    rnorm = arma::sqrt(arma::sum(arma::square(R), 0)).t();
    if (arma::all(rnorm <= tol * bnorm)) { ++iter; break; }
    arma::mat Znew = R.each_col() / d;
    arma::mat RtZnew = R.t() * Znew;
    arma::mat beta = arma::pinv(RtZ) * RtZnew;
    P = Znew + P * beta;
    Z = Znew;
    RtZ = RtZnew;
  }

  Rcpp::LogicalVector converged(s);
  for (arma::uword k = 0; k < s; ++k)
    converged[k] = rnorm(k) <= tol * bnorm(k);

  return Rcpp::List::create(
    Rcpp::Named("solution") = X,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("residual_norm") = rnorm);
}
