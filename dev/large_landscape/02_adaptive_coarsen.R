# =============================================================================
# Prototype B: adaptive, current-density-driven multiresolution coarsening
# -----------------------------------------------------------------------------
# Tier-1 STRUCTURAL win from SCALING_DESIGN.md. Attacks N itself. The refinement
# indicator IS the information-density map: current density between focal pairs
# plus conductance-gradient magnitude. Resolution follows where the data are
# informative about theta; homogeneous low-current regions are coarsened.
#
# This prototype demonstrates the core loop on a regular grid:
#   1. fine solve -> per-edge current density J_e = w_e * |phi_i - phi_j|
#   2. block-aggregate ONLY homogeneous, low-current, non-focal blocks
#   3. keep focal cells + high-current + high-gradient cells at native resolution
#   4. re-solve coarse graph and report ||E_fine - E_coarse|| / ||E_fine||
#
# KEY CORRECTNESS RULES (see design doc, approach B risks):
#   * never coarsen across a one-cell corridor/barrier (protect high-gradient)
#   * conductance aggregation is LINEAR here (block mean) so dc_coarse/dtheta
#     stays a sparse matrix and the adjoint gradient of prototype A still applies.
#     NOTE: arithmetic mean models parallel flow; for barrier-like (series) flow
#     a harmonic mean is more accurate -- a production version should choose per
#     block based on within-block conductance contrast.
#
# INTEGRATION: produces a coarse (conductance, edges, demes) triple consumable by
# 01_block_solver_adjoint.R / terradish_graph. The aggregation matrix P
# (n_coarse x N) is the object terradish would differentiate through.
#
# Depends: Matrix; sources 01_block_solver_adjoint.R for the solver/forward.
# =============================================================================

suppressMessages(library(Matrix))
# source the Tier-1 engine (same folder). Adjust the path if running elsewhere.
if (!exists("forward_E")) source("01_block_solver_adjoint.R")

## ---- regular-grid helpers ----------------------------------------------------
grid_edges <- function(nr, nc, directions = 8L) {
  cid <- function(r, c) (r - 1L) * nc + c
  E_ <- list()
  for (r in 1:nr) for (cc in 1:nc) {
    if (cc < nc) E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r, cc + 1L))
    if (r < nr)  E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc))
    if (directions == 8L) {
      if (r < nr && cc < nc) E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc + 1L))
      if (r < nr && cc > 1)  E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc - 1L))
    }
  }
  do.call(rbind, E_)
}

## ---- per-edge current density from a forward solve --------------------------
# Sums |current| over all focal-site solves: J_e = w_e * sum_k |phi_ik - phi_jk|.
# phi columns are the (padded) potentials G with the grounded vertex set to 0.
edge_current_density <- function(c_vec, edges, G) {
  N <- length(c_vec)
  phi <- rbind(G, 0)                      # pad grounded vertex
  w <- c_vec[edges[, 1]] + c_vec[edges[, 2]]
  dphi <- abs(phi[edges[, 1], , drop = FALSE] - phi[edges[, 2], , drop = FALSE])
  w * rowSums(dphi)
}

## ---- build a coarsening (aggregation) matrix P : n_coarse x N ----------------
# block_id: length-N integer label; cells sharing a label merge into one coarse
# vertex. Protected cells must each carry a UNIQUE label so they never merge.
aggregation_matrix <- function(block_id) {
  labs <- sort(unique(block_id))
  remap <- match(block_id, labs)
  N <- length(block_id); nc <- length(labs)
  P <- sparseMatrix(i = remap, j = seq_len(N), x = 1, dims = c(nc, N))
  # row-normalise -> block MEAN (linear, differentiable)
  P / rowSums(P)
}

## ---- decide which native cells to coarsen -----------------------------------
# Returns a length-N block_id. Strategy: tile the grid into k x k super-cells;
# a super-cell collapses to one coarse vertex IFF none of its cells are focal,
# its max node current density is below the q-quantile, AND its internal
# conductance contrast is low (no barrier/corridor inside). Otherwise its cells
# keep native (unique) labels.
make_block_id <- function(nr, nc, k, focal_cells, node_current, c_vec,
                          current_q = 0.6, contrast_tol = 0.25) {
  N <- nr * nc
  block_id <- integer(N)
  next_lab <- 0L
  cur_thresh <- as.numeric(quantile(node_current, current_q))
  focal_set <- focal_cells
  for (br in seq(1L, nr, by = k)) for (bc in seq(1L, nc, by = k)) {
    rr <- br:min(br + k - 1L, nr); ccs <- bc:min(bc + k - 1L, nc)
    cells <- as.vector(outer((rr - 1L) * nc, ccs, "+"))
    has_focal   <- any(cells %in% focal_set)
    hi_current  <- max(node_current[cells]) > cur_thresh
    cc_vals     <- c_vec[cells]
    hi_contrast <- (max(cc_vals) - min(cc_vals)) / (mean(cc_vals) + 1e-9) > contrast_tol
    if (has_focal || hi_current || hi_contrast) {
      # keep native: unique label per cell
      block_id[cells] <- next_lab + seq_along(cells); next_lab <- next_lab + length(cells)
    } else {
      next_lab <- next_lab + 1L; block_id[cells] <- next_lab     # merge
    }
  }
  block_id
}

## ---- coarse graph from aggregation ------------------------------------------
# Coarse conductance = P %*% c (block mean). Coarse edges connect coarse vertices
# that contained adjacent fine cells; coarse focal vertex = label of focal cell.
coarsen_graph <- function(c_vec, edges, demes, P) {
  block_of <- apply(P, 2, function(col) which(col > 0)[1])  # fine cell -> coarse id
  c_coarse <- as.numeric(P %*% c_vec)
  ce <- cbind(block_of[edges[, 1]], block_of[edges[, 2]])
  ce <- ce[ce[, 1] != ce[, 2], , drop = FALSE]              # drop intra-block
  ce <- unique(t(apply(ce, 1, sort)))
  demes_coarse <- block_of[demes]
  list(c = c_coarse, edges = ce, demes = demes_coarse, block_of = block_of)
}

## ============================ self-test ======================================
if (sys.nframe() == 0L) {
  set.seed(7)
  nr <- 40L; nc <- 40L; N <- nr * nc
  edges <- grid_edges(nr, nc, 8L)
  # covariate field with a low-conductance vertical barrier (a corridor test)
  rr <- ((seq_len(N) - 1L) %/% nc) + 1L
  ccs <- ((seq_len(N) - 1L) %% nc) + 1L
  barrier <- ifelse(ccs == 20L, -3, 0)            # one-cell barrier at column 20
  smooth_field <- 0.3 * scale(rr) + 0.2 * scale(ccs)
  X <- cbind(as.numeric(smooth_field), barrier)
  cond_fun <- function(th) as.numeric(exp(X %*% th))
  theta <- c(0.5, 1.0)
  c_vec <- cond_fun(theta)
  demes <- c(1L * nc + 5L, 10L * nc + 35L, 30L * nc + 8L, 38L * nc + 33L)

  # fine solve
  fwd <- forward_E(theta, edges, demes, cond_fun)
  E_fine <- fwd$E

  # current density -> node current -> coarsening decision
  Je <- edge_current_density(c_vec, edges, fwd$G)
  node_current <- rowsum_fast(Je, edges[, 1], N) + rowsum_fast(Je, edges[, 2], N)
  block_id <- make_block_id(nr, nc, k = 4L, focal_cells = demes,
                            node_current = node_current, c_vec = c_vec,
                            current_q = 0.6, contrast_tol = 0.25)
  P <- aggregation_matrix(block_id)
  cg <- coarsen_graph(c_vec, edges, demes, P)

  # coarse solve (use the coarse conductance directly)
  Lc <- build_reduced_laplacian(cg$c, cg$edges)
  Nc <- length(cg$c)
  Gc <- block_pcg_jacobi(Lc, make_Z(Nc, cg$demes))$X
  E_coarse <- make_A(Nc, cg$demes) %*% Gc

  rel <- norm(as.matrix(E_fine - E_coarse), "F") / norm(as.matrix(E_fine), "F")
  cat(sprintf("fine vertices   : %d\n", N))
  cat(sprintf("coarse vertices : %d  (%.1f%% reduction)\n", Nc, 100 * (1 - Nc / N)))
  cat(sprintf("||E_fine - E_coarse||_F / ||E_fine||_F : %.4f\n", rel))
  cat("Interpretation: reduction comes from merging homogeneous low-current\n")
  cat("blocks; the column-20 barrier and focal neighbourhoods stay native, so\n")
  cat("E error stays small. Tighten current_q/contrast_tol to trade size vs error;\n")
  cat("a production loop refines until rel < tol (resolution convergence gate).\n")
}
