# Scaling `terradish` to Large, High-Resolution Landscapes

A ranked menu of approaches for estimating resistance/conductance over rasters
much larger than ~1M cells, while maximally leveraging the information content
of the data.

Author note: this document extends `future_dev_ideas.md` (2026-05-28). That
memo surveyed the option space well. This one does three additional things: it
states the bottleneck as a cost model so the options can be ranked on a common
scale; it separates and then re-couples the two distinct meanings of "maximal
information content"; and it backs the top-tier recommendations with verified
math and runnable prototypes (`dev/large_landscape/*.R`, `*.cpp`).

---

## 0. The bottleneck, stated as a cost model

Every likelihood evaluation in `terradish_algorithm()` does the following, on a
reduced graph Laplacian `L` of size `(N-1) x (N-1)` where `N` is the number of
active raster cells:

1. one **forward solve** `L G = Z`, where `Z` has `m` columns (one per focal
   site); `G` is `(N-1) x m`;
2. reads the focal covariance `E = A G` (`m x m`, the inferential target);
3. feeds `E` to the measurement model to get the objective and `dl/dE`;
4. for the **gradient**, one reverse-mode pass through `L`;
5. for the **exact Hessian**, `p` *additional* block solves (one per conductance
   parameter), assembled in `.terradish_algorithm_derivative_chunk()`.

Per outer optimisation step the dominant costs are therefore:

| Resource | Direct (CHOLMOD) | Cause |
|---|---|---|
| Factorisation time | `O(N^1.5)` (2-D nested dissection) | fill-in of the Cholesky factor |
| Factor memory | `O(N log N)` | fill-in; the factor is far denser than `L` |
| Solve time | `(1 + p) * m` back-substitutions | gradient + exact Hessian |
| Potential memory | `O(N m)` dense | `G` (and `tG`) must be materialised |

There are in fact **three** distinct walls, and they fail at different scales:

- **W1 — factorisation fill-in.** The classic ~1M-cell wall. This is what the
  AMG/PCG backends you already added (`src/amg_solver.cpp`, `solver = "amg"`,
  `auto` switching at 1.5M vertices, preconditioner reuse) are designed to
  break: an AMG V-cycle solve is `O(N)` and the hierarchy is reusable across
  nearby `theta`. **W1 is mostly solved in principle; what remains is making the
  iterative path the default, robust, and diagnostic.**
- **W2 — dense potential memory.** `G` is `(N-1) x m`. At `N = 5e7`, `m = 100`
  that is ~40 GB in double precision — independent of solver. This becomes the
  binding constraint once W1 is removed.
- **W3 — the `(1+p)`-solve exact Hessian.** Linear in the number of conductance
  parameters. Painful for spline / Gaussian-scale / many-covariate models.

A scaling solution has to name which wall each technique attacks. Most of the
survey items in `future_dev_ideas.md` attack W1 only; the highest-leverage moves
attack `N` itself (and thereby all three) or attack W2/W3 directly.

---

## 1. Two meanings of "maximally leveraging the information content"

This phrase has a numerical reading and a statistical reading. The strongest
solution serves both, and they turn out to reinforce each other.

**Numerical information — where resolution matters.** The data do not require
exact conductance at every cell; they require an accurate focal-site electrical
geometry `E`. Information about `E` is spatially concentrated: it lives where
current flows between focal pairs and where conductance changes sharply (narrow
corridors, hard barriers). Uniform high resolution spends compute where the
solution is smooth and the data are uninformative. *Maximally leveraging
information, numerically, means letting resolution follow the current density
and the conductance gradient* — refine where the solution carries information,
coarsen where it does not. This is an a-posteriori, solution-driven mesh, not an
arbitrary coarsening.

**Statistical information — using all of the genetic signal.** The current
default workflow (MLPE on a lower-triangular vector of pairwise distances with a
working correlation) discards information twice: it collapses the `m x m`
covariance to pairwise distances, and it usually operates on *population*-level
summaries. The `generalized_wishart` / `wishart_covariance` models you already
have are the full-information likelihood: they model the entire covariance with
an explicit effective degrees-of-freedom `nu`, and the Fisher information they
imply scales with `nu`. *Maximally leveraging information, statistically, means
(a) preferring the Wishart likelihood when marker counts are known, (b) using
the largest defensible `m` — including individual-level focal points, not just
populations — and (c) reporting the Fisher information explicitly.*

**The coupling.** These are not independent. A faster solver is not just about
fitting a bigger raster at fixed `m`; it is what lets you *raise* `m` (more focal
points = more statistical information) and afford the full-information Wishart
likelihood. And the Fisher-information curvature (Section 3, approach D) is
simultaneously the cheap replacement for the exact Hessian (attacks W3) *and*
the object that quantifies statistical information for inference. The numerical
refinement indicator (approach B) is, in the same spirit, an
information-density map. Scaling and information are the same problem viewed
from two sides.

---

## 2. Ranking criteria

Each approach below is scored on:

1. **Leverage** — which wall(s) it attacks (`N` itself > W1 > W2/W3).
2. **Target fidelity** — does it preserve `E` *and its exact derivatives*?
3. **Information retention** — does it keep fine-scale covariate signal and the
   full genetic signal, or discard them?
4. **Feasibility** — distance from the current architecture (graph + reverse-mode
   adjoint + SPD solvers).
5. **Inferential defensibility** — can approximation error be bounded/reported?

---

## 3. The ranked menu

### Tier 1 — build now

---

#### A. Matrix-free, block, preconditioner-reused iterative engine + streamed adjoint gradient  *(attacks W1, W2)*

**Idea.** Make the iterative solver the default for large graphs, and finish the
job so that no dense `O(N m)` object is ever held longer than necessary:

- Solve `L G = Z` as a **block** (multiple-RHS) PCG/AMG, not column-by-column.
  The `m` focal columns share one preconditioner and one Krylov space; block
  Krylov methods converge in fewer matvecs than `m` independent solves.
- **Reuse** the AMG hierarchy (already supported) *and* **warm-start** from the
  previous outer step's `G` (the optimizer moves `theta` slightly each step, so
  the previous potentials are excellent initial guesses — `solver_warm_start` is
  already plumbed through `terradish_algorithm()`).
- Compute the gradient by the **adjoint-state method**: it needs exactly **one**
  extra solve `L Lambda = A^T (dl/dE)`, *independent of the number of
  parameters*, and then a single streamed pass over edges to accumulate
  `dl/dc`. Tile that pass by raster block so peak memory is `O(N) + O(N m)` for
  the two potential blocks rather than several copies.

**Math (verified).** With `G = L^{-1} Z`, `E = A G`, and scalar loss
`l = f(E)`, the reverse-mode gradient is

```
dl/dG   = A^T (df/dE)
Lambda  = L^{-1} (dl/dG)               # one adjoint block solve
dl/dL   = -Lambda G^T                  # never formed densely; contracted per edge
dl/dw_e = (dl/dL)_{ii} + (dl/dL)_{jj} - (dl/dL)_{ij} - (dl/dL)_{ji}
dl/dc   = sum over incident edges (since w_e = c_i + c_j)
dl/dtheta = (dc/dtheta)^T dl/dc
```

`dev/large_landscape/verify.py` confirms this adjoint gradient matches a central
finite-difference gradient to a max relative error of **3e-9** on an 8-neighbour
grid, and that block Jacobi-PCG matches the direct solve to **1e-12**. The R
prototype `01_block_solver_adjoint.R` implements the streamed contraction.

**What it buys.** Removes the factorisation wall (W1) and caps memory at the two
potential blocks (W2 partially). The gradient cost stops scaling with `p`. This
is the floor everything else builds on, and most of it already exists in the
package — the missing pieces are block Krylov, cross-step warm starts wired into
the optimiser loop, the streamed edge contraction, and solver-error diagnostics.

**Cost / risk.** Medium. Likelihood and gradient become *approximate* at the
solver tolerance; tolerance must be tied to the optimiser's gradient tolerance
(loose early, tight near the optimum) or the line search can stall. Report
solver residuals as model diagnostics (your memo's "approximate large-landscape
mode" — agreed, and it should be mandatory, not optional).

**Scores.** Leverage W1/W2 · fidelity high (exact up to tolerance) · information
full · feasibility high · defensibility high (residuals reportable).

---

#### B. Adaptive, current-density–weighted multiresolution graph  *(attacks N itself → W1+W2+W3)*

**Idea.** The single largest win, and the most literal realisation of "leverage
the information content." Replace one-cell-one-vertex with a graph whose
resolution is driven by the solution:

1. Build a coarse base graph (e.g. aggregate native cells by 8x or 16x).
2. Solve once (Tier-1 engine). From the forward potentials, compute per-edge
   **current density** `J_e = w_e * |phi_i - phi_j|` summed over focal-pair
   solves, and the **conductance gradient magnitude** from the covariates.
3. **Refine** (split back toward native resolution) cells where `J_e` or the
   gradient is large, and **always** keep focal cells and one-cell-wide
   corridors/barriers at native resolution; keep homogeneous, low-current
   regions coarse.
4. Re-solve, re-estimate `theta`, and check **resolution convergence** of `E`:
   stop when `||E_fine - E_coarse|| / ||E_fine||` is below tolerance.

This is standard adaptive mesh refinement with an a-posteriori error indicator,
but the indicator is exactly the information-density map of Section 1: current
density is where the focal geometry is sensitive to conductance, i.e. where the
data are informative about `theta`.

**What it buys.** Order-of-magnitude reduction in vertex count for typical
landscapes (sparse focal sites, mostly-smooth covariates with localised
barriers). Because it shrinks `N`, it relieves all three walls at once. Fine
covariate detail is *retained where it matters*, so it does not throw away
numerical information the way uniform coarsening does.

**Cost / risk.** High implementation cost. Conductance aggregation over a coarse
cell is nontrivial (harmonic vs arithmetic mean of `c` matters: harmonic for
series/barrier-like flow, arithmetic for parallel — get this wrong and you bias
barrier resistance). Coarsening can sever a one-cell corridor — hence the
hard rule to protect corridors/barriers and to validate against a fine solve on
cropped regions. Derivatives must propagate through the aggregation operator
(prototype `02_adaptive_coarsen.R` keeps aggregation linear so `dc_coarse/dtheta`
stays a sparse matrix product).

**Scores.** Leverage highest (N) · fidelity controllable via refinement ·
information retained adaptively · feasibility medium-low · defensibility high
(resolution-convergence is a built-in diagnostic).

---

### Tier 2 — high value, compose with Tier 1

---

#### C. Focal-site Schur complement / domain decomposition  *(attacks W2, enables out-of-core)*

**Idea.** The target is `m x m` but exact computation touches all `N` cells.
Eliminate non-focal cells. Partition the raster into tiles; within each tile
eliminate interior cells to a boundary (interface) operator; assemble a small
global interface + focal problem; recover local quantities only where needed.
The exact reduced operator is the Schur complement
`S = L_ff - L_fu L_uu^{-1} L_uf`.

**What it buys.** Naturally tiled and out-of-core (pairs with `terra` block I/O
and `crop_to_focal_buffer()` which you already ship), parallel across tiles,
and memory scales with the interface, not the interior. Best when focal sites
are sparse in a huge extent — a common real case.

**Cost / risk.** Exact Schur complements are dense on the interface; you need an
approximate/iterative interface solve, and the gradient must propagate through
it (doable: the Schur complement is differentiable, but bookkeeping is heavier
than approach A). Poor partitioning can bias resistance paths across cuts.

**Scores.** Leverage W2 + parallelism · fidelity high if interface solve tight ·
information full · feasibility medium · defensibility medium.

---

#### D. Gauss–Newton / Fisher-information curvature instead of the exact Hessian  *(attacks W3, and is the inference object)*

**Idea.** Replace the `(1+p)`-solve exact Hessian with the Gauss–Newton /
Fisher-information approximation. For the Wishart and least-squares measurement
models the objective is (a transform of) a Gaussian/Wishart negative
log-likelihood, so the expected Hessian is the **Fisher information**

```
I(theta) = J^T  W  J ,   J = dE/dtheta (vectorised), W = measurement-model weight
```

where the columns of `J` are obtained by the *same* adjoint machinery (one solve
each, but you can also estimate `tr`-type terms by **Hutchinson** stochastic
probing rather than forming `J` fully). `verify.py` confirms the Hutchinson
estimator of `tr(L^{-1})` to ~1% at 2000 probes; for curvature you typically
need far fewer probes because you only need a few leading directions.

**What it buys.** Three things at once: (i) drops the parameter-linear Hessian
cost (W3); (ii) is positive semidefinite by construction, so the
Newton/box-constrained step is always a descent direction (more robust than the
exact Hessian, which can be indefinite far from the optimum); (iii) **is** the
asymptotic information matrix, so `vcov`, standard errors, and the `nu`-scaling
of confidence intervals come directly from it — the literal "information content
of the data." This is the cheapest high-value item on the list and it composes
with everything above.

**Cost / risk.** Low. Gauss–Newton curvature differs from the true Hessian away
from the optimum (slightly more outer iterations possible); negligible for
well-specified models. Keep the exact Hessian available for small problems as a
check (it already exists).

**Scores.** Leverage W3 · fidelity exact at optimum (GN ≈ Hessian) · information
*is* the Fisher matrix · feasibility high · defensibility high.

---

### Tier 3 — research / longer-term

- **E. Continuum FEM/FVM reformulation** of `∇·(c∇u)=b` with a genuine adaptive
  mesh. The cleanest long-term answer (decouples mesh from raster entirely; huge
  numerical-methods literature) but the deepest rewrite and needs careful
  validation against the graph formulation. Best framed as a separate backend.
- **F. Randomized / landmark / low-rank resistance embeddings** (Nyström,
  Johnson–Lindenstrauss, Hutchinson effective-resistance estimators). Excellent
  for *model screening* and for initialising the optimiser; defensible for final
  inference only after substantial validation. Note: Hutchinson effective-
  resistance estimation is the same primitive as approach D — build it once.
- **G. Hierarchical / H-matrix compression** of the focal Green's function.
  High implementation complexity; benefit depends on exploitable low-rank
  structure under heterogeneous conductance. Longest-term.
- **H. GPU / compiled smoothing.** Real wins for the FFT Gaussian-scale
  smoothing and for iterative-solver matvecs, but secondary: it accelerates
  kernels, not the structural `N` problem. Do it after A–D.

---

## 4. How they compose (a suggested build order, not a single mandate)

The menu is not mutually exclusive; the practical recommendation is a stack:

```
                +-----------------------------------------------+
  inference     |  D. Gauss–Newton / Fisher curvature           |  → vcov, SEs, model selection
                +-----------------------------------------------+
  structure     |  B. adaptive current-density multiresolution  |  → shrinks N
                +-----------------------------------------------+
  solve         |  A. matrix-free block iterative + adjoint      |  → breaks W1/W2 on whatever graph B hands it
                +-----------------------------------------------+
  data /        |  C. Schur / tiling for out-of-core when needed |
  out-of-core   |  + Wishart likelihood, individual-level m      |  → full statistical information
                +-----------------------------------------------+
```

Recommended sequence by effort-to-payoff:

1. **D first** (days, not weeks): cheapest, improves robustness *and* inference
   immediately, and is needed by everything as the curvature/uncertainty object.
2. **A**: finish the iterative engine (block Krylov + warm start + streamed
   adjoint + mandatory solver diagnostics). This is mostly hardening work on
   code you already have.
3. **B**: the structural prize. Prototype on cropped regions with a
   resolution-convergence gate before integrating.
4. **C** only when focal sites are sparse in extents too large to hold even
   coarsened — then go out-of-core.

If only one thing is built this cycle, build **D**; if the goal is genuinely
continental-scale rasters, **B** is the one that changes the asymptotics.

---

## 5. Validation plan

Any approximation must be judged on ecological/statistical accuracy, not just
runtime. The decisive experiments:

1. **Exact-vs-approximate `E`** on cropped regions small enough for a direct
   CHOLMOD solve — the ground truth already in the package.
2. **Gradient check**: adjoint vs finite difference on small graphs (the
   `verify.py` harness; port to `testthat`).
3. **Parameter recovery under simulation.** This is where your local **SLiM**
   install earns its keep: simulate genotypes under a *known* conductance
   surface on a large landscape, then test whether the scalable pipeline
   recovers the generating `theta` and whether its confidence intervals have
   nominal coverage. SLiM gives a defensible, mechanism-based truth that an
   analytic resistance simulator cannot.
4. **Resolution convergence** (approach B): `theta` and `E` as the adaptive mesh
   refines; report the convergence curve.
5. **Model-ranking stability**: does AIC / CV ranking change under
   approximation? (ranking stability matters more than absolute likelihood).
6. **Barrier / corridor stress tests**: one-cell corridors and hard barriers are
   where coarsening fails; test them explicitly.
7. **Focal-density tests**: sparse, clustered, and broadly distributed `m`.

---

## 6. Prototypes in this folder

- `verify.py` — standalone numerical proof of the three core identities
  (adjoint gradient ≈ FD to 3e-9; block-PCG ≈ direct to 1e-12; Hutchinson
  trace ≈ truth). Run: `python3 verify.py`.
- `01_block_solver_adjoint.R` — matrix-free block PCG on the reduced Laplacian +
  streamed adjoint gradient, with a self-test against a dense pseudo-inverse on a
  small grid. Approach **A**.
- `02_adaptive_coarsen.R` — current-density / gradient-driven adaptive
  coarsening prototype with focal- and corridor-preservation and an
  `E`-resolution-convergence checker. Approach **B**.
- `03_gauss_newton_fisher.R` — Gauss–Newton / Fisher-information curvature and a
  Hutchinson trace primitive, returning the information matrix and standard
  errors. Approach **D**.
- `hutchinson_diag.cpp` — optional Rcpp stochastic diagonal/trace estimator for
  `diag(L^{-1})` and curvature probing, for when the R prototype's per-probe
  solve is the bottleneck.

These are research prototypes: correct on the math (verified) and structured to
mirror `terradish` conventions (`w_e = c_i + c_j`, last-vertex grounding,
`E = A G`), but not yet wired into `terradish_algorithm()`'s dispatch. Each file
documents its integration points.
