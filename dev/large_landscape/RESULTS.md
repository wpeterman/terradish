# Validation results: D, A, B on SLiM-simulated data

What was built and run end-to-end (all reproducible in this folder), and what
the numbers say. Methods refer to the ranked menu in `SCALING_DESIGN.md`.

## Toolchain (stood up and executed)
- **SLiM 5.2** (built from source) runs the user's nonWF continuous-space
  template (`slim/landscape.slim`) with a habitat/conductance map.
- **tree sequences + msprime/pyslim/tskit** recapitate, overlay neutral
  mutations, sample focal individuals, and emit genetic distance/covariance
  matrices (`slim/process_trees.py`).
- **R 4.5** (Matrix) runs the reference implementation and the prototypes.
- Three scenarios generated: `smooth`, `barrier`, `corridor`
  (64x64 covariate grids, K≈6000 individuals, N≈1000–1600 alive, ν=4000 SNPs,
  m=20 focal sites). Data live in `slim/scenarios/<name>/`.

All three are run by `validate_DAB.R` against a faithful reference model
(reduced Laplacian `w_e=c_i+c_j`, last-vertex grounding, `E = A G`, log-linear
conductance, least-squares resistance-distance measurement model).

## Parameter recovery (sanity that the data carry the signal)
On the **barrier** scenario, the fitted conductance coefficients align with the
SLiM ground truth: cosine(θ̂, θ_true) = **0.96** at native resolution and
**0.96** at 2× resolution. Magnitudes differ (log-conductance is identifiable
only up to scale, which the LS slope absorbs), but the ecology is recovered:
gradient effect positive, barrier effect strongly negative. Recovery is stable
across resolution — a resolution-convergence check in itself.

## A — matrix-free block iterative solve + streamed adjoint
- **Adjoint gradient** matches finite differences to **1.9e-8** (R) / **3.3e-9**
  (Python); the full gradient is one forward + one adjoint solve, independent of
  the number of parameters.
- **Block-PCG reproduces the direct-solve E** to **2–3e-11** — identical up to
  tolerance. The iterative path is correct.
- **Why the direct solver hits a wall (W1), shown directly.** Cholesky factor
  fill-in grows with resolution on the barrier landscape:

  | factor | cells | direct solve (s) | factor nnz | fill / nnz(L) |
  |---|---|---|---|---|
  | 1× | 4,096 | 0.003 | 108,693 | 5.4 |
  | 2× | 16,384 | 0.023 | 555,737 | 6.9 |
  | 3× | 36,864 | 0.076 | 1,412,262 | 7.7 |
  | 4× | 65,536 | 0.167 | 2,696,048 | 8.3 |

  The factor is already 5–8× denser than the Laplacian and the ratio climbs with
  N — the mechanism behind the ~1M-cell ceiling. The iterative engine avoids
  forming this factor at all.
- **Cached symbolic factor (terradish's `reused_factor_template`)** cut the
  native-resolution fit from **5.8 s → 0.4 s** (~15×) by reusing the
  fill-reducing ordering across optimiser steps. This is a free, immediate win
  in the current direct path.

## D — Gauss-Newton / Fisher curvature
- **GN converges to the exact Hessian as the model fits well**, on a clean
  sum-of-squares check with model-generated data:
  rel‖I_GN − H‖/‖H‖ = **0.038** at noise sd 0.05 and **0.0036** at sd 0.005.
  This validates the curvature claim.
- **Cost**: GN curvature is essentially free (it reuses the gradient Jacobian),
  versus the exact Hessian's extra per-parameter solves. The gap grows with the
  number of conductance parameters.
- **Inference**: Fisher information yields θ standard errors directly
  (e.g. barrier fit SEs ≈ 0.017 and 0.97 at ν=4000), and is positive
  semidefinite by construction → robust box-constrained Newton steps.
- **Honest caveat**: on the real SLiM data the GN-vs-exact-Hessian gap is ~0.7.
  That is **model misspecification** (a linear log-conductance model cannot
  perfectly reproduce SLiM's genetic distances, so residuals are non-negligible
  at the optimum), not a defect of GN. GN still gives a valid descent direction
  and valid asymptotic SEs; the gap itself is a useful misspecification signal.

## B — adaptive current-density multiresolution (the cautionary result)
This is where validation changed the recommendation. Naive node-aggregation —
even with correct Galerkin edge-weight lumping — **biases resistance distance**,
because collapsing an aggregate to a single node removes the internal path
resistance that resistance distance integrates over. On the barrier scenario:

| current_q (protect) | coarse N | reduction | rel E error |
|---|---|---|---|
| 0.85 | 1,516 | 63.0% | 0.46 |
| 0.60 | 2,236 | 45.4% | 0.34 |
| 0.35 | 3,046 | 25.6% | 0.21 |
| 0.15 | 3,721 |  9.2% | 0.07 |

A clean, monotone accuracy/size tradeoff (the resolution-convergence gate works),
but you only reach <10% E error at <10% node reduction. The same pattern holds
on the smooth scenario, so the bias is **not** primarily a barrier artifact — it
is intrinsic to lumping.

**Revised conclusion for B.** Pure coarsening is a *screening / coarse-stage*
tool with a mandatory convergence gate, not a path to low-error large-N
inference on its own. For accurate scalable inference, prefer **A** (iterative
solve on the full graph, exact up to tolerance) and, when N is too large to hold
even coarsened, **tiled Schur-complement elimination (C)**, which removes
non-focal cells *exactly* rather than lumping them. In the ranking, this nudges
C up relative to B: B feeds A/C a good starting point and identifies where
refinement is needed, but the exact resistance must come from A or C.

## Reproduce
```bash
# 1) generate + run a scenario (SLiM + tree seq + genetics)
cd slim
python3 make_landscape.py --scenario barrier --nx 64 --ny 64 --n_focal 20 --out scenarios/barrier
slim -d "MAP='scenarios/barrier/habitat.png'" -d K=6000 -d SIGMA_C=0.035 -d SIGMA_M=0.035 \
     -d SIGMA_D=0.035 -d NTICKS=200 -d L=100000000 -d "OUTPATH='/tmp/out.trees'" landscape.slim
python3 process_trees.py --trees /tmp/out.trees --truth scenarios/barrier/truth.npz \
     --mu 5e-9 --recomb 1e-8 --Ne 2000 --n_per_site 8 --out scenarios/barrier/genetics.npz
# 2) validate D / A / B against the resulting data
cd ..
Rscript validate_DAB.R slim/scenarios/barrier
```

## Bottom line
- **A** and **D** are validated and ready to harden into terradish: A removes
  the fill-in wall and the parameter-linear gradient cost; D removes the
  parameter-linear Hessian cost and delivers the inference (SEs) as a byproduct.
  Both are low-risk, high-value, and largely build on code already in the
  package (AMG/PCG backends, reused factor template, reverse-mode adjoint).
- **B** is real but limited: keep it as an adaptive screening/refinement stage
  with a convergence gate; route final inference through A (exact-iterative) or
  C (tiled Schur). The validation moved C up the priority list.
- SLiM-based parameter-recovery and CI-coverage testing is the right standing
  validation harness and is now wired up; extend it to larger landscapes and to
  the `generalized_wishart` likelihood when integrating into the package core.
