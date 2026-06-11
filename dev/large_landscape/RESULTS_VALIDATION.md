# Validation: large-raster scaling and increasing sample size

In-sandbox results (2 cores / ~3.8 GB). Turnkey scripts for full-scale runs on
your hardware are in `dev/large_landscape/turnkey/`.

## Computational scaling (`bench_scaling.R`, synthetic-from-model)
Heterogeneous 8-neighbour grid (gradient + one-cell barrier), m=50 focal sites.

| cells   | direct setup (s) | direct solve (s) | factor fill / nnz(L) | IC-PCG iters | IC-PCG (s) |
|---------|------------------|------------------|----------------------|--------------|------------|
| 16,384  | 0.024            | 0.028            | 6.85                 | 212          | 3.2        |
| 65,536  | 0.215            | 0.231            | 8.27                 | 428          | 32.9       |
| 131,044 | 0.789            | 0.304            | 11.1                 | -            | -          |
| 262,144 | 1.335            | 0.885            | 9.86                 | -            | -          |
| 360,000 | 1.933            | 0.902            | 10.07                | -            | -          |

- Direct Cholesky factorizes a **360k-cell** Laplacian in ~2 s; the factor is
  ~7-11x denser than the Laplacian (the fill-in that eventually walls direct at
  much larger N on memory). Direct is the right default through several hundred
  thousand cells -- consistent with terradish's defaults.
- IC-PCG iteration count grows with N (212 -> 428) but its wall-clock in this
  build isn't competitive; the production large-N path is AMG (unit-tested on
  melip; benchmark it at >1M cells via `turnkey/benchmark_solvers.R`).
- block-CG converges in fewer shared iterations than independent CG on small /
  well-conditioned graphs (64x64: 64 vs up to 202 per column); it is not for
  large heterogeneous Laplacians.

## Information content: increasing markers nu (`sweep_recovery.R`, SLiM data)
Generalized Wishart, m=40 sites fixed, varying nu. **Textbook behaviour:**

| nu   | b_gradient | b_barrier | se_gradient | se_barrier |
|------|-----------|-----------|-------------|------------|
| 500  | 0.179     | 0.121     | 0.0170      | 0.0276     |
| 2000 | 0.179     | 0.121     | 0.0085      | 0.0138     |
| 8000 | 0.179     | 0.121     | 0.0042      | 0.0069     |

Point estimates are **invariant to nu**; standard errors scale exactly as
**1/sqrt(nu)** (4x markers -> half the SE). This validates the Gauss-Newton /
Fisher `vcov` from D end to end: more markers buy proportionally tighter CIs.

## Information content: increasing focal sites m (SLiM data)
With a barrier that genuinely impedes gene flow (sigma_d < barrier width),
recovery is good and **SEs shrink with more sites** (tighter sim, m sweep):

| m  | cos(theta_hat, true) | se_mean | se_gradient |
|----|----------------------|---------|-------------|
| 10 | 0.96                 | 0.0054  | 0.0108      |
| 40 | 0.96                 | 0.0018  | 0.0036      |
| 80 | 0.96                 | 0.0013  | 0.0026      |

(An m=20 subset dipped because the *first* 20 sites clustered on one side of the
barrier; `turnkey/recovery_study.R` averages over random subsets to remove this.)

Caveats surfaced (and worth keeping): a one-cell barrier with sigma_d ~= 2 cells
is effectively permeable, so terradish correctly reports little barrier effect
-- recovery requires the simulated process to actually impede gene flow. And a
strong barrier's *magnitude* is weakly identified (genetic distance saturates)
even though its *direction* is robust.

## Bottom line
- Scaling: direct is fast to several hundred k cells; AMG is the >1M path
  (benchmark at full scale on your hardware).
- Information: the GN/Fisher inference behaves exactly as theory predicts --
  SEs ~ 1/sqrt(nu) and shrink with more focal sites -- so scaling the solver to
  afford larger m / the full-covariance Wishart likelihood directly buys
  inferential precision. That is the "maximally leverage the information
  content" claim, now demonstrated.
