# Turnkey full-scale validation (run on your hardware)

Two parameterized scripts to validate terradish at scales beyond what the dev
sandbox can reach. Both require an installed `terradish` (so `library(terradish)`
works) and `terra`.

## 1. Computational scaling -- `benchmark_solvers.R`
Times one likelihood+gradient evaluation across raster sizes and solvers on a
synthetic heterogeneous landscape (no genetic simulation needed). This is where
you push past 1M cells and compare the direct factorization against AMG.

```bash
Rscript benchmark_solvers.R --sides 256,512,1024,1448,2048 --m 50 \
        --solvers direct,amg,pcg,block_cg --out solver_scaling.csv
```
- `direct` is fastest until Cholesky fill-in exhausts RAM (watch memory past
  ~1M cells); `amg` is the intended large-N path; `block_cg`/`pcg` are
  iterative fallbacks for small/well-conditioned graphs.
- Raise `--m` to chart cost vs the number of focal sites at fixed N.
- Output columns: n, cells, m, solver, eval_s, iters, converged.

## 2. Statistical recovery + information -- `recovery_study.R`
Simulates genotypes under a known conductance surface (SLiM + tree sequences),
then fits terradish while increasing the number of focal sites `m` (averaged
over random site subsets) and the marker count `nu`.

Set the external tools first (Windows examples):
```bash
# PowerShell
$env:SLIM_BIN   = "C:\msys64\mingw64\bin\slim.exe"
$env:PYTHON_BIN = "C:\Users\peterman.73\AppData\Local\anaconda3\envs\slim\python.exe"
$env:SLIM_DIR   = "<...>/terradish/dev/large_landscape/slim"
Rscript recovery_study.R --scenario barrier --nx 64 --n_focal 80 \
        --K 8000 --sigma_d 0.012 --nticks 300 \
        --m 10,20,40,80 --nreps 10 --nu 500,1000,2000,4000 --out recovery.csv
```
- `--sigma_d` should be **smaller than the barrier width** or the dispersal
  kernel jumps the barrier and it leaves little genetic imprint (we observed
  this: a 1-cell barrier with sigma_d ~= 2 cells is effectively permeable).
- `--nreps` averages over random site subsets so the m curve isn't dominated by
  which sites happen to be drawn.
- Expected pattern: `se_mean` falls with `m` and as `~1/sqrt(nu)`; point
  estimates are invariant to `nu`; `cos_mean` (recovery vs truth) is stable and
  high when the simulated effect genuinely impedes gene flow.

The Python genetics deps (in your `slim` conda env): `msprime pyslim tskit numpy pillow`.
