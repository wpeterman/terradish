"""
Scenario 4 (DRAGON rung-3 validation): asymmetric migration + a density/Ne
gradient correlated with elevation, via the recapitation pipeline.

  SLiM (tree-seq: downhill migration + size gradient)
    -> pyslim.recapitate -> msprime.sim_mutations -> sample per deme
    -> per-deme covariance
    -> DRAGON fits: symmetric vs direction-only vs direction+drift.

Expected (the confound test on real coalescent data, mirroring the rung-2 sweep):
  - direction is detected (direction beats symmetric);
  - a drift/coalescence-rate gradient is detected (joint beats direction-only);
  - the DIRECTION-ONLY coefficient is biased relative to the JOINT coefficient,
    because the omitted density gradient runs along the directional covariate.

Run from the terradish repo root with the conda `slim` Python, e.g.
  & "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe" \
       dev/slim/run_scen4_recap.py
Environment overrides: SLIM (slim.exe path), SLIM_FORCE=TRUE (re-run SLiM),
  DIM, N0, MIG, BETA_DIR, KAPPA_N, BURNIN, NSAMP, SEED.
"""
import os, sys, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))          # .../dev/slim
REPO = os.path.dirname(os.path.dirname(HERE))               # repo root
PY = sys.executable                                         # this (conda) python
SLIM = os.environ.get("SLIM", "C:/msys64/mingw64/bin/slim.exe")

OUTREL = "dev/slim/scen4ts"
OUTPRE = os.path.join(REPO, OUTREL)
SLIM_SCRIPT = "dev/slim/scen4_ts.slim"
RECAP = "dev/slim/recap_sample.py"

DIM      = int(os.environ.get("DIM", 6))
N0       = int(os.environ.get("N0", 200))
MIG      = os.environ.get("MIG", "0.004")
BETA_DIR = os.environ.get("BETA_DIR", "0.8")
KAPPA_N  = os.environ.get("KAPPA_N", "0.6")
BURNIN   = os.environ.get("BURNIN", "6000")
NSAMP    = int(os.environ.get("NSAMP", 25))
SEED     = int(os.environ.get("SEED", 4004))
FORCE    = os.environ.get("SLIM_FORCE", "FALSE").upper() == "TRUE"

sys.path.insert(0, HERE)
from dragon_engine import fit_dragon

def run(cmd):
    print("  $", " ".join(str(c) for c in cmd))
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    print("   " + "\n   ".join(out.strip().splitlines()[-3:]))
    if r.returncode != 0:
        raise SystemExit(f"command failed ({r.returncode})")

if FORCE or not os.path.exists(OUTPRE + "_Y.csv"):
    print("1) SLiM tree-sequence (asymmetric migration + density gradient) ...")
    run([SLIM, "-seed", SEED,
         "-d", f"DIM={DIM}", "-d", f"N0={N0}", "-d", f"MIG={MIG}",
         "-d", f"BETA_DIR={BETA_DIR}", "-d", f"KAPPA_N={KAPPA_N}",
         "-d", "RHO=1e-8", "-d", "L=1000000", "-d", f"BURNIN={BURNIN}",
         "-d", f"OUTPRE='{OUTREL}'", SLIM_SCRIPT])

    demes = np.genfromtxt(OUTPRE + "_demes.csv", delimiter=",", names=True)
    anc = int(np.sum(demes["size"]))                       # ancestral Ne ~ total census size
    print(f"2) recapitate (ancestral_Ne={anc}) + mutations + sample ...")
    run([PY, RECAP, "--trees", OUTPRE + ".trees", "--out", OUTPRE,
         "--ancestral_Ne", anc, "--recombination_rate", "1e-8",
         "--mutation_rate", "1e-8", "--nsample", NSAMP, "--seed", SEED,
         "--mode", "population"])

Y = np.loadtxt(OUTPRE + "_Y.csv", delimiter=",")
N = np.loadtxt(OUTPRE + "_N.csv", delimiter=",")
demeids = np.loadtxt(OUTPRE + "_demeids.csv", delimiter=",").astype(int)
demes = np.genfromtxt(OUTPRE + "_demes.csv", delimiter=",", names=True)
elev_by = demes["elev"][np.searchsorted(demes["deme"], demeids)]

print("3) DRAGON fits (symmetric / direction-only / direction+drift) ...")
o = fit_dragon(Y, N, demeids, DIM, elev_by)

print("\n=============== Scenario 4: DRAGON rung-3 ===============")
print(f"observed demes = {o['n_obs']}   markers (MAF>=0.05) = {o['p']}")
print(f"  NLL  symmetric        = {o['nll_sym']:.2f}")
print(f"  NLL  direction-only   = {o['nll_dir']:.2f}")
print(f"  NLL  direction+drift  = {o['nll_joint']:.2f}")
print(f"  LRT direction vs sym       = {o['LRT_dir_vs_sym']:6.2f}  (>3.84 => direction detected)")
print(f"  LRT drift  (joint vs dir)  = {o['LRT_joint_vs_dir']:6.2f}  (>3.84 => drift gradient detected)")
print(f"  gamma_dir  direction-only  = {o['gdir_dironly']:+.3f}")
print(f"  gamma_dir  joint           = {o['gdir_joint']:+.3f}   <- compare to direction-only")
print(f"  kappa      joint (drift)   = {o['kappa_joint']:+.3f}")
bias_shift = abs(o['gdir_dironly']) - abs(o['gdir_joint'])
print(f"\n  confound check: |gdir| shifts by {bias_shift:+.3f} when the drift surface is added.")
print("  expected: direction detected; drift detected; direction-only |gdir| inflated vs joint.")
print("========================================================")
