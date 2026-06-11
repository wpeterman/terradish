#!/usr/bin/env python
"""Scenario 5: the realistic end-to-end fair test for the FRAME-coupled model.

Effective size is set to genuinely track the migration-stationary distribution:
K_i proportional to pi_i^k, where pi is the stationary distribution of the directed
generator M_ab = exp(BETA*(elev_a - elev_b)) on the 5x5 rook lattice. This is the
coupled model's own assumption (gamma_i = c * pi_i^(-alpha)), so a coupled fit
should recover BOTH the directional coefficient (~BETA) and alpha (~k) -- unlike
scen4, where size was an independent elevation gradient and alpha pinned at the
bound.

  pi (from M) -> K_i ~ pi_i^k -> SLiM tree-seq (directional migration + those sizes)
    -> pyslim.recapitate -> msprime.sim_mutations -> sample per deme
    -> Y/N/demeids/demes CSVs consumed by validate_scen5.R

Run from the terradish repo root with the conda `slim` python (no scipy needed):
  & ".../envs/slim/python.exe" dev/slim/run_scen5_recap.py
"""
import argparse, json, os, subprocess, sys, time


def _bootstrap_conda_dlls():
    if os.name != "nt":
        return
    dirs = [sys.prefix, os.path.join(sys.prefix, "Library", "bin"),
            os.path.join(sys.prefix, "Library", "mingw-w64", "bin"),
            os.path.join(sys.prefix, "Library", "usr", "bin"),
            os.path.join(sys.prefix, "Scripts")]
    dirs = [d for d in dirs if os.path.isdir(d)]
    os.environ["PATH"] = os.pathsep.join(dirs) + os.pathsep + os.environ.get("PATH", "")
    for d in dirs:
        try:
            os.add_dll_directory(d)
        except OSError:
            pass


_bootstrap_conda_dlls()
import numpy as np
import msprime, pyslim, tskit

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SLIM = r"C:\msys64\mingw64\bin\slim.exe"


def stationary_pi(grid, beta):
    """pi of the directed generator M_ab = exp(beta*(elev_a-elev_b)), elev = column."""
    d = grid * grid
    cx = np.arange(d) % grid
    cy = np.arange(d) // grid
    elev = cx.astype(float)
    M = np.zeros((d, d))
    for a in range(d):
        for b in range(d):
            if a != b and abs(cx[a]-cx[b]) + abs(cy[a]-cy[b]) == 1:
                M[a, b] = np.exp(beta * (elev[a] - elev[b]))
    L = np.diag(M.sum(1)) - M
    pi = np.linalg.solve(L.T + np.ones((d, d)), np.ones(d))   # (L'+11') pi = 1
    pi = pi / pi.sum()
    return pi, elev, cx, cy


def recap_mutate_sample(trees, grid, mu, rho, n_dip, seed, outpre):
    ts = tskit.load(trees)
    ndeme = grid * grid
    anc = int(round(ts.num_samples / 2))
    rts = pyslim.recapitate(ts, ancestral_Ne=anc, recombination_rate=rho, random_seed=seed)
    mts = msprime.sim_mutations(rts, rate=mu, model=msprime.BinaryMutationModel(),
                                random_seed=seed)
    print(f"  [recap] anc_Ne={anc}  sites={mts.num_sites}", flush=True)
    rng = np.random.default_rng(seed)
    by = {i: [] for i in range(ndeme)}
    for iid in pyslim.individuals_alive_at(mts, 0):
        ind = mts.individual(iid)
        by[mts.node(ind.nodes[0]).population].append(ind.nodes)
    nodes, col_deme = [], []
    for dm in range(ndeme):
        inds = by[dm]
        if len(inds) < n_dip:
            raise RuntimeError(f"deme {dm}: {len(inds)} inds < n_dip {n_dip}")
        for k in rng.choice(len(inds), size=n_dip, replace=False):
            nodes.extend(inds[k]); col_deme.extend([dm, dm])
    nodes = np.asarray(nodes, np.int32); col_deme = np.asarray(col_deme)
    G = mts.genotype_matrix(samples=nodes)
    Y = np.zeros((ndeme, G.shape[0]), np.int64)
    for dm in range(ndeme):
        Y[dm] = G[:, col_deme == dm].sum(1)
    Nhap = np.full(ndeme, 2 * n_dip, np.int64)
    tot = Y.sum(0); keep = (tot > 0) & (tot < Nhap.sum())
    Y = Y[:, keep]
    print(f"  [sample] {n_dip} dip/deme, {Y.shape[1]} segregating loci", flush=True)
    np.savetxt(outpre + "_Y.csv", Y, fmt="%d", delimiter=",")
    np.savetxt(outpre + "_N.csv", Nhap, fmt="%d", delimiter=",")
    np.savetxt(outpre + "_demeids.csv", np.arange(ndeme), fmt="%d", delimiter=",")
    return Y, Nhap


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slim", default=os.environ.get("SLIM_BIN", DEFAULT_SLIM))
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--grid", type=int, default=5)
    ap.add_argument("--beta", type=float, default=0.5, help="directional truth (gamma_dir)")
    ap.add_argument("--k", type=float, default=0.5, help="coupling exponent: K ~ pi^k (alpha truth)")
    ap.add_argument("--k0", type=int, default=300, help="carrying capacity at pi=mean")
    ap.add_argument("--kmin", type=int, default=30)
    ap.add_argument("--m0", type=float, default=0.05)
    ap.add_argument("--ngen", type=int, default=10000)
    ap.add_argument("--L", type=float, default=1e6)
    ap.add_argument("--rho", type=float, default=1e-8)
    ap.add_argument("--mu", type=float, default=1e-7)
    ap.add_argument("--ndip", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.k0, args.ngen, args.L, args.mu, args.ndip = 120, 1500, 5e5, 2e-7, 15
    os.makedirs(args.out, exist_ok=True)
    outpre = os.path.join(args.out, "scen5_ts")

    pi, elev, cx, cy = stationary_pi(args.grid, args.beta)
    K = np.maximum(np.round(args.k0 * (pi / pi.mean()) ** args.k).astype(int), args.kmin)
    print(f"truth: BETA(gamma_dir)={args.beta}  k(alpha)={args.k}", flush=True)
    print(f"pi range [{pi.min():.4g},{pi.max():.4g}]  K range [{K.min()},{K.max()}]  "
          f"cor(log pi, elev)={np.corrcoef(np.log(pi), elev)[0,1]:.3f}", flush=True)

    sizes_str = "c(" + ",".join(str(int(x)) for x in K) + ")"
    cmd = [args.slim, "-seed", str(args.seed), "-d", f"GRID={args.grid}",
           "-d", f"M0={args.m0}", "-d", f"BETA={args.beta}", "-d", f"SIZES={sizes_str}",
           "-d", f"NGEN={args.ngen}", "-d", f"L={args.L}", "-d", f"RHO={args.rho}",
           "-d", f"OUTPRE='{outpre.replace(os.sep, '/')}'",
           os.path.join(HERE, "scen5_ts.slim")]
    print("1) SLiM ...", flush=True)
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    for ln in r.stdout.splitlines():
        if ln.startswith("scen5:"):
            print("   " + ln, flush=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + "\n" + r.stderr + "\n"); raise SystemExit("SLiM failed")
    print(f"   done in {time.time()-t0:.1f}s", flush=True)

    print("2) recapitate + mutate + sample ...", flush=True)
    recap_mutate_sample(outpre + ".trees", args.grid, args.mu, args.rho,
                        args.ndip, args.seed, outpre)

    # demes table with the true pi (for reference / the coupled cross-check)
    with open(outpre + "_demes.csv", "w") as fh:   # overwrite SLiM's with pi + dens cols
        fh.write("deme,gx,gy,elev,dens,pi,size\n")
        for i in range(args.grid * args.grid):
            fh.write(f"{i},{cx[i]},{cy[i]},{elev[i]:.0f},{elev[i]:.0f},{pi[i]:.6f},{K[i]}\n")
    with open(outpre + "_truth.json", "w") as fh:
        json.dump(dict(beta=args.beta, k=args.k, pi=pi.tolist(), K=K.tolist()), fh, indent=2)
        fh.write("\n")
    print(f"\nOutputs in {args.out} (prefix scen5_ts). Now: Rscript dev/slim/validate_scen5.R",
          flush=True)


if __name__ == "__main__":
    main()
