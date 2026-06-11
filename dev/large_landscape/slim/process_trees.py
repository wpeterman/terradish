#!/usr/bin/env python3
"""
process_trees.py — turn a SLiM tree sequence into terradish-ready genetic data.

Steps:
  1. load the .trees, recapitate (msprime) so all lineages coalesce,
  2. overlay neutral mutations at rate --mu,
  3. take individuals alive at the end; assign the --n_per_site nearest to each
     focal coordinate (from truth.npz),
  4. compute, among focal sites:
        * allele-frequency covariance matrix  COV   (for wishart_covariance)
        * a genetic DISTANCE matrix           S     (mean per-SNP squared
          allele-frequency difference; an Fst-like distance),
        * nu = number of polymorphic SNPs retained (effective d.f.),
  5. write genetics.npz AND plain CSVs (S, COV, X, focal_cell, dims, theta_true)
     so the R harness can read them without an npz reader.

Usage:
  python3 process_trees.py --trees OUT/out.trees --truth OUT/truth.npz \
      --mu 5e-9 --recomb 1e-8 --Ne 2000 --n_per_site 8 --out OUT/genetics.npz
"""
import argparse, json, os
import numpy as np
import tskit, msprime, pyslim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--mu", type=float, default=5e-9)
    ap.add_argument("--recomb", type=float, default=1e-8)
    ap.add_argument("--Ne", type=float, default=2000)
    ap.add_argument("--n_per_site", type=int, default=8)
    ap.add_argument("--max_snps", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    ts = tskit.load(a.trees)
    ts = pyslim.recapitate(ts, recombination_rate=a.recomb,
                           ancestral_Ne=a.Ne, random_seed=a.seed)
    ts = msprime.sim_mutations(ts, rate=a.mu, random_seed=a.seed, keep=False)
    print("after mutation: sites=%d" % ts.num_sites)

    truth = np.load(a.truth, allow_pickle=True)
    focal_xy = truth["focal_xy"]
    m = focal_xy.shape[0]

    alive = pyslim.individuals_alive_at(ts, 0)
    locs = ts.individuals_location[alive][:, :2]
    if len(alive) < m * 2:
        raise SystemExit("too few alive individuals (%d) for %d sites" % (len(alive), m))

    chosen, site_of, used = [], [], set()
    for k in range(m):
        order = np.argsort(np.sum((locs - focal_xy[k]) ** 2, axis=1))
        cnt = 0
        for idx in order:
            if idx in used:
                continue
            used.add(idx); chosen.append(alive[idx]); site_of.append(k); cnt += 1
            if cnt >= a.n_per_site:
                break
    chosen = np.array(chosen); site_of = np.array(site_of)

    ind_nodes = [ts.individual(i).nodes for i in chosen]
    sample_nodes = np.array([n for nodes in ind_nodes for n in nodes])
    G = ts.genotype_matrix(samples=sample_nodes)
    dos = G[:, 0::2] + G[:, 1::2]
    p = dos.mean(axis=1) / 2.0
    poly = (p > 0) & (p < 1)
    dos = dos[poly]
    if dos.shape[0] > a.max_snps:
        dos = dos[rng.choice(dos.shape[0], a.max_snps, replace=False)]
    n_snp = dos.shape[0]
    print("polymorphic SNPs retained (nu): %d" % n_snp)

    P = np.zeros((m, n_snp))
    for k in range(m):
        P[k] = dos[:, np.where(site_of == k)[0]].mean(axis=1) / 2.0

    Pc = P - P.mean(axis=0, keepdims=True)
    COV = (Pc @ Pc.T) / n_snp

    S = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            S[i, j] = np.mean((P[i] - P[j]) ** 2)
    S = 0.5 * (S + S.T); np.fill_diagonal(S, 0.0)

    np.savez(a.out, COV=COV, S=S, focal_xy=focal_xy, nu=n_snp,
             theta=truth["theta"], names=truth["names"],
             nx=truth["nx"], ny=truth["ny"], conductance=truth["conductance"],
             X=truth["X"], focal_cell=truth["focal_cell"])

    d = os.path.dirname(a.out) or "."
    names = list(map(str, truth["names"]))
    np.savetxt(os.path.join(d, "S.csv"), S, delimiter=",")
    np.savetxt(os.path.join(d, "COV.csv"), COV, delimiter=",")
    np.savetxt(os.path.join(d, "X.csv"), np.asarray(truth["X"]),
               delimiter=",", header=",".join(names), comments="")
    fc = np.asarray(truth["focal_cell"])                 # (m,2) 0-based fx,fy
    np.savetxt(os.path.join(d, "focal_cell.csv"),
               np.column_stack([fc[:, 0] + 1, fc[:, 1] + 1]),
               delimiter=",", header="col,row", comments="", fmt="%d")
    with open(os.path.join(d, "dims.csv"), "w") as fh:
        fh.write("nx,ny,nu\n%d,%d,%d\n" % (int(truth["nx"]), int(truth["ny"]), n_snp))
    np.savetxt(os.path.join(d, "theta_true.csv"), np.asarray(truth["theta"])[None, :],
               delimiter=",", header=",".join(names), comments="")

    print(json.dumps(dict(n_sites=int(ts.num_sites), nu=int(n_snp), m=int(m),
                          n_alive=int(len(alive)),
                          S_offdiag_mean=float(S[np.triu_indices(m, 1)].mean())),
                     indent=2))


if __name__ == "__main__":
    main()
