"""
Recapitate a SLiM tree sequence, overlay neutral mutations, sample individuals
(by subpopulation or by nearest focal point), and write per-deme derived-allele
counts Y (demes x loci) and haploid sample sizes N for terradish::cov_from_biallelic.

Usage (driven from run_*_recap.R):
  python recap_sample.py --trees x.trees --out x \
      --ancestral_Ne 5000 --recombination_rate 1e-8 --mutation_rate 1e-8 \
      --nsample 25 --seed 1 --mode population
  ... --mode spatial --focal focal_xy.csv   (focal CSV: header + columns fx,fy)
"""
import argparse
import warnings
import numpy as np
import tskit
import pyslim
import msprime

# 1 SLiM tick == 1 generation for the WF / discrete-generation models here.
warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ancestral_Ne", type=float, required=True)
    ap.add_argument("--recombination_rate", type=float, default=1e-8)
    ap.add_argument("--mutation_rate", type=float, default=1e-8)
    ap.add_argument("--nsample", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", choices=["population", "spatial"], default="population")
    ap.add_argument("--focal", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    ts = tskit.load(args.trees)
    ts = pyslim.recapitate(ts, ancestral_Ne=args.ancestral_Ne,
                           recombination_rate=args.recombination_rate,
                           random_seed=args.seed)
    ts = msprime.sim_mutations(ts, rate=args.mutation_rate,
                               model=msprime.BinaryMutationModel(),
                               random_seed=args.seed, keep=False)

    alive = np.array(pyslim.individuals_alive_at(ts, 0))
    pops = np.array([ts.node(ts.individual(int(i)).nodes[0]).population for i in alive])
    locs = np.array([ts.individual(int(i)).location[:2] for i in alive])

    if args.mode == "population":
        groups_all = pops
        group_ids = np.array(sorted(np.unique(pops)))
    else:
        focal = np.loadtxt(args.focal, delimiter=",", skiprows=1)
        fxy = np.atleast_2d(focal)[:, :2]
        d2 = ((locs[:, None, :] - fxy[None, :, :]) ** 2).sum(2)
        groups_all = d2.argmin(1)
        group_ids = np.arange(fxy.shape[0])

    sel_nodes = []
    sel_node_group = []
    nsamp_deme = {}
    for g in group_ids:
        idx = alive[groups_all == g]
        if len(idx) == 0:
            nsamp_deme[int(g)] = 0
            continue
        take = min(args.nsample, len(idx))
        chosen = rng.choice(idx, size=take, replace=False)
        nsamp_deme[int(g)] = int(take)
        for ind in chosen:
            for nd in ts.individual(int(ind)).nodes:
                sel_nodes.append(int(nd))
                sel_node_group.append(int(g))

    sel_nodes = np.array(sel_nodes)
    sel_node_group = np.array(sel_node_group)

    tss = ts.simplify(sel_nodes, keep_input_roots=True)
    G = tss.genotype_matrix()  # sites x nodes, columns in sel_nodes order

    demes_order = np.array(sorted(np.unique(sel_node_group)))
    n_sites = G.shape[0]
    Y = np.zeros((len(demes_order), n_sites), dtype=np.int64)
    N = np.zeros(len(demes_order), dtype=np.int64)
    for di, g in enumerate(demes_order):
        cols = sel_node_group == g
        Y[di, :] = G[:, cols].sum(1)
        N[di] = int(cols.sum())

    np.savetxt(args.out + "_Y.csv", Y, delimiter=",", fmt="%d")
    np.savetxt(args.out + "_N.csv", N, delimiter=",", fmt="%d")
    np.savetxt(args.out + "_demeids.csv", demes_order, delimiter=",", fmt="%d")
    print(f"recap+mut done: {len(demes_order)} demes, {n_sites} sites, "
          f"{len(sel_nodes)} nodes; nsample/deme min={N.min()} max={N.max()}")


if __name__ == "__main__":
    main()
