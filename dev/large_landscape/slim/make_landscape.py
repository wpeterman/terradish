#!/usr/bin/env python3
"""
make_landscape.py — build a known conductance surface and its SLiM habitat map.

Produces, for a scenario, the GROUND TRUTH that terradish must recover:
  * covariate raster(s) X         (what terradish regresses conductance on)
  * true coefficients theta_true  (the answer)
  * conductance c = exp(X theta)  (per cell)
  * habitat.png                   (grayscale map SLiM reads via Image().floatK;
                                   habitat in [hab_min,1], monotone in conductance)
  * focal coordinates (in SLiM's unit [0,1]x[0,1] space)

Habitat <-> conductance link: SLiM gene flow is impeded where carrying capacity
is low and where dispersal must cross low-occupancy cells, so we set habitat as a
monotone increasing function of conductance. This makes high-conductance cells
both more occupied and easier to traverse — the resistance-distance assumption.

Usage:
  python3 make_landscape.py --scenario barrier --nx 64 --ny 64 --out OUTDIR
Scenarios: smooth | barrier | corridor
"""
import argparse, json, os
import numpy as np

def build_covariates(scenario, nx, ny, rng):
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    xn = xx / (nx - 1); yn = yy / (ny - 1)          # normalised 0..1
    # smooth background gradient covariate, standardised
    grad = (xn - 0.5) + 0.4 * (yn - 0.5)
    grad = (grad - grad.mean()) / grad.std()
    layers = {"gradient": grad}
    theta = {"gradient": 0.8}                        # gradient raises conductance
    if scenario == "smooth":
        pass
    elif scenario == "barrier":
        # a near-vertical low-conductance barrier with a narrow gap (corridor)
        bar = np.zeros((ny, nx))
        col = nx // 2
        bar[:, col-0:col+1] = 1.0
        gap = slice(int(0.45*ny), int(0.55*ny))      # one passable gap
        bar[gap, col-0:col+1] = 0.0
        layers["barrier"] = bar
        theta["barrier"] = -2.5                       # strong resistance
    elif scenario == "corridor":
        # a high-conductance horizontal corridor through a low background
        cor = np.zeros((ny, nx))
        row = slice(int(0.45*ny), int(0.55*ny))
        cor[row, :] = 1.0
        layers["corridor"] = cor
        theta["corridor"] = 2.0
        theta["gradient"] = 0.3
    else:
        raise ValueError(scenario)
    return layers, theta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="barrier")
    ap.add_argument("--nx", type=int, default=64)
    ap.add_argument("--ny", type=int, default=64)
    ap.add_argument("--n_focal", type=int, default=25)
    ap.add_argument("--hab_min", type=float, default=0.08)  # floor so barriers stay sparsely occupied, not empty
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    layers, theta = build_covariates(a.scenario, a.nx, a.ny, rng)
    names = list(layers.keys())
    X = np.stack([layers[n].ravel() for n in names], axis=1)   # (N, p)
    th = np.array([theta[n] for n in names])
    logc = X @ th
    logc -= logc.mean()                                       # identifiability: centre log-conductance
    c = np.exp(logc).reshape(a.ny, a.nx)

    # habitat map: monotone in conductance, scaled to [hab_min, 1]
    cn = (c - c.min()) / (c.max() - c.min() + 1e-12)
    habitat = a.hab_min + (1 - a.hab_min) * cn

    # write grayscale PNG (row 0 = top). SLiM Image().floatK returns [0,1] with
    # image row 0 at the TOP; we store so that habitat[y,x] maps to SLiM (x,y).
    try:
        from PIL import Image
        img = (np.clip(habitat, 0, 1) * 255).astype(np.uint8)
        Image.fromarray(img, mode="L").save(os.path.join(a.out, "habitat.png"))
        png_ok = True
    except Exception:
        png_ok = False
    # always write the matrix as CSV fallback (SLiM can readCSV -> matrix)
    np.savetxt(os.path.join(a.out, "habitat.csv"), habitat, delimiter=",", fmt="%.6f")

    # focal coordinates in SLiM unit space, drawn preferentially from occupied cells
    probs = (habitat / habitat.sum()).ravel()
    idx = rng.choice(a.nx * a.ny, size=a.n_focal, replace=False, p=probs)
    fy, fx = np.divmod(idx, a.nx)
    # cell centre -> unit coords; SLiM y increases upward, image row 0 at top
    focal_xy = np.column_stack([(fx + 0.5) / a.nx, 1.0 - (fy + 0.5) / a.ny])

    np.savez(os.path.join(a.out, "truth.npz"),
             X=X, theta=th, names=np.array(names), conductance=c,
             habitat=habitat, nx=a.nx, ny=a.ny,
             focal_xy=focal_xy, focal_cell=np.column_stack([fx, fy]))
    meta = dict(scenario=a.scenario, nx=a.nx, ny=a.ny, names=names,
                theta_true={n: float(v) for n, v in theta.items()},
                n_focal=a.n_focal, hab_min=a.hab_min, png=png_ok)
    with open(os.path.join(a.out, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))
    print("conductance range: %.3f .. %.3f" % (c.min(), c.max()))

if __name__ == "__main__":
    main()
