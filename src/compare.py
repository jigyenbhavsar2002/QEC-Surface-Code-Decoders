import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymatching

from analysis import crossing
from simulate import make_circuit
from union_find import UnionFindDecoder


def run(distances, ps, shots):
    rows = []
    for d in distances:
        for p in ps:
            circuit = make_circuit(d, p)
            dem = circuit.detector_error_model(decompose_errors=True)
            mw = pymatching.Matching.from_detector_error_model(dem)
            uf = UnionFindDecoder(mw)
            det, obs = circuit.compile_detector_sampler().sample(shots, separate_observables=True)
            t = time.perf_counter(); pm = mw.decode_batch(det); t_mw = time.perf_counter() - t
            t = time.perf_counter(); pu = uf.decode_batch(det); t_uf = time.perf_counter() - t
            rows.append(dict(distance=d, p=float(p), shots=shots,
                             err_mwpm=int(np.sum(np.any(pm != obs, axis=1))),
                             err_uf=int(np.sum(np.any(pu != obs, axis=1))),
                             t_mwpm=t_mw, t_uf=t_uf, mean_defects=float(det.sum(1).mean())))
            print(rows[-1], flush=True)
    return pd.DataFrame(rows)


def plots(df, fig_dir):
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ths = {}
    for name, col, ls in (("MWPM", "err_mwpm", "-"), ("Union-find", "err_uf", "--")):
        sub = df.assign(errors=df[col], pL=df[col] / df.shots)
        for d, g in sub.groupby("distance"):
            g = g[g.errors > 0]  # zero-error points can't be drawn on a log axis
            ax.errorbar(g.p, g.pL, yerr=np.sqrt(g.errors) / g.shots, fmt=ls, marker="o", ms=3,
                        capsize=2, label=f"{name} d={d}",
                        color=f"C{sorted(df.distance.unique()).index(d)}")
        ds = sorted(df.distance.unique())
        cs = [c for c in (crossing(sub.sort_values(["distance", "p"]), a, b)
                          for a, b in zip(ds[:-1], ds[1:])) if c]
        ths[name] = (float(np.mean(cs)), float(np.std(cs))) if cs else (None, None)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("physical error rate p"); ax.set_ylabel("logical error rate")
    ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3, which="both")
    ax.set_title("MWPM (solid) vs union-find (dashed)")
    fig.tight_layout(); fig.savefig(f"{fig_dir}/mwpm_vs_uf.png", dpi=200); plt.close(fig)

    # runtime vs distance at a fixed p near the middle of the sweep
    p0 = df.p.iloc[(df.p - 0.004).abs().argmin()]
    g = df[df.p == p0]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(g.distance, 1e6 * g.t_uf / g.shots, "o--", label="Union-find (pure Python)")
    ax.plot(g.distance, 1e6 * g.t_mwpm / g.shots, "s-", label="MWPM (PyMatching, C++)")
    ax.set_yscale("log"); ax.set_xlabel("distance d"); ax.set_ylabel("decode time per shot (µs)")
    ax.set_xticks(sorted(df.distance.unique())); ax.set_title(f"Runtime at p = {p0:.4f}"); ax.legend(); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(f"{fig_dir}/runtime.png", dpi=200); plt.close(fig)
    return ths


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", type=int, nargs="+", default=[3, 5, 7])
    ap.add_argument("--p-min", type=float, default=0.001)
    ap.add_argument("--p-max", type=float, default=0.012)
    ap.add_argument("--n-p", type=int, default=10)
    ap.add_argument("--shots", type=int, default=20_000)
    ap.add_argument("--out", default="data/compare.csv")
    ap.add_argument("--fig-dir", default="figures")
    a = ap.parse_args()
    df = run(a.distances, np.geomspace(a.p_min, a.p_max, a.n_p), a.shots)
    df.to_csv(a.out, index=False)
    for k, (m, s) in plots(df, a.fig_dir).items():
        print(f"{k}: threshold ~ {m:.4f} +/- {s:.4f}" if m else f"{k}: no crossing found")
