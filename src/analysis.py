import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load(path):
    df = pd.read_csv(path)
    df["pL"] = df.errors / df.shots
    # binomial standard error (floor at 1 error to avoid zero-width bars)
    df["err"] = np.sqrt(np.maximum(df.errors, 1)) / df.shots
    return df.sort_values(["distance", "p"])


def crossing(df, d1, d2):
    """Crossing of the log-log curves of two distances, via linear interpolation."""
    a, b = df[df.distance == d1], df[df.distance == d2]
    common = set(a[a.errors > 0].p) & set(b[b.errors > 0].p)  # drop zero-error points
    a, b = a[a.p.isin(common)], b[b.p.isin(common)]
    p = a.p.values
    diff = np.log(b.pL.values) - np.log(a.pL.values)  # <0 below threshold, >0 above
    for i in range(len(p) - 1):
        if diff[i] < 0 <= diff[i + 1]:
            t = -diff[i] / (diff[i + 1] - diff[i])
            return float(np.exp(np.log(p[i]) + t * (np.log(p[i + 1]) - np.log(p[i]))))
    return None


def threshold_plot(df, path):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for d, g in df.groupby("distance"):
        ax.errorbar(g.p, g.pL, yerr=g.err, marker="o", ms=4, capsize=2, label=f"d = {d}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("physical error rate p"); ax.set_ylabel("logical error rate (d rounds)")
    ax.set_title("Rotated surface code, MWPM, circuit-level noise")
    ax.legend(); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def subthreshold(df, p_th, path):
    """Fit log pL = c - s*d at fixed p; Lambda = exp(2s) so pL ~ Lambda^-(d+1)/2."""
    rows = []
    for p, g in df.groupby("p"):
        if p >= 0.7 * p_th or (g.errors < 20).any() or g.distance.nunique() < 3:
            continue
        s, c = np.polyfit(g.distance, -np.log(g.pL), 1)
        rows.append((p, np.exp(2 * s)))
    if not rows:
        print("Not enough sub-threshold data with >=20 errors for a fit; increase --max-shots.")
        return None
    ps, lam = map(np.array, zip(*rows))
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(ps, lam, "o", label="measured $\\Lambda(p)$")
    ax.loglog(ps, p_th / ps, "--", label="$p_{th}/p$ (ideal scaling)")
    ax.set_xlabel("physical error rate p"); ax.set_ylabel("error suppression factor $\\Lambda$")
    ax.legend(); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(path, dpi=200); plt.close(fig)
    return ps, lam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/mwpm_sweep.csv")
    ap.add_argument("--fig-dir", default="figures")
    a = ap.parse_args()
    df = load(a.data)
    threshold_plot(df, f"{a.fig_dir}/threshold_mwpm.png")
    ds = sorted(df.distance.unique())
    ests = []
    for d1, d2 in zip(ds[:-1], ds[1:]):
        c = crossing(df, d1, d2)
        print(f"crossing d={d1} vs d={d2}: {c}")
        if c: ests.append(c)
    if ests:
        p_th = float(np.mean(ests))
        print(f"Threshold estimate: {p_th:.4f} (spread of pairwise crossings: +/- {np.std(ests):.4f})")
        out = subthreshold(df, p_th, f"{a.fig_dir}/subthreshold_lambda.png")
        if out is not None:
            for p, l in zip(*out):
                print(f"p={p:.4f}  Lambda={l:.2f}  (p_th/p={p_th / p:.2f})")


if __name__ == "__main__":
    main()
