import argparse
import csv
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pymatching

from simulate import make_circuit


def run_point(args):
    d, p, max_shots, max_errors, batch = args
    circuit = make_circuit(d, p)
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()
    shots = errors = 0
    while shots < max_shots and errors < max_errors:
        det, obs = sampler.sample(batch, separate_observables=True)
        pred = matching.decode_batch(det)
        errors += int(np.sum(np.any(pred != obs, axis=1)))
        shots += batch
    return d, p, shots, errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", type=int, nargs="+", default=[3, 5, 7, 9])
    ap.add_argument("--p-min", type=float, default=0.001)
    ap.add_argument("--p-max", type=float, default=0.012)
    ap.add_argument("--n-p", type=int, default=12)
    ap.add_argument("--max-shots", type=int, default=200_000)
    ap.add_argument("--max-errors", type=int, default=1000)
    ap.add_argument("--batch", type=int, default=10_000)
    ap.add_argument("--out", default="data/mwpm_sweep.csv")
    a = ap.parse_args()

    ps = np.geomspace(a.p_min, a.p_max, a.n_p)
    jobs = [(d, float(p), a.max_shots, a.max_errors, a.batch) for d in a.distances for p in ps]
    t0 = time.time()
    with ProcessPoolExecutor() as ex:
        results = list(ex.map(run_point, jobs))
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["distance", "p", "shots", "errors"])
        w.writerows(results)
    print(f"Wrote {len(results)} points to {a.out} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
