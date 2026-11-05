import numpy as np
import stim
import pymatching


def make_circuit(d: int, p: float, rounds: int | None = None) -> stim.Circuit:
    """Rotated surface code Z-memory experiment with circuit-level depolarizing noise."""
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d,
        rounds=rounds or d,
        after_clifford_depolarization=p,
        after_reset_flip_probability=p,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p,
    )


def logical_error_rate(circuit: stim.Circuit, shots: int):
    """Sample detection events, decode with MWPM, return (rate, n_errors)."""
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()
    detections, observables = sampler.sample(shots, separate_observables=True)
    predictions = matching.decode_batch(detections)
    n_errors = int(np.sum(np.any(predictions != observables, axis=1)))
    return n_errors / shots, n_errors


if __name__ == "__main__":
    d, p, shots = 3, 0.005, 100_000
    circuit = make_circuit(d, p)
    print(str(circuit.diagram("timeline-text"))[:1500])  # read this until you understand it
    dem = circuit.detector_error_model(decompose_errors=True)
    print("detectors:", dem.num_detectors, "| error mechanisms:", len(dem))
    rate, n = logical_error_rate(circuit, shots)
    print(f"d={d}, p={p}: logical error rate = {rate:.4f} ({n}/{shots})")
