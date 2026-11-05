import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pymatching
import pytest

from simulate import make_circuit
from union_find import UnionFindDecoder


@pytest.mark.parametrize("d", [3, 5])
def test_single_fault_is_corrected(d):
    """Every single error mechanism must be decoded correctly (code distance >= 3)."""
    circuit = make_circuit(d, 0.001)
    matching = pymatching.Matching.from_detector_error_model(
        circuit.detector_error_model(decompose_errors=True))
    uf = UnionFindDecoder(matching)
    for u, v, attr in matching.edges():
        defects = [x for x in (u, v) if x is not None]
        expected = int(0 in attr["fault_ids"])
        assert uf.decode_mask(defects) == expected


def test_no_defects_no_correction():
    uf = UnionFindDecoder.from_circuit(make_circuit(3, 0.001))
    assert uf.decode_mask([]) == 0


def test_agrees_with_mwpm_at_low_noise():
    d, shots = 5, 5000
    circuit = make_circuit(d, 0.001)
    matching = pymatching.Matching.from_detector_error_model(
        circuit.detector_error_model(decompose_errors=True))
    uf = UnionFindDecoder(matching)
    det, obs = circuit.compile_detector_sampler(seed=1).sample(shots, separate_observables=True)
    mw = matching.decode_batch(det)
    u = uf.decode_batch(det)
    assert np.mean(np.any(mw != u, axis=1)) < 0.01


def test_logical_error_rate_below_threshold_improves_with_distance():
    rates = []
    for d in (3, 5):
        circuit = make_circuit(d, 0.003)
        uf = UnionFindDecoder.from_circuit(circuit)
        det, obs = circuit.compile_detector_sampler(seed=2).sample(20000, separate_observables=True)
        rates.append(np.mean(np.any(uf.decode_batch(det) != obs, axis=1)))
    assert rates[1] < rates[0]
