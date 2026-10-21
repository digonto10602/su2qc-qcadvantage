"""M4: noisy-emulation outputs (`results/m4/noise_runs.json`, `reports/m4_noise.md`).

API used: `su2torus/noise.py` (specified in test_m4_noise_api.py).

Schema of results/m4/noise_runs.json (produced by scripts/run_m4.py):
{
  "N": 12, "L1": 3, "L2": 2, "g": 1.25, "dt": 0.703125 (= 0.45 g^2), "order": 1,
  "state": "stripe", "p1": 3e-5, "readout": 1e-3,
  "fidelity_formula": "<plain-text statement of the predicted-F formula and eps convention>",
  "mitigation_method": "rescaled" | "zne",      # ONE method, used for every configuration
  "configs": [                                  # must contain every (steps, p2) in CONFIGS
    {"steps": int, "p2": float, "shots": int | null (null = exact density matrix), "seed": int,
     "n_2q": int, "n_1q": int,                  # gate counts of the evolution circuit (no prep)
     "ideal_electric_energy": float,            # noiseless circuit <H_E>
     "raw_electric_energy": float,              # noisy <H_E>, no mitigation
     "F_predicted": float,                      # predicted_fidelity(evolution circuit, p2, p1)
     "F_echo": float,                           # echo_fidelity(...)
     "F_rescale": float,                        # factor used in the rescaling: F_echo, or an
                                                #   observable-specific echo decay factor
                                                #   (the report must say which)
     "rescaled_electric_energy": float,         # mitigate_rescale(raw, F_rescale, O_inf)
     "zne_values": [v1, v3, v5],                # noisy <H_E> at fold factors 1, 3, 5
     "zne_electric_energy": float},             # zne(zne_values, (1, 3, 5))
    ...],
  "repro": {"steps": 2, "p2": float, "shots": int <= 4000, "seed": int,
            "raw_electric_energy": float}       # = run_noisy(trotter_circuit(lat, g, dt, 2, 1,
                                                #   False, stripe), noise_model(p2), shots, seed)
                                                #   @ diag(H_E); recomputed live below
}
<H_E> is the full electric energy of METHODS sec. 2 (constant included); O_inf = 5.2734375.
The "ideal" reference is the noiseless Trotter circuit (mitigation targets noise, not Trotter error).
Physics fixed by this test: N = 12 torus, g = 1.25, stripe state, first-order Trotter,
dt = 0.45 g^2 (2 dt/g^2 = 0.9, the frontier-run choice of PLAN M6).
Ideal values are recomputed here from the brute-force reference, not read from the code.
"""
import functools
import json
import pathlib
import numpy as np
import pytest
from su2torus.lattice import HoneycombTorus
from su2torus.trotter import trotter_circuit
from scipy.sparse.linalg import expm_multiply
from su2torus.noise import noise_model, run_noisy
from tests.acceptance._ref import ref_neighbors, ref_hamiltonian, basis_vec

G = 1.25
DT = 0.45 * G ** 2          # 2 dt / g^2 = 0.9, as in the frontier run (PLAN M6)


def stripe(L1, L2):
    return [2 * (i * L2 + j) for i in range(0, L1, 2) for j in range(L2)]


@functools.lru_cache(maxsize=None)
def ref_parts(L1, L2, g):
    return ref_hamiltonian(ref_neighbors(L1, L2, True), g, parts=True)


def ref_trotter_state(L1, L2, g, dt, steps, init):
    """Order-1 Trotter state from brute-force H parts (METHODS section 3); U_E acts first."""
    HE, HA, HB = ref_parts(L1, L2, g)
    dE = HE.diagonal()
    v = basis_vec(2 * L1 * L2, init)
    for _ in range(steps):
        v = np.exp(-1j * dE * dt) * v
        v = expm_multiply(-1j * dt * HA, v)    # commuting terms: exact sublattice exponential
        v = expm_multiply(-1j * dt * HB, v)
    return v, dE


def n_gates(circ):
    ops = circ.count_ops()
    n2 = sum(v for k, v in ops.items() if k in ("cx", "rzz"))
    n1 = sum(v for k, v in ops.items() if k not in ("cx", "rzz", "barrier"))
    return n2, n1

ROOT = pathlib.Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "m4" / "noise_runs.json"
L1, L2 = 3, 2
P1, READOUT = 3e-5, 1e-3
O_INF = 0.5 * G ** 2 * 0.375 * 18             # METHODS sec. 4, 18 bonds on the N=12 torus
CONFIGS = [(2, 7.9e-4), (4, 1.0e-3), (6, 1.5e-3), (10, 1.0e-3)]


@pytest.fixture(scope="module")
def data():
    assert RES.exists(), "run scripts/run_m4.py first"
    d = json.loads(RES.read_text())
    assert (d["N"], d["L1"], d["L2"], d["order"], d["state"]) == (12, L1, L2, 1, "stripe")
    assert abs(d["g"] - G) < 1e-12 and abs(d["dt"] - DT) < 1e-12
    assert abs(d["p1"] - P1) < 1e-15 and abs(d["readout"] - READOUT) < 1e-15
    return d


def find(d, steps, p2):
    hits = [c for c in d["configs"] if c["steps"] == steps and abs(c["p2"] - p2) < 1e-12]
    assert len(hits) == 1, f"missing or duplicate config {(steps, p2)}"
    return hits[0]


@pytest.mark.parametrize("steps,p2", CONFIGS)
def test_config_consistency(data, steps, p2):
    c = find(data, steps, p2)
    assert isinstance(c["seed"], int)
    psi, dE = ref_trotter_state(L1, L2, G, DT, steps, stripe(L1, L2))
    ideal = float(np.abs(psi) ** 2 @ dE)
    # 1e-8: exact statevector values; leaves room for the code's own evolution route
    assert abs(c["ideal_electric_energy"] - ideal) < 1e-8
    circ = trotter_circuit(HoneycombTorus(L1, L2), G, DT, steps, 1, False, None)
    n2, n1 = n_gates(circ)
    assert (c["n_2q"], c["n_1q"]) == (n2, n1)
    Fp = (1 - 15 * p2 / 16) ** n2 * (1 - 3 * P1 / 4) ** n1
    assert abs(c["F_predicted"] / Fp - 1) < 1e-9          # same closed form, round-off only
    assert 0 < c["F_rescale"] <= 1
    # 1e-9: the stored mitigated values must be the stated formulas applied to stored inputs
    resc = O_INF + (c["raw_electric_energy"] - O_INF) / c["F_rescale"]
    assert abs(c["rescaled_electric_energy"] - resc) < 1e-9
    v = c["zne_values"]
    assert len(v) == 3 and abs(v[0] - c["raw_electric_energy"]) < 1e-12
    assert abs(c["zne_electric_energy"] - (15 / 8 * v[0] - 5 / 4 * v[1] + 3 / 8 * v[2])) < 1e-9


@pytest.mark.parametrize("steps,p2", CONFIGS)
def test_echo_matches_prediction(data, steps, p2):
    c = find(data, steps, p2)
    # 25%: PLAN M4 acceptance criterion (relative to the predicted F)
    assert abs(c["F_echo"] - c["F_predicted"]) <= 0.25 * c["F_predicted"], c


def test_mitigation_halves_error(data):
    key = {"rescaled": "rescaled_electric_energy", "zne": "zne_electric_energy"}[
        data["mitigation_method"]]
    wins, table = 0, []
    for steps, p2 in CONFIGS:
        c = find(data, steps, p2)
        psi, dE = ref_trotter_state(L1, L2, G, DT, steps, stripe(L1, L2))
        ideal = float(np.abs(psi) ** 2 @ dE)
        raw_err = abs(c["raw_electric_energy"] - ideal)
        mit_err = abs(c[key] - ideal)
        table.append((steps, p2, raw_err, mit_err))
        wins += mit_err <= 0.5 * raw_err                  # 0.5: PLAN M4 acceptance criterion
    assert wins >= 3, f"only {wins}/4 configs halve the error: {table}"


def test_reproducible(data):
    r = data["repro"]
    assert r["steps"] == 2 and r["shots"] <= 4000 and isinstance(r["seed"], int)
    lat = HoneycombTorus(L1, L2)
    circ = trotter_circuit(lat, G, DT, r["steps"], 1, False, stripe(L1, L2))
    _, dE = ref_trotter_state(L1, L2, G, DT, 0, [])
    p = np.asarray(run_noisy(circ, noise_model(r["p2"], P1, READOUT), shots=r["shots"],
                             seed=r["seed"]))
    # 1e-12: same seed must reproduce the sampled frequencies exactly (round-off in the sum only)
    assert abs(float(p @ dE) - r["raw_electric_energy"]) < 1e-12


def test_report():
    rep = ROOT / "reports" / "m4_noise.md"
    assert rep.exists()
    text = rep.read_text()
    assert "\\prod" in text, "the report must state the predicted-fidelity formula"
