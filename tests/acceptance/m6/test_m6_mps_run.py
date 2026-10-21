"""M6: the N = 48 MPS convergence run (`results/m6/n48_stripe_g1.25.json`).

PLAN M6: "Run MPS on the N=48 torus (4x6), stripe state, g=1.25, 1-8 first-order steps with dt
chosen so that 2 dt/g^2 ~ 0.9, at chi in {32, 64, 128, 256}."

Produced by (no new API; M5 entry point):
    python3 scripts/run_classical.py --config configs/n48_stripe_g1.25.json
Schema (tests/acceptance/m5/test_m5_outputs.py):
    {"config": {"name", "L1": 4, "L2": 6, "g": 1.25, "dt": 0.703125, "steps": 8, "order": 1,
                "abelian": false, "state": "stripe", "chi": [32, 64, 128, 256], "cutoff",
                "backend", "output"},
     "runs": [{"chi": int, "times": [9], "electric_energy": [9], "z": [9 x 48] (plaquette order),
               "error_estimate": float, "error_estimate_per_step": [9], "max_bond": int,
               "wall_time_s": float}, ...]}           runs in the order 32, 64, 128, 256

Independent physics check (no su2torus code): steps 0 and 1 are known in closed form.
  step 0: the stripe product state, <Z_p> = -1 on excited plaquettes, +1 otherwise.
  step 1 (order 1, U = U_B U_A U_E, METHODS sec. 3): U_E only adds phases to a product state;
    U_A rotates every A plaquette independently (its B neighbours are definite);
    U_B rotates every B plaquette b by exp(-i dt a_b X_b) with a_b = -(1/g^2) 0.5^m,
    m = number of excited A neighbours, which are now independent superpositions. Hence
    <Z_b> = z_b(0) sum_cfg P(cfg) cos(2 dt a_b(cfg)) and <Z_a Z_b> likewise (Z_a commutes with
    U_B). This formula was checked by the reviewer against the brute-force Trotter statevector
    of tests/acceptance/m5/_exact_m5.py on the 3x2, 3x3 and 2x4 tori (agreement 1e-13).
"""
import itertools
import json
import math
import pathlib

import numpy as np
import pytest

from tests.acceptance._ref import ref_bonds, ref_neighbors

ROOT = pathlib.Path(__file__).resolve().parents[3]
RUN = ROOT / "results" / "m6" / "n48_stripe_g1.25.json"
L1, L2, N, G, STEPS = 4, 6, 48, 1.25, 8
CHIS = [32, 64, 128, 256]


def stripe(a, b):
    """METHODS sec. 4 'stripe': up(i, j) = 2 (i L2 + j) for even i."""
    return sorted(2 * (i * b + j) for i in range(a) for j in range(b) if i % 2 == 0)


def step1_reference(a, b, g, dt, excited):
    """Closed-form <Z_p> and <H_E> after one first-order Trotter step (see module docstring).
    Valid for any initial product state; sublattice A = even p (METHODS sec. 1)."""
    nb = ref_neighbors(a, b, True)
    n = len(nb)
    exc = set(excited)
    z0 = np.array([-1.0 if p in exc else 1.0 for p in range(n)])
    occ0 = (1 - z0) / 2
    z = z0.copy()
    for p in range(0, n, 2):                       # U_A: controls (B) definite
        amp = -(0.5 ** sum(occ0[q] for q in nb[p])) / g ** 2
        z[p] = z0[p] * math.cos(2 * dt * amp)
    p_exc = (1 - z) / 2                            # A excitation probabilities after U_A
    zz = {}
    for p in range(1, n, 2):                       # U_B: controls (A) independent superpositions
        zb, zab = 0.0, {q: 0.0 for q in nb[p]}
        for cfg in itertools.product((0, 1), repeat=len(nb[p])):
            pr = np.prod([p_exc[q] if c else 1 - p_exc[q] for q, c in zip(nb[p], cfg)])
            v = z0[p] * math.cos(2 * dt * (-(0.5 ** sum(cfg)) / g ** 2))
            zb += pr * v
            for q, c in zip(nb[p], cfg):
                zab[q] += pr * v * (1 - 2 * c)
        z[p] = zb
        for q in nb[p]:
            zz[(min(p, q), max(p, q))] = zab[q]
    bonds = ref_bonds(nb)
    he = (3 * g ** 2 / 16) * sum(1 - zz[bd] for bd in bonds)
    return z, he, len(bonds)


@pytest.fixture(scope="module")
def run():
    assert RUN.exists(), f"missing {RUN}"
    return json.loads(RUN.read_text())


def test_run_config(run):
    """PLAN M6: N=48 (4x6), stripe, g=1.25, first order, 8 steps, 2 dt/g^2 ~ 0.9, chi set."""
    c = run["config"]
    assert (c["L1"], c["L2"]) == (L1, L2)
    assert c["state"] == "stripe" and c["abelian"] is False and c["order"] == 1
    assert abs(c["g"] - G) < 1e-12 and c["steps"] == STEPS
    # "~0.9": 1% band; the committed config has exactly 0.9 (dt = 0.703125)
    assert abs(2 * c["dt"] / G ** 2 - 0.9) <= 0.009, c["dt"]
    assert c["chi"] == CHIS
    assert [r["chi"] for r in run["runs"]] == CHIS


def test_run_records(run):
    """Every chi run covers steps 0..8 with the M5 schema; error estimate recorded per step."""
    dt = run["config"]["dt"]
    for r in run["runs"]:
        chi = r["chi"]
        np.testing.assert_allclose(r["times"], [k * dt for k in range(STEPS + 1)], atol=1e-12)
        assert len(r["electric_energy"]) == STEPS + 1
        assert len(r["z"]) == STEPS + 1 and all(len(zk) == N for zk in r["z"])
        eps = r["error_estimate_per_step"]
        assert len(eps) == STEPS + 1 and abs(eps[0]) <= 1e-12
        assert all(-1e-12 <= e <= 1 + 1e-12 for e in eps), chi
        assert abs(r["error_estimate"] - eps[-1]) <= 1e-12, chi
        assert 1 <= r["max_bond"] <= chi
        assert r["wall_time_s"] > 0
        # physical ranges: |<Z>| <= 1, 0 <= <H_E> <= (3 g^2 / 8) n_bonds (n_bonds = 72); 1e-6 slack
        # for floating-point/normalisation round-off of the MPS contraction
        assert np.all(np.abs(np.asarray(r["z"])) <= 1 + 1e-6), chi
        he_max = 3 * G ** 2 / 8 * 72
        assert all(-1e-6 <= e <= he_max + 1e-6 for e in r["electric_energy"]), chi


def test_step0_exact(run):
    """Step 0 is the stripe product state for every chi (exact; 1e-9 = round-off only)."""
    z_ref = np.array([-1.0 if p in set(stripe(L1, L2)) else 1.0 for p in range(N)])
    _, _, nbonds = step1_reference(L1, L2, G, run["config"]["dt"], stripe(L1, L2))
    he_ref = (3 * G ** 2 / 16) * sum(
        1 - z_ref[p] * z_ref[q] for p, q in ref_bonds(ref_neighbors(L1, L2, True)))
    assert nbonds == 72
    for r in run["runs"]:
        np.testing.assert_allclose(r["z"][0], z_ref, atol=1e-9, err_msg=f"chi={r['chi']}")
        assert abs(r["electric_energy"][0] - he_ref) <= 1e-9, r["chi"]


def test_step1_matches_closed_form(run):
    """Step 1 vs the closed form above, for every chi.
    Tolerance = trace-distance bound with the run's own discarded-weight estimate eps1:
    |<O>_mps - <O>_exact| <= 2 ||O - c|| sqrt(eps1); times 2 because eps1 is an estimate of
    1 - F, not a bound; ||Z|| = 1, ||H_E - c|| = (3 g^2/16) n_bonds. Plus 1e-6 round-off.
    The chi = 256 run must have eps1 <= 1e-3 (Z tolerance <= 0.13), which is loose enough for any
    sane MPS yet still separates the abelian or factor-2-in-dt mistakes (they move <Z_b> by > 0.3)."""
    dt = run["config"]["dt"]
    z_ref, he_ref, nbonds = step1_reference(L1, L2, G, dt, stripe(L1, L2))
    for r in run["runs"]:
        eps1 = max(r["error_estimate_per_step"][1], 0.0)
        tol_z = 1e-6 + 4 * math.sqrt(eps1)
        tol_he = 1e-6 + 4 * (3 * G ** 2 / 16) * nbonds * math.sqrt(eps1)
        np.testing.assert_allclose(r["z"][1], z_ref, atol=tol_z, rtol=0,
                                   err_msg=f"chi={r['chi']} eps1={eps1}")
        assert abs(r["electric_energy"][1] - he_ref) <= tol_he, (r["chi"], eps1)
    top = run["runs"][-1]
    assert top["chi"] == 256 and top["error_estimate_per_step"][1] <= 1e-3
