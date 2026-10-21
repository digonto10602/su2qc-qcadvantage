"""M5 (a): MPS simulation of the Trotter circuits -- API-level acceptance tests (live, N = 12).

New public API: module `su2torus/mps.py`  (quimb CircuitMPS; docs/METHODS.md sections 2-4)
------------------------------------------------------------------------------------------
snake_order(lat) -> list[int]
    Permutation of range(lat.n); entry k is the plaquette held by MPS site k. Boustrophedon
    ("snake") over the cells, lines running along the SHORTER side of the torus (smaller MPS cut;
    reviewer check at N=24, dt=0.2: Schmidt tail beyond 256 is 8e-14 for lines along the short
    side vs 2e-10 along the long side), the two plaquettes of a cell on consecutive sites,
    direction reversed on alternate lines. Example for L1 <= L2 (lines along i):
        for j in range(L2): for i in (range(L1) if j even else reversed(range(L1))):
            [up(i,j), down(i,j)]  ([down(i,j), up(i,j)] on odd j)
    and the same with the roles of i and j exchanged when L1 > L2.
    Tested property: at most max(L1, L2) - 1 consecutive site pairs are NOT honeycomb neighbours.

simulate_mps(lat, g, dt, steps, order=1, abelian=False, initial=None, chi=256, cutoff=1e-12,
             backend="numpy") -> dict
    Applies, gate by gate, exactly the circuit su2torus.trotter.trotter_circuit(lat, g, dt, steps,
    order, abelian, initial) to a quimb CircuitMPS (max_bond=chi, cutoff=cutoff) whose site k is
    plaquette snake_order(lat)[k]. Measures observables after every Trotter step.
    backend: "numpy" (CPU) or a GPU array backend ("cupy" / "torch"); any other string ->
    ValueError; a GPU backend whose library is missing -> ImportError or RuntimeError
    (never a silent fallback to CPU).
    Returned dict (JSON-serialisable except "mps"):
      "n"               : int
      "chi"             : int (the max_bond used)
      "site_order"      : list[int] = snake_order(lat)
      "times"           : [k * dt for k = 0..steps]
      "electric_energy" : list (steps+1) of <H_E>, H_E = (3 g^2/16) sum_bonds (1 - Z_p Z_q)
      "z"               : list (steps+1) of list n, <Z_p> in PLAQUETTE index order (not site order)
      "fidelity_estimate"       : float, quimb CircuitMPS.fidelity_estimate() at the end
                                  (product of the kept norms of all truncations)
      "error_estimate"          : float = 1 - fidelity_estimate, in [0, 1]  (the discarded-weight
                                  error estimate; quimb CircuitMPS.error_estimate())
      "error_estimate_per_step" : list (steps+1), error_estimate after each step (0 at k = 0)
      "max_bond"        : int, largest bond dimension reached (<= chi)
      "wall_time_s"     : float
      "mps"             : the final quimb MatrixProductState (sites in snake order)

mps_statevector(res) -> np.ndarray, complex, shape (2**n,)
    Dense, normalised final state of a simulate_mps result in the METHODS basis order
    b = sum_p n_p 2**p (i.e. the snake permutation undone). Needed only for n <= 24.

Physics point used throughout: g = 1.25, stripe initial state (METHODS section 4).
"""
import inspect

import numpy as np
import pytest

from tests.acceptance._ref import ref_neighbors
from tests.acceptance.m5._exact_m5 import exact_trotter, stripe, z_and_he

G = 1.25
DT = 0.2           # 4 steps reach t = 0.8; same dt as results/m5/mps_check.json (chi=256 suffices at N=24)
STEPS = 4
L1, L2 = 3, 2      # N = 12 torus (METHODS section 4)


def _lat(a, b):
    from su2torus.lattice import HoneycombTorus
    return HoneycombTorus(a, b)


def _fid(a, b):
    return abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real)


@pytest.fixture(scope="module")
def mps():
    import su2torus.mps as m
    return m


@pytest.fixture(scope="module")
def run256(mps):
    return mps.simulate_mps(_lat(L1, L2), G, DT, STEPS, order=1, abelian=False,
                            initial=stripe(L1, L2), chi=256)


@pytest.mark.parametrize("dims", [(3, 2), (4, 3), (4, 6), (4, 7), (4, 8), (6, 8)])
def test_snake_order(mps, dims):
    a, b = dims
    order = list(mps.snake_order(_lat(a, b)))
    n = 2 * a * b
    assert sorted(order) == list(range(n))
    nb = ref_neighbors(a, b, True)
    jumps = sum(1 for k in range(n - 1) if order[k + 1] not in nb[order[k]])
    assert jumps <= max(a, b) - 1, f"snake order has {jumps} non-adjacent steps"


def test_signature_and_backend_flag(mps):
    params = inspect.signature(mps.simulate_mps).parameters
    for name in ("lat", "g", "dt", "steps", "order", "abelian", "initial", "chi", "cutoff",
                 "backend"):
        assert name in params, name
    with pytest.raises(ValueError):
        mps.simulate_mps(_lat(L1, L2), G, DT, 1, initial=stripe(L1, L2), chi=8,
                         backend="not-a-backend")


@pytest.mark.parametrize("order,abelian", [(1, False), (2, False), (1, True)])
def test_mps_equals_exact_statevector_n12(mps, run256, order, abelian):
    """PLAN M5: state fidelity > 1 - 1e-6 at chi = 256, 4 steps (threshold taken from PLAN)."""
    ex = stripe(L1, L2)
    res = run256 if (order, abelian) == (1, False) else mps.simulate_mps(
        _lat(L1, L2), G, DT, STEPS, order=order, abelian=abelian, initial=ex, chi=256)
    psi = np.asarray(mps.mps_statevector(res))
    assert psi.shape == (2 ** 12,)
    ref = exact_trotter(L1, L2, G, DT, STEPS, ex, order=order, abelian=abelian)
    assert _fid(ref, psi) > 1 - 1e-6


def test_observables_per_step_n12(run256):
    res = run256
    ex = stripe(L1, L2)
    ref = exact_trotter(L1, L2, G, DT, STEPS, ex, observe=lambda v: z_and_he(v, L1, L2, G))
    assert res["n"] == 12 and res["chi"] == 256
    assert list(res["site_order"]) and sorted(res["site_order"]) == list(range(12))
    np.testing.assert_allclose(res["times"], [k * DT for k in range(STEPS + 1)], atol=1e-12)
    assert len(res["electric_energy"]) == STEPS + 1 and len(res["z"]) == STEPS + 1
    # chi = 256 exceeds the maximal Schmidt rank 2**6 = 64 at n = 12, so only the 1e-12 cutoff
    # and round-off remain; 1e-6 is a wide margin.
    for k, (z_ref, he_ref) in enumerate(ref):
        np.testing.assert_allclose(res["z"][k], z_ref, atol=1e-6, err_msg=f"step {k}")
        assert abs(res["electric_energy"][k] - he_ref) < 1e-6, k
    assert res["error_estimate"] < 1e-8      # no truncation beyond the cutoff (see above)
    assert res["max_bond"] <= 64


def test_error_estimate_decreases_with_chi(mps):
    """PLAN M5: the discarded-weight error estimate is reported and decreases with chi."""
    ex = stripe(L1, L2)
    ref = exact_trotter(L1, L2, G, DT, STEPS, ex)
    chis = [2, 4, 8, 16, 32, 64]
    err, infid = [], []
    for chi in chis:
        res = mps.simulate_mps(_lat(L1, L2), G, DT, STEPS, initial=ex, chi=chi)
        e = float(res["error_estimate"])
        assert 0.0 <= e <= 1.0
        assert abs(e - (1.0 - float(res["fidelity_estimate"]))) < 1e-12
        assert res["max_bond"] <= chi
        per = np.asarray(res["error_estimate_per_step"], dtype=float)
        assert per.shape == (STEPS + 1,) and per[0] < 1e-12
        # quimb's estimate is 1 - prod(kept norms): cumulative, so non-decreasing in time
        # (1e-12 = round-off allowance)
        assert np.all(np.diff(per) >= -1e-12) and abs(per[-1] - e) < 1e-12
        err.append(e)
        infid.append(1.0 - _fid(ref, np.asarray(mps.mps_statevector(res))))
    # non-increasing in chi (1e-9 = round-off/cutoff allowance), strictly lower at the top
    assert all(err[k + 1] <= err[k] + 1e-9 for k in range(len(chis) - 1)), err
    assert err[0] > 1e-6 and err[0] > 10 * err[-1], err
    assert err[-1] < 1e-8, err          # chi = 64 = 2**6 is exact at n = 12
    # the estimate is meaningful: never reports ~0 while the state is actually truncated
    for e, f in zip(err, infid):
        if f > 1e-6:
            assert e > 1e-9, (err, infid)
