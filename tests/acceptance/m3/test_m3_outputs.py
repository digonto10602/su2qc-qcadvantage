"""M3 / G9-G10: small-torus physics outputs (LOCKED). Schema: docs/METHODS.md section 6."""
import json
import pathlib
import numpy as np
import pytest
from tests.acceptance._ref import ref_neighbors, ref_bonds

ROOT = pathlib.Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "m3"
TORI = {12: (3, 2), 24: (4, 3)}
GRID = [(12, g, s, m) for g in (1.1, 1.25, 1.4) for s in ("vacuum", "neel_half", "stripe")
        for m in ("su2", "abelian")] + \
       [(24, 1.25, s, m) for s in ("vacuum", "neel_half", "stripe") for m in ("su2", "abelian")]


def excited(L1, L2, state):
    if state == "vacuum":
        return []
    if state == "neel_half":
        return [2 * (i * L2 + j) for i in range(L1 // 2) for j in range(L2)]
    return [2 * (i * L2 + j) for i in range(0, L1, 2) for j in range(L2)]


def e_electric_t0(L1, L2, periodic, exc, g):
    bonds = ref_bonds(ref_neighbors(L1, L2, periodic))
    s = set(exc)
    return 0.5 * g ** 2 * sum(0.75 for p, q in bonds if (p in s) != (q in s))


@pytest.fixture(scope="module")
def summary():
    f = RES / "summary.json"
    assert f.exists(), "run scripts/run_m3.py first"
    return json.loads(f.read_text())


def find(summary, N, g, s, m):
    hits = [r for r in summary["runs"] if r["N"] == N and abs(r["g"] - g) < 1e-9
            and r["state"] == s and r["model"] == m]
    assert len(hits) == 1, f"missing or duplicate run {N, g, s, m}"
    return hits[0]


@pytest.mark.parametrize("N,g,s,m", GRID)
def test_run(summary, N, g, s, m):
    r = find(summary, N, g, s, m)
    L1, L2 = TORI[N]
    t = np.array(r["times"])
    EE, EB = np.array(r["electric_energy"]), np.array(r["magnetic_energy"])
    assert t[0] == 0 and t[-1] >= 4 * g ** 2 - 1e-9 and len(t) >= 21
    assert abs(EE[0] - e_electric_t0(L1, L2, True, excited(L1, L2, s), g)) < 1e-9
    assert abs(EB[0]) < 1e-12
    assert np.abs(EE + EB - (EE[0] + EB[0])).max() < 1e-7, "energy not conserved"
    assert len(r["link_energies"]) == len(t) and len(r["link_energies"][0]) == 3 * N // 2
    th = r["thermal"]
    assert abs(th["electric_energy"] + th["magnetic_energy"] - EE[0]) < 2e-2 * max(1.0, abs(EE[0]))
    for k in ("link_energy", "electric_energy", "magnetic_energy"):
        assert k in r["infinite_temperature"]
    assert isinstance(r["flags"], list)


def test_figures():
    figs = list((RES / "figures").glob("*.png"))
    assert len(figs) >= 6


def test_open_string():
    f = RES / "string_open_3x4_g1.4.json"
    assert f.exists()
    r = json.loads(f.read_text())
    L1, L2, g = 3, 4, 1.4
    exc = [2 * (i * L2 + 2) + k for i in range(L1) for k in (0, 1)]
    t = np.array(r["times"])
    EE, EB = np.array(r["electric_energy"]), np.array(r["magnetic_energy"])
    assert t[-1] >= 4 * g ** 2 - 1e-9
    assert abs(EE[0] - e_electric_t0(L1, L2, False, exc, g)) < 1e-9
    assert np.abs(EE + EB - EE[0]).max() < 1e-7
    assert len(r["z"]) == len(t) and len(r["z"][0]) == 2 * L1 * L2
