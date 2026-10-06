"""M6: per-step chi-convergence record (`results/m6/convergence.json`).

PLAN M6: "Record per step where observables stop converging in chi."

New script: scripts/analyze_convergence.py [--in PATH] [--out PATH]
    defaults: --in results/m6/n48_stripe_g1.25.json, --out results/m6/convergence.json.
    Pure post-processing of the MPS run file (no MPS simulation), deterministic.

Convergence criterion (fixed here; the JSON must state it in words in "criterion"):
    chi list c_0 < c_1 < c_2 < c_3 = [32, 64, 128, 256]; successive pairs i = 0, 1, 2 are
    (c_i, c_{i+1}). At step k (k = 0..8):
      d_electric_energy[i] = |E_{c_{i+1}}(k) - E_{c_i}(k)|                  (E = <H_E>, full)
      d_z_max[i]           = max_p |<Z_p>_{c_{i+1}}(k) - <Z_p>_{c_i}(k)|     (all 48 plaquettes)
      converged_electric_energy[i] = d_electric_energy[i] <= TOL_E
      converged_z[i]               = d_z_max[i] <= TOL_Z
      converged[i]                 = both
    TOL_Z = 0.01: about the statistical error 1/sqrt(10^4) of a hardware <Z_p> estimate with
            10^4 shots, i.e. the precision at which the classical value is needed.
    TOL_E = 0.01 * E_inf, E_inf = (g^2/2)(3/8) n_bonds = 21.09375 (decohered value, METHODS
            sec. 4): 1% of the decohered value (M3 flags use 10% of the same scale).
    last_converged_step*[i] = largest K in 0..8 with converged*[i] true at every step 0..K
            (-1 if step 0 fails); first_unconverged_step*[i] = that K + 1, or null if K = 8.
    Headline: last_converged_step[2] = last step at which chi = 256 agrees with chi = 128.

Schema:
{
  "source": "results/m6/n48_stripe_g1.25.json",
  "N": 48, "L1": 4, "L2": 6, "g": 1.25, "dt": float, "order": 1, "state": "stripe",
  "chi": [32, 64, 128, 256],
  "pairs": [[32, 64], [64, 128], [128, 256]],
  "tol_electric_energy": 0.2109375, "tol_z": 0.01,
  "criterion": "<plain-English statement of the rule above>",
  "steps": [{"step": k, "time": k * dt,
             "electric_energy": [E per chi], "z_mean": [mean_p <Z_p> per chi],
             "error_estimate": [error_estimate_per_step[k] per chi],
             "d_electric_energy": [3], "d_z_max": [3],
             "converged_electric_energy": [3 bools], "converged_z": [3 bools],
             "converged": [3 bools]} for k = 0..8],
  "last_converged_step_electric_energy": [3 ints], "last_converged_step_z": [3 ints],
  "last_converged_step": [3 ints],
  "first_unconverged_step": [3 ints or null]
}
All quantities are recomputed below directly from the MPS run file (independent of the script).
"""
import json
import math
import pathlib
import subprocess
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[3]
SRC = "results/m6/n48_stripe_g1.25.json"
CONV = ROOT / "results" / "m6" / "convergence.json"
CHIS = [32, 64, 128, 256]
G, N_BONDS, STEPS = 1.25, 72, 8
E_INF = G ** 2 / 2 * 3 / 8 * N_BONDS
TOL_E, TOL_Z = 0.01 * E_INF, 0.01
FTOL = 1e-12          # recomputed floats: same arithmetic on the same JSON numbers, round-off only


def _load(path):
    assert path.exists(), f"missing {path}"
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def raw():
    return _load(ROOT / SRC)


@pytest.fixture(scope="module")
def conv():
    return _load(CONV)


def _last(flags):
    k = -1
    for f in flags:
        if not f:
            break
        k += 1
    return k


def reference(raw):
    runs = raw["runs"]
    assert [r["chi"] for r in runs] == CHIS
    steps = []
    for k in range(STEPS + 1):
        E = [r["electric_energy"][k] for r in runs]
        Z = [np.asarray(r["z"][k], dtype=float) for r in runs]
        dE = [float(abs(E[i + 1] - E[i])) for i in range(3)]
        dZ = [float(np.max(np.abs(Z[i + 1] - Z[i]))) for i in range(3)]
        cE = [bool(d <= TOL_E) for d in dE]
        cZ = [bool(d <= TOL_Z) for d in dZ]
        steps.append({"electric_energy": E, "z_mean": [float(z.mean()) for z in Z],
                      "error_estimate": [r["error_estimate_per_step"][k] for r in runs],
                      "d_electric_energy": dE, "d_z_max": dZ,
                      "converged_electric_energy": cE, "converged_z": cZ,
                      "converged": [a and b for a, b in zip(cE, cZ)]})
    out = {}
    for key, flag in (("last_converged_step_electric_energy", "converged_electric_energy"),
                      ("last_converged_step_z", "converged_z"), ("last_converged_step", "converged")):
        out[key] = [_last([s[flag][i] for s in steps]) for i in range(3)]
    out["first_unconverged_step"] = [None if K == STEPS else K + 1 for K in out["last_converged_step"]]
    return steps, out


def test_header(conv, raw):
    assert conv["source"] == SRC
    assert (conv["N"], conv["L1"], conv["L2"]) == (48, 4, 6)
    assert conv["state"] == "stripe" and conv["order"] == 1
    assert abs(conv["g"] - G) < 1e-12 and abs(conv["dt"] - raw["config"]["dt"]) < 1e-12
    assert conv["chi"] == CHIS and conv["pairs"] == [[32, 64], [64, 128], [128, 256]]
    assert abs(conv["tol_electric_energy"] - TOL_E) < 1e-12 and abs(conv["tol_z"] - TOL_Z) < 1e-12
    assert isinstance(conv["criterion"], str) and len(conv["criterion"]) >= 40


def _near_threshold(d, tol):
    return abs(d - tol) <= 1e-9          # a decision this close to the threshold is round-off


def test_per_step_record(conv, raw):
    """Every step 0..8 recorded; differences and flags equal the independent recomputation."""
    ref_steps, _ = reference(raw)
    dt = raw["config"]["dt"]
    assert [s["step"] for s in conv["steps"]] == list(range(STEPS + 1))
    for s, r in zip(conv["steps"], ref_steps):
        k = s["step"]
        assert abs(s["time"] - k * dt) < 1e-12
        for key in ("electric_energy", "z_mean", "error_estimate", "d_electric_energy", "d_z_max"):
            np.testing.assert_allclose(s[key], r[key], rtol=0, atol=FTOL, err_msg=f"step {k} {key}")
        for key, dkey, tol in (("converged_electric_energy", "d_electric_energy", TOL_E),
                               ("converged_z", "d_z_max", TOL_Z)):
            for i in range(3):
                if not _near_threshold(r[dkey][i], tol):
                    assert s[key][i] is r[key][i], (k, key, i)
        assert s["converged"] == [a and b for a, b in
                                  zip(s["converged_electric_energy"], s["converged_z"])], k


def test_last_converged_steps(conv, raw):
    """Where the observables stop converging in chi (per pair; headline = pair (128, 256))."""
    steps = conv["steps"]
    for key, flag in (("last_converged_step_electric_energy", "converged_electric_energy"),
                      ("last_converged_step_z", "converged_z"), ("last_converged_step", "converged")):
        assert conv[key] == [_last([s[flag][i] for s in steps]) for i in range(3)], key
    assert conv["first_unconverged_step"] == [None if K == STEPS else K + 1
                                              for K in conv["last_converged_step"]]
    ref_steps, ref = reference(raw)
    flags_safe = all(not _near_threshold(r[d][i], t) for r in ref_steps
                     for d, t in (("d_electric_energy", TOL_E), ("d_z_max", TOL_Z)) for i in range(3))
    if flags_safe:
        for key in ref:
            assert conv[key] == ref[key], key


def test_analyze_script_reproduces(tmp_path):
    """scripts/analyze_convergence.py regenerates the committed file (deterministic)."""
    out = tmp_path / "conv.json"
    p = subprocess.run([sys.executable, "scripts/analyze_convergence.py", "--in", SRC, "--out", str(out)],
                       cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr[-2000:]
    new, old = json.loads(out.read_text()), _load(CONV)

    def same(a, b):
        if isinstance(a, float) or isinstance(b, float):
            return a is not None and b is not None and math.isclose(a, b, rel_tol=0, abs_tol=FTOL)
        if isinstance(a, list):
            return isinstance(b, list) and len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
        if isinstance(a, dict):
            return isinstance(b, dict) and a.keys() == b.keys() and all(same(a[x], b[x]) for x in a)
        return a == b

    assert same(new, old)
