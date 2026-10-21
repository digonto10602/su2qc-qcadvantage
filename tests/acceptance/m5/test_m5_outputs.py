"""M5 (a)-(c): committed outputs, resource table, config entry point and Slurm templates.

Required files and schemas (in addition to su2torus/mps.py, see test_m5_mps.py)
--------------------------------------------------------------------------------
scripts/run_m5_check.py  -> results/m5/mps_check.json   (run once with nohup; the N = 24, chi = 256
    MPS took > 30 min on the loaded 4-vCPU cloud machine in the reviewer's prototype, too slow for CI)
    {"chi": 256, "steps": 4, "order": 1, "model": "su2", "g": 1.25, "dt": 0.2, "state": "stripe",
     "cutoff": float,
     "runs": [{"N": 12, "L1": 3, "L2": 2, ...}, {"N": 24, "L1": 4, "L2": 3, ...}]}
    each run: "fidelity"        : |<exact|mps>|^2 / (<exact|exact><mps|mps>) of the final states, the
                                  exact one being the statevector of the same Trotter circuit
              "error_estimate"  : float (simulate_mps), "max_bond": int, "wall_time_s": float,
              "times"           : steps+1 floats,
              "electric_energy" : steps+1 floats (MPS),  "z": steps+1 lists of N floats (MPS,
                                  plaquette order)

scripts/resource_table.py [--out PATH]   (default PATH = results/m5/resources.json; deterministic)
    Compiles first-order su2 Trotter circuits (su2torus.trotter) for the tori below and writes
    {"order": 1, "model": "su2",
     "p2_values": [7.9e-4, 1e-3, 1.5e-3],
     "steps": [4, 5, 6, 7, 8, 9, 10],
     "tori": [{"N": 48, "L1": 4, "L2": 6, "per_step": {"native_zz": int, "cz": int}}, ... N = 56 (4x7),
              64 (4x8), 96 (6x8)],
     "rows": [{"N": int, "steps": int, "gate_set": "native_zz" | "cz", "n_2q": int, "p2": float,
               "raw_fidelity": exp(-p2 * n_2q)}, ...]}      one row per (torus, steps, gate_set, p2)
    Gate-count convention: METHODS section 3 (native ZZ: cx = rzz = 1; CZ: cx = 1, rzz = 2).

scripts/run_classical.py --config PATH [--out PATH] [--dry-run]
    Config JSON (configs/*.json; at least one committed):
    {"name": str, "L1": int, "L2": int, "g": float, "dt": float, "steps": int, "order": 1 | 2,
     "abelian": bool, "state": "vacuum" | "neel_half" | "stripe", "chi": [int, ...],
     "cutoff": float, "backend": "numpy" | "cupy" | "torch",
     "output": "results/.../name.json"}            (relative paths are relative to the repo root)
    Runs su2torus.mps.simulate_mps once per chi and writes to --out (else config "output"):
    {"config": {...the config...},
     "runs": [{"chi": int, "times": [...], "electric_energy": [...], "z": [[...]],
               "error_estimate": float, "error_estimate_per_step": [...], "max_bond": int,
               "wall_time_s": float}, ...]}            (runs in the order of config "chi")
    --dry-run: validate the config, print the planned runs, write nothing, exit 0.
    A config missing a required key -> non-zero exit status.

scripts/slurm/*.sbatch  (at least one)
    Site-specific values only as placeholders <PARTITION>, <ACCOUNT>, <ENV> (other <UPPER_CASE>
    placeholders allowed); calls scripts/run_classical.py --config; no real cluster names.
"""
import json
import math
import pathlib
import re
import subprocess
import sys

import numpy as np
import pytest

from tests.acceptance._ref import ref_bonds, ref_neighbors
from tests.acceptance.m5._exact_m5 import exact_trotter, stripe, z_and_he

ROOT = pathlib.Path(__file__).resolve().parents[3]
M5 = ROOT / "results" / "m5"
G, DT, STEPS, CHI = 1.25, 0.2, 4, 256
TORI = {48: (4, 6), 56: (4, 7), 64: (4, 8), 96: (6, 8)}
P2 = [7.9e-4, 1e-3, 1.5e-3]
STEP_RANGE = list(range(4, 11))


def _load(name):
    path = M5 / name
    assert path.exists(), f"missing {path}"
    return json.loads(path.read_text())


def _run(args, timeout=600):
    return subprocess.run([sys.executable, *args], cwd=ROOT, capture_output=True, text=True,
                          timeout=timeout)


# ---------------------------------------------------------------- (a) MPS vs exact, N = 12 and 24
@pytest.fixture(scope="module")
def check():
    d = _load("mps_check.json")
    assert d["chi"] == CHI and d["steps"] == STEPS and d["order"] == 1
    assert d["model"] == "su2" and d["state"] == "stripe"
    assert abs(d["g"] - G) < 1e-12 and abs(d["dt"] - DT) < 1e-12
    return {r["N"]: r for r in d["runs"]}


def test_mps_check_fidelity(check):
    """PLAN M5: MPS = exact statevector at N = 12 and 24, 4 steps, chi = 256, F > 1 - 1e-6."""
    assert (ROOT / "scripts" / "run_m5_check.py").exists()
    assert set(check) >= {12, 24}
    for N, (a, b) in ((12, (3, 2)), (24, (4, 3))):
        r = check[N]
        assert (r["L1"], r["L2"]) == (a, b)
        assert r["fidelity"] > 1 - 1e-6, (N, r["fidelity"])
        assert 0.0 <= r["error_estimate"] <= 1.0 and r["max_bond"] <= CHI
        assert len(r["times"]) == STEPS + 1 and len(r["electric_energy"]) == STEPS + 1
        assert len(r["z"]) == STEPS + 1 and all(len(z) == N for z in r["z"])


@pytest.mark.parametrize("N", [12, 24])
def test_mps_check_observables_independent(check, N):
    """Independent cross-check of the committed MPS observables against an exact statevector
    built here from METHODS (no su2torus code). If F > 1 - eps then for any observable O,
    |<O>_mps - <O>_exact| <= 2 ||O - c|| sqrt(eps) (trace-distance bound), eps = 1e-6:
      Z_p : ||Z|| = 1                     -> 2e-3
      H_E : centred norm (3g^2/16) n_bonds -> 2e-3 (3g^2/16) n_bonds."""
    a, b = (3, 2) if N == 12 else (4, 3)
    r = check[N]
    ref = exact_trotter(a, b, G, DT, STEPS, stripe(a, b), observe=lambda v: z_and_he(v, a, b, G))
    nbonds = len(ref_bonds(ref_neighbors(a, b, True)))
    tol_z, tol_he = 2e-3, 2e-3 * (3 * G ** 2 / 16) * nbonds
    np.testing.assert_allclose(r["times"], [k * DT for k in range(STEPS + 1)], atol=1e-12)
    for k, (z_ref, he_ref) in enumerate(ref):
        np.testing.assert_allclose(r["z"][k], z_ref, atol=tol_z, err_msg=f"N={N} step {k}")
        assert abs(r["electric_energy"][k] - he_ref) <= tol_he, (N, k)


# ---------------------------------------------------------------- (b) resource table
@pytest.fixture(scope="module")
def resources():
    return _load("resources.json")


def _m2_per_plaquette():
    """First-order su2 torus counts per plaquette per step from results/m2/gate_counts.json."""
    m2 = json.loads((ROOT / "results" / "m2" / "gate_counts.json").read_text())
    vals = {(e["per_plaquette_per_step"]["native_zz"], e["per_plaquette_per_step"]["cz"])
            for e in m2 if e["model"] == "su2" and e["order"] == 1 and e["periodic"]}
    assert len(vals) == 1, f"M2 per-plaquette counts differ between tori: {vals}"
    return vals.pop()


def test_resource_table_per_step(resources):
    """PLAN M5: the table matches the M2 per-plaquette counts x N (x steps, below)."""
    from su2torus.lattice import HoneycombTorus
    from su2torus.trotter import trotter_step, two_qubit_counts
    zz_pp, cz_pp = _m2_per_plaquette()
    assert resources["order"] == 1 and resources["model"] == "su2"
    np.testing.assert_allclose(resources["p2_values"], P2, rtol=0, atol=1e-15)
    assert list(resources["steps"]) == STEP_RANGE
    tori = {t["N"]: t for t in resources["tori"]}
    assert set(tori) == set(TORI)
    for N, (a, b) in TORI.items():
        t = tori[N]
        assert (t["L1"], t["L2"]) == (a, b)
        assert t["per_step"]["native_zz"] == zz_pp * N
        assert t["per_step"]["cz"] == cz_pp * N
        # independent count (METHODS section 3): 2**k cx per flip (k = 3 neighbours on the
        # torus) and one rzz per bond
        nbonds = len(ref_bonds(ref_neighbors(a, b, True)))
        assert t["per_step"]["native_zz"] == 8 * N + nbonds
        assert t["per_step"]["cz"] == 8 * N + 2 * nbonds
        # the table comes from the actual compiled circuit of this torus
        assert two_qubit_counts(trotter_step(HoneycombTorus(a, b), G, 0.1, order=1)) == t["per_step"]


def test_resource_table_rows(resources):
    """Totals = per-step count x steps; raw fidelity = exp(-p2 N_2q), for every combination."""
    per = {t["N"]: t["per_step"] for t in resources["tori"]}
    rows = {(r["N"], r["steps"], r["gate_set"], round(r["p2"], 10)): r for r in resources["rows"]}
    want = {(N, s, gs, round(p, 10)) for N in TORI for s in STEP_RANGE
            for gs in ("native_zz", "cz") for p in P2}
    assert set(rows) == want and len(resources["rows"]) == len(want)
    for (N, s, gs, p), r in rows.items():
        assert r["n_2q"] == per[N][gs] * s
        # 1e-9 relative: JSON float round-trip only
        assert math.isclose(r["raw_fidelity"], math.exp(-p * r["n_2q"]), rel_tol=1e-9, abs_tol=0)


def test_resource_table_script_reproduces(tmp_path, resources):
    out = tmp_path / "resources.json"
    proc = _run(["scripts/resource_table.py", "--out", str(out)], timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    new = json.loads(out.read_text())
    for key in ("order", "model", "p2_values", "steps", "tori", "rows"):
        assert new[key] == resources[key], key


# ---------------------------------------------------------------- (c) config-driven entry point
def _configs():
    return sorted((ROOT / "configs").glob("*.json"))


REQUIRED = ("name", "L1", "L2", "g", "dt", "steps", "order", "abelian", "state", "chi", "cutoff",
            "backend", "output")


def test_committed_configs_dry_run(tmp_path):
    cfgs = _configs()
    assert cfgs, "no configs/*.json committed"
    for cfg in cfgs:
        c = json.loads(cfg.read_text())
        for key in REQUIRED:
            assert key in c, (cfg.name, key)
        assert c["state"] in ("vacuum", "neel_half", "stripe") and c["order"] in (1, 2)
        assert isinstance(c["chi"], list) and all(int(x) > 0 for x in c["chi"])
        out = tmp_path / (cfg.stem + ".json")
        proc = _run(["scripts/run_classical.py", "--config", str(cfg), "--out", str(out),
                     "--dry-run"], timeout=120)
        assert proc.returncode == 0, (cfg.name, proc.stderr[-2000:])
        assert not out.exists(), "--dry-run must not write output"


def test_run_classical_small_config(tmp_path):
    a, b = 3, 2
    cfg = {"name": "acceptance_n12", "L1": a, "L2": b, "g": G, "dt": DT, "steps": 2, "order": 1,
           "abelian": False, "state": "stripe", "chi": [4, 64], "cutoff": 1e-12,
           "backend": "numpy", "output": "results/m5/should_not_be_written.json"}
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(cfg))
    out = tmp_path / "out.json"
    proc = _run(["scripts/run_classical.py", "--config", str(path), "--out", str(out)])
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert not (ROOT / cfg["output"]).exists(), "--out must override the config output path"
    d = json.loads(out.read_text())
    assert d["config"]["name"] == cfg["name"]
    assert [r["chi"] for r in d["runs"]] == [4, 64]
    ref = exact_trotter(a, b, G, DT, 2, stripe(a, b), observe=lambda v: z_and_he(v, a, b, G))
    exact_run = d["runs"][1]          # chi = 64 >= 2**6: exact at n = 12 (1e-6 margin as in test_m5_mps)
    for k, (z_ref, he_ref) in enumerate(ref):
        np.testing.assert_allclose(exact_run["z"][k], z_ref, atol=1e-6)
        assert abs(exact_run["electric_energy"][k] - he_ref) < 1e-6
    assert d["runs"][0]["error_estimate"] >= d["runs"][1]["error_estimate"] - 1e-9


def test_run_classical_rejects_bad_config(tmp_path):
    assert (ROOT / "scripts" / "run_classical.py").exists()
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"name": "bad", "L1": 3, "L2": 2, "steps": 1}))
    proc = _run(["scripts/run_classical.py", "--config", str(path), "--out",
                 str(tmp_path / "o.json"), "--dry-run"], timeout=120)
    assert proc.returncode != 0


# ---------------------------------------------------------------- Slurm templates
FACILITIES = r"\b(nersc|olcf|alcf|tacc|ncsa|perlmutter|cori|summit|polaris|aurora|bridges2?" \
             r"|expanse|stampede\d?|frontera|lumi|juwels|jureca|leonardo|archer2?|marenostrum" \
             r"|isambard|crusher|sunspot)\b"
SITE_OPTS = re.compile(r"^#SBATCH\s+(--partition|-p|--account|-A|--qos|-q|--constraint|-C"
                       r"|--reservation|--clusters|-M|--mail-user)(=|\s+)(\S+)")


def test_slurm_templates_placeholders_only():
    """PLAN M5: Slurm templates contain only placeholders (<PARTITION>, <ACCOUNT>, <ENV>)."""
    files = sorted((ROOT / "scripts" / "slurm").glob("*.sbatch"))
    assert files, "no scripts/slurm/*.sbatch"
    for f in files:
        text = f.read_text()
        for ph in ("<PARTITION>", "<ACCOUNT>", "<ENV>"):
            assert ph in text, (f.name, ph)
        assert re.search(r"^#SBATCH\s+(--partition|-p)(=|\s+)<PARTITION>\s*$", text, re.M), f.name
        assert re.search(r"^#SBATCH\s+(--account|-A)(=|\s+)<ACCOUNT>\s*$", text, re.M), f.name
        assert any("<ENV>" in ln for ln in text.splitlines() if not ln.startswith("#")), f.name
        assert "run_classical.py" in text and "--config" in text, f.name
        for ln in text.splitlines():
            m = SITE_OPTS.match(ln.strip())
            if m:
                assert re.fullmatch(r"<[A-Z_]+>", m.group(3)), (f.name, ln)
        assert not re.search(FACILITIES, text, re.I), (f.name, re.search(FACILITIES, text, re.I))
        assert not re.search(r"[\w.+-]+@[\w-]+\.[\w.]+", text), (f.name, "e-mail address")
