import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]


def _run(args):
    return subprocess.run([sys.executable, "scripts/run_classical.py", *args], cwd=ROOT,
                          capture_output=True, text=True, timeout=600)


def test_resume_keeps_finished_chi(tmp_path):
    cfg = {"name": "resume", "L1": 2, "L2": 2, "g": 1.25, "dt": 0.3, "steps": 1, "order": 1,
           "abelian": False, "state": "stripe", "chi": [2, 4], "cutoff": 1e-12,
           "backend": "numpy", "output": "unused.json"}
    c = tmp_path / "c.json"
    c.write_text(json.dumps(cfg))
    out = tmp_path / "o.json"
    assert _run(["--config", str(c), "--out", str(out)]).returncode == 0
    full = json.loads(out.read_text())
    partial = dict(full, runs=[r for r in full["runs"] if r["chi"] == 4])   # "chi=2 lost"
    partial["runs"][0]["wall_time_s"] = -1.0                                   # marker
    out.write_text(json.dumps(partial))
    p = _run(["--config", str(c), "--out", str(out)])
    assert p.returncode == 0 and "resuming" in p.stdout
    d = json.loads(out.read_text())
    assert [r["chi"] for r in d["runs"]] == [2, 4]
    assert d["runs"][1]["wall_time_s"] == -1.0          # chi=4 kept, not recomputed
    cfg2 = dict(cfg, steps=2)
    c.write_text(json.dumps(cfg2))
    assert _run(["--config", str(c), "--out", str(out)]).returncode != 0   # config mismatch
    assert _run(["--config", str(c), "--out", str(out), "--restart"]).returncode == 0
