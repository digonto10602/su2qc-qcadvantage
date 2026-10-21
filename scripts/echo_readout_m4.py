"""M4 follow-up: readout-corrected echo fidelity for the four M4 configurations (exact density
matrix), so the echo is compared with predicted_fidelity on the same footing (readout excluded)
-> results/m4/echo_readout_corrected.json.  Usage: python3 scripts/echo_readout_m4.py"""
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus, initial_state  # noqa: E402
from su2torus.trotter import trotter_circuit  # noqa: E402
from su2torus.noise import noise_model, echo_fidelity_readout_corrected  # noqa: E402


def main():
    d = json.loads((ROOT / "results" / "m4" / "noise_runs.json").read_text())
    lat = HoneycombTorus(d["L1"], d["L2"])
    init = initial_state(lat, d["state"])
    rows = []
    for c in d["configs"]:
        t0 = time.time()
        evo = trotter_circuit(lat, d["g"], d["dt"], c["steps"], 1, False, None)
        fc = echo_fidelity_readout_corrected(evo, init, noise_model(c["p2"], d["p1"], d["readout"]),
                                             None, c["seed"])
        rows.append({"steps": c["steps"], "p2": c["p2"], "F_predicted": c["F_predicted"],
                     "F_echo": c["F_echo"], "F_echo_readout_corrected": fc,
                     "rel_dev_raw": c["F_echo"] / c["F_predicted"] - 1,
                     "rel_dev_corrected": fc / c["F_predicted"] - 1})
        print(time.strftime("%H:%M:%S"), c["steps"], c["p2"], f"Fp={c['F_predicted']:.4f} "
              f"Fe={c['F_echo']:.4f} Fe_corr={fc:.4f} ({time.time() - t0:.0f}s)", flush=True)
    (ROOT / "results" / "m4" / "echo_readout_corrected.json").write_text(json.dumps(
        {"note": "echo return probability divided by (1-r)^n, r = readout flip rate, n = 12",
         "configs": rows}, indent=1))


if __name__ == "__main__":
    main()
