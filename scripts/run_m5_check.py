"""M5 check: MPS (chi=256) vs exact statevector of the same Trotter circuit, N=12 and N=24,
4 first-order steps, dt=0.2, g=1.25, stripe -> results/m5/mps_check.json. Run with nohup."""
import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus, initial_state  # noqa: E402
from su2torus.trotter import trotter_circuit  # noqa: E402
from su2torus.mps import simulate_mps, mps_statevector  # noqa: E402

G, DT, STEPS, CHI, CUTOFF = 1.25, 0.2, 4, 256, 1e-12


def exact_state(lat, init):
    from qiskit_aer import AerSimulator
    qc = trotter_circuit(lat, G, DT, STEPS, 1, False, init)
    qc.save_statevector()
    res = AerSimulator(method="statevector").run(qc).result()
    return np.asarray(res.get_statevector(), dtype=np.complex128)


def main():
    runs = []
    for L1, L2 in ((3, 2), (4, 3)):
        lat = HoneycombTorus(L1, L2)
        init = initial_state(lat, "stripe")
        res = simulate_mps(lat, G, DT, STEPS, initial=init, chi=CHI, cutoff=CUTOFF)
        psi = mps_statevector(res)
        res.pop("mps")
        ref = exact_state(lat, init)
        fid = abs(np.vdot(ref, psi)) ** 2 / (np.vdot(ref, ref).real * np.vdot(psi, psi).real)
        del psi, ref
        runs.append({"N": lat.n, "L1": L1, "L2": L2, "fidelity": float(fid),
                     "error_estimate": res["error_estimate"], "max_bond": res["max_bond"],
                     "times": res["times"], "electric_energy": res["electric_energy"],
                     "z": res["z"], "wall_time_s": res["wall_time_s"]})
        print(time.strftime("%H:%M:%S"), f"N={lat.n}: 1-F={1 - fid:.3e} "
              f"err_est={res['error_estimate']:.3e} max_bond={res['max_bond']} "
              f"({res['wall_time_s']:.0f}s)", flush=True)
    out = {"chi": CHI, "steps": STEPS, "order": 1, "model": "su2", "g": G, "dt": DT,
           "state": "stripe", "cutoff": CUTOFF, "runs": runs}
    (ROOT / "results" / "m5").mkdir(parents=True, exist_ok=True)
    (ROOT / "results" / "m5" / "mps_check.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
