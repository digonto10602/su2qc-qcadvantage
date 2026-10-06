"""M4: noisy emulation + mitigation on the N=12 torus -> results/m4/noise_runs.json.

Exact Aer density-matrix runs (shots=None) for every (steps, p2) configuration; the 'repro'
block is a seeded shot-sampling run. Physics: g = 1.25, stripe state, first-order Trotter,
dt = 0.45 g^2 (2 dt / g^2 = 0.9). Usage: python3 scripts/run_m4.py
"""
import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus, initial_state  # noqa: E402
from su2torus.trotter import trotter_circuit  # noqa: E402
from su2torus import hamiltonian as hm, observables as ob  # noqa: E402
from su2torus.noise import (noise_model, run_noisy, predicted_fidelity, echo_fidelity,  # noqa: E402
                            mitigate_rescale, zne, fold_gates, _gate_counts)

G, L1, L2 = 1.25, 3, 2
DT = 0.45 * G ** 2
P1, READOUT = 3e-5, 1e-3
CONFIGS = [(2, 7.9e-4), (4, 1.0e-3), (6, 1.5e-3), (10, 1.0e-3)]
EXTRA = []        # a full (steps x p2) grid costs ~4.5 h of density-matrix time on 4 vCPUs
SEED = 20261006


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def main():
    lat = HoneycombTorus(L1, L2)
    init = initial_state(lat, "stripe")
    dE = hm.electric_diagonal(lat, G)
    O_inf = ob.infinite_temperature(lat, G)["electric_energy"]
    out = {"N": lat.n, "L1": L1, "L2": L2, "g": G, "dt": DT, "order": 1, "state": "stripe",
           "p1": P1, "readout": READOUT,
           "fidelity_formula": "F_pred = prod_gates (1 - eps), eps = p (1 - 1/d^2): "
                               "15 p2/16 per cx/rzz, 3 p1/4 per 1-qubit gate (readout excluded)",
           "mitigation_method": "zne", "configs": [], "extra_configs": []}
    for k, (steps, p2) in enumerate(CONFIGS + EXTRA):
        t0 = time.time()
        evo = trotter_circuit(lat, G, DT, steps, 1, False, None)
        full = trotter_circuit(lat, G, DT, steps, 1, False, init)
        nm = noise_model(p2, P1, READOUT)
        n2, n1 = _gate_counts(evo)
        ideal = float(run_noisy(full, None) @ dE)
        zv = [float(run_noisy(fold_gates(full, f), nm, None, SEED) @ dE) for f in (1, 3, 5)]
        raw = zv[0]
        Fp = predicted_fidelity(evo, p2, P1)
        Fe = echo_fidelity(evo, init, nm, None, SEED)
        c = {"steps": steps, "p2": p2, "shots": None, "seed": SEED, "n_2q": n2, "n_1q": n1,
             "ideal_electric_energy": ideal, "raw_electric_energy": raw,
             "F_predicted": Fp, "F_echo": Fe, "F_rescale": Fe,
             "rescaled_electric_energy": float(mitigate_rescale(raw, Fe, O_inf)),
             "zne_values": zv, "zne_electric_energy": zne(zv, (1, 3, 5))}
        (out["configs"] if k < len(CONFIGS) else out["extra_configs"]).append(c)
        log(f"steps={steps} p2={p2}: ideal={ideal:.4f} raw={raw:.4f} "
            f"resc={c['rescaled_electric_energy']:.4f} zne={c['zne_electric_energy']:.4f} "
            f"Fp={Fp:.4f} Fe={Fe:.4f} ({time.time() - t0:.0f}s)")
    full = trotter_circuit(lat, G, DT, 2, 1, False, init)
    p = run_noisy(full, noise_model(1.0e-3, P1, READOUT), shots=4000, seed=SEED)
    out["repro"] = {"steps": 2, "p2": 1.0e-3, "shots": 4000, "seed": SEED,
                    "raw_electric_energy": float(p @ dE)}
    res = ROOT / "results" / "m4"
    res.mkdir(parents=True, exist_ok=True)
    (res / "noise_runs.json").write_text(json.dumps(out, indent=1))
    figure(out, res)
    log("done")


def figure(out, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = out["configs"] + out["extra_configs"]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    for j, p2 in enumerate((7.9e-4, 1.0e-3, 1.5e-3)):
        rs = sorted([r for r in rows if abs(r["p2"] - p2) < 1e-12], key=lambda r: r["steps"])
        s = [r["steps"] for r in rs]
        for key, lab, st in (("ideal_electric_energy", "ideal (noiseless circuit)", "k-"),
                             ("raw_electric_energy", "raw noisy", "C3o-"),
                             ("rescaled_electric_energy", "echo rescaled", "C1s--"),
                             ("zne_electric_energy", "ZNE (1,3,5)", "C2^-")):
            axs[j].plot(s, [r[key] for r in rs], st, label=lab)
        axs[j].axhline(5.2734375, color="k", ls=":", lw=0.8, label="decohered")
        axs[j].set_title(f"N=12, stripe, g=1.25, p2={p2:g}")
        axs[j].set_xlabel("Trotter steps"); axs[j].set_ylabel("<H_E>")
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    (res / "figures").mkdir(exist_ok=True)
    fig.savefig(res / "figures" / "noise_mitigation.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
