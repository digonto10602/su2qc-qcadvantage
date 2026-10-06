"""M6: per-step chi convergence of the N=48 MPS run -> results/m6/convergence.json.
Usage: python3 scripts/analyze_convergence.py [--in PATH] [--out PATH]   (deterministic)"""
import argparse
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
TOL_Z = 0.01


def last_true_prefix(flags):
    K = -1
    for f in flags:
        if not f:
            break
        K += 1
    return K


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="results/m6/n48_stripe_g1.25.json")
    ap.add_argument("--out", default=str(ROOT / "results" / "m6" / "convergence.json"))
    a = ap.parse_args()
    src = pathlib.Path(a.inp)
    d = json.loads((src if src.is_absolute() else ROOT / src).read_text())
    cfg = d["config"]
    runs = sorted(d["runs"], key=lambda r: r["chi"])
    chis = [r["chi"] for r in runs]
    n = 2 * cfg["L1"] * cfg["L2"]
    e_inf = cfg["g"] ** 2 / 2 * 3 / 8 * (3 * n // 2)
    tol_e = 0.01 * e_inf
    pairs = [[chis[i], chis[i + 1]] for i in range(len(chis) - 1)]
    steps = []
    for k in range(cfg["steps"] + 1):
        E = [r["electric_energy"][k] for r in runs]
        Z = [np.array(r["z"][k]) for r in runs]
        dE = [abs(E[i + 1] - E[i]) for i in range(len(pairs))]
        dZ = [float(np.abs(Z[i + 1] - Z[i]).max()) for i in range(len(pairs))]
        ce = [x <= tol_e for x in dE]
        cz = [x <= TOL_Z for x in dZ]
        steps.append({"step": k, "time": k * cfg["dt"], "electric_energy": E,
                      "z_mean": [float(z.mean()) for z in Z],
                      "error_estimate": [r["error_estimate_per_step"][k] for r in runs],
                      "d_electric_energy": dE, "d_z_max": dZ,
                      "converged_electric_energy": ce, "converged_z": cz,
                      "converged": [x and y for x, y in zip(ce, cz)]})
    lc = {key: [last_true_prefix([s[key][i] for s in steps]) for i in range(len(pairs))]
          for key in ("converged_electric_energy", "converged_z", "converged")}
    out = {"source": str(a.inp), "N": n, "L1": cfg["L1"], "L2": cfg["L2"], "g": cfg["g"],
           "dt": cfg["dt"], "order": cfg["order"], "state": cfg["state"], "chi": chis,
           "pairs": pairs, "tol_electric_energy": tol_e, "tol_z": TOL_Z,
           "criterion": ("At each Trotter step k, successive bond dimensions (chi_i, chi_i+1) "
                         "agree if |<H_E>| differs by at most 1% of the decohered value "
                         f"({tol_e:g}) and every <Z_p> differs by at most {TOL_Z}. "
                         "last_converged_step[i] is the last step K such that the pair agrees "
                         "at all steps 0..K (-1 if step 0 fails)."),
           "steps": steps,
           "last_converged_step_electric_energy": lc["converged_electric_energy"],
           "last_converged_step_z": lc["converged_z"],
           "last_converged_step": lc["converged"],
           "first_unconverged_step": [K + 1 if K < cfg["steps"] else None
                                      for K in lc["converged"]]}
    p = pathlib.Path(a.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1))
    print("last converged step per pair", pairs, out["last_converged_step"])


if __name__ == "__main__":
    main()
