"""Jackknife error bars for the N=24 thermal anchors in results/m3/summary.json.
Re-evaluates typicality at each stored beta with the SAME 4 random vectors (seed 1234), which
must reproduce the stored values, and reports leave-one-vector-out errors
-> results/m3/thermal_errors_n24.json.  Usage: OPENBLAS_NUM_THREADS=1 python3 scripts/thermal_errors_n24.py"""
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus  # noqa: E402
from su2torus import hamiltonian as hm, thermal  # noqa: E402


def main():
    runs = json.loads((ROOT / "results" / "m3" / "summary.json").read_text())["runs"]
    lat = HoneycombTorus(4, 3)
    out = []
    for r in runs:
        if r["N"] != 24:
            continue
        ab = r["model"] == "abelian"
        g, beta = r["g"], r["thermal"]["beta"]
        bounds = json.loads((ROOT / "results" / "m3" / "parts" /
                             f"bounds_n24_{r['model']}.json").read_text())
        diag = hm.electric_diagonal(lat, g)
        t0 = time.time()
        res = thermal.thermal_typicality(lambda v: hm.apply_h(lat, g, v, ab), 2 ** lat.n, beta,
                                         {"electric_energy": lambda v: diag * v}, 4, 1234,
                                         bounds=bounds, recur=hm.cheb_recur(lat, g, ab))
        e = res["errors"]
        rec = {"state": r["state"], "model": r["model"], "beta": beta,
               "electric_energy": res["electric_energy"], "electric_energy_err": e["electric_energy"],
               "energy": res["energy"], "energy_err": e["energy"],
               "stored_electric_energy": r["thermal"]["electric_energy"],
               "reproduces_stored": abs(res["electric_energy"] - r["thermal"]["electric_energy"]) < 1e-8,
               "target_energy": r["electric_energy"][0] + r["magnetic_energy"][0]}
        out.append(rec)
        print(time.strftime("%H:%M:%S"), r["model"], r["state"],
              f"E_E={rec['electric_energy']:.4f}+-{rec['electric_energy_err']:.4f} "
              f"E={rec['energy']:.4f}+-{rec['energy_err']:.4f} target={rec['target_energy']:.4f} "
              f"reproduces={rec['reproduces_stored']} ({time.time() - t0:.0f}s)", flush=True)
    (ROOT / "results" / "m3" / "thermal_errors_n24.json").write_text(json.dumps(
        {"method": "jackknife over the 4 typicality vectors (seed 1234) at the stored beta",
         "runs": out}, indent=1))


if __name__ == "__main__":
    main()
