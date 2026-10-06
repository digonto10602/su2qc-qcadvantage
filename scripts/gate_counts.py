"""M2: true two-qubit counts per plaquette per Trotter step -> results/m2/gate_counts.json."""
import json
import pathlib
from su2torus.lattice import HoneycombTorus
from su2torus.trotter import trotter_step, trotter_circuit, two_qubit_counts

ROOT = pathlib.Path(__file__).resolve().parents[1]
rows = []
for (L1, L2, per) in [(3, 2, True), (4, 3, True), (4, 6, True), (3, 4, False)]:
    lat = HoneycombTorus(L1, L2, per)
    for model in ("su2", "abelian"):
        for order in (1, 2):
            c1 = two_qubit_counts(trotter_step(lat, 1.25, 0.3, order, model == "abelian"))
            row = {"L1": L1, "L2": L2, "periodic": per, "n": lat.n, "model": model,
                   "order": order, "per_step": c1,
                   "per_plaquette_per_step": {k: v / lat.n for k, v in c1.items()}}
            if order == 2 and model == "su2":
                # ...U_A U_E | U_E U_A...: only the two electric half-layers merge
                fused = {"native_zz": c1["native_zz"] - len(lat.bonds),
                         "cz": c1["cz"] - 2 * len(lat.bonds)}
                row["fused_per_plaquette_per_step"] = {k: v / lat.n for k, v in fused.items()}
            rows.append(row)
out = ROOT / "results" / "m2" / "gate_counts.json"
out.write_text(json.dumps(rows, indent=1))
for r in rows:
    print(r["L1"], r["L2"], r["periodic"], r["model"], r["order"], r["per_plaquette_per_step"],
          r.get("fused_per_plaquette_per_step", ""))
