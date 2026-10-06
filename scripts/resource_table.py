"""M5: two-qubit gate resources of first-order SU(2) Trotter circuits on the frontier tori
-> results/m5/resources.json (deterministic). Usage: python3 scripts/resource_table.py [--out P]"""
import argparse
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus  # noqa: E402
from su2torus.trotter import trotter_step, two_qubit_counts  # noqa: E402

TORI = [(48, 4, 6), (56, 4, 7), (64, 4, 8), (96, 6, 8)]
P2 = [7.9e-4, 1e-3, 1.5e-3]          # Helios, H2, Heron (PLAN M4)
STEPS = list(range(4, 11))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results" / "m5" / "resources.json"))
    args = ap.parse_args()
    tori, rows = [], []
    for N, L1, L2 in TORI:
        lat = HoneycombTorus(L1, L2)
        assert lat.n == N
        per = two_qubit_counts(trotter_step(lat, 1.25, 0.1, order=1))   # counts independent of g, dt
        tori.append({"N": N, "L1": L1, "L2": L2, "per_step": per})
        for s in STEPS:
            for gs in ("native_zz", "cz"):
                for p in P2:
                    n2 = per[gs] * s
                    rows.append({"N": N, "steps": s, "gate_set": gs, "n_2q": n2, "p2": p,
                                 "raw_fidelity": math.exp(-p * n2)})
    out = {"order": 1, "model": "su2", "p2_values": P2, "steps": STEPS, "tori": tori,
           "rows": rows}
    path = pathlib.Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    for t in tori:
        print(t["N"], t["per_step"])


if __name__ == "__main__":
    main()
