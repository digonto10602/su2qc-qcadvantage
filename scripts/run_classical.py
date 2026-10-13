"""Config-driven MPS runs (cluster entry point).
Usage: python3 scripts/run_classical.py --config configs/X.json [--out PATH] [--dry-run]
Config keys: name, L1, L2, g, dt, steps, order (1|2), abelian, state, chi (list), cutoff,
backend ('numpy'|'cupy'|'torch'), output (path, relative to the repo root)."""
import argparse
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

REQUIRED = {"name": str, "L1": int, "L2": int, "g": (int, float), "dt": (int, float),
            "steps": int, "order": int, "abelian": bool, "state": str, "chi": list,
            "cutoff": (int, float), "backend": str, "output": str}


def validate(cfg):
    errs = [f"missing key {k!r}" for k in REQUIRED if k not in cfg]
    errs += [f"key {k!r} has wrong type" for k, t in REQUIRED.items()
             if k in cfg and not isinstance(cfg[k], t)]
    if not errs:
        if cfg["state"] not in ("vacuum", "neel_half", "stripe"):
            errs.append("state must be vacuum, neel_half or stripe")
        if cfg["order"] not in (1, 2):
            errs.append("order must be 1 or 2")
        if cfg["backend"] not in ("numpy", "cupy", "torch"):
            errs.append("backend must be numpy, cupy or torch")
        if not cfg["chi"] or not all(isinstance(c, int) and c > 0 for c in cfg["chi"]):
            errs.append("chi must be a non-empty list of positive ints")
        if min(cfg["L1"], cfg["L2"]) < 2 or cfg["steps"] < 0:
            errs.append("need L1, L2 >= 2 and steps >= 0")
    return errs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--restart", action="store_true",
                    help="ignore an existing output file instead of resuming from it")
    args = ap.parse_args()
    cfg = json.loads(pathlib.Path(args.config).read_text())
    errs = validate(cfg)
    if errs:
        print("invalid config:", "; ".join(errs), file=sys.stderr)
        sys.exit(2)
    out = pathlib.Path(args.out) if args.out else ROOT / cfg["output"]
    n = 2 * cfg["L1"] * cfg["L2"]
    print(f"{cfg['name']}: N={n} g={cfg['g']} dt={cfg['dt']} steps={cfg['steps']} "
          f"order={cfg['order']} state={cfg['state']} backend={cfg['backend']} chi={cfg['chi']} "
          f"-> {out}", flush=True)
    if args.dry_run:
        return
    from su2torus.lattice import HoneycombTorus, initial_state
    from su2torus.mps import simulate_mps
    lat = HoneycombTorus(cfg["L1"], cfg["L2"])
    init = initial_state(lat, cfg["state"])
    result = {"config": cfg, "runs": []}
    if out.exists() and not args.restart:
        prev = json.loads(out.read_text())
        if prev.get("config") == cfg:                 # resume: keep finished chi values
            result["runs"] = [r for r in prev.get("runs", []) if r.get("chi") in cfg["chi"]]
            print(f"resuming: chi already done {[r['chi'] for r in result['runs']]}", flush=True)
        else:
            sys.exit(f"{out} exists with a different config; use --restart or another --out")
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {r["chi"] for r in result["runs"]}
    for chi in cfg["chi"]:
        if chi in done:
            continue
        r = simulate_mps(lat, cfg["g"], cfg["dt"], cfg["steps"], order=cfg["order"],
                         abelian=cfg["abelian"], initial=init, chi=chi, cutoff=cfg["cutoff"],
                         backend=cfg["backend"])
        result["runs"].append({k: r[k] for k in (
            "chi", "times", "electric_energy", "z", "error_estimate", "error_estimate_per_step",
            "max_bond", "wall_time_s")})
        result["runs"].sort(key=lambda x: cfg["chi"].index(x["chi"]))   # keep config order
        out.write_text(json.dumps(result))        # checkpoint after every chi
        print(time.strftime("%H:%M:%S"), f"chi={chi}: err={r['error_estimate']:.3e} "
              f"max_bond={r['max_bond']} ({r['wall_time_s']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
