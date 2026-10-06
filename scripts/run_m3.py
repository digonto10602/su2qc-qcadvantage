"""M3: small-torus physics (docs/METHODS.md sections 4-7).

Parts (each writes results/m3/parts/<name>.json; `assemble` builds summary.json + figures):
    python3 scripts/run_m3.py n12                 # all N=12 runs (exact diag thermal)
    python3 scripts/run_m3.py n24 STATE MODEL     # one N=24 run (Chebyshev + typicality)
    python3 scripts/run_m3.py string              # open 3x4 string, g=1.4
    python3 scripts/run_m3.py assemble            # summary.json, flags, figures
"""
import json
import pathlib
import sys
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from su2torus.lattice import HoneycombTorus, initial_state  # noqa: E402
from su2torus import exact, hamiltonian as hm, observables as ob, thermal  # noqa: E402

RES = ROOT / "results" / "m3"
PARTS = RES / "parts"
TORI = {12: (3, 2), 24: (4, 3)}
STATES = ("vacuum", "neel_half", "stripe")
MODELS = ("su2", "abelian")
FLAG_REL = 0.10          # flag if late-time value within 10% of the decohered value
LATE_FRAC = 0.25         # late time = last 25% of the time grid


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def time_grid(g, npts):
    return np.linspace(0.0, 4.0 * g ** 2, npts)


def measure(lat, g, psi, abelian, hexagon, Hpsi=None):
    """Observables of one state. Hpsi = H psi if already available."""
    if Hpsi is None:
        Hpsi = hm.apply_h(lat, g, psi, abelian)
    links = ob.link_energies(lat, psi)
    ee = 0.5 * g ** 2 * links.sum()
    z = ob.z_expectations(lat, psi)
    p, q = lat.bonds[0]
    zz = ob.zz_expectations(lat, psi, [(p, q)])[0] - z[p] * z[q]
    return {"electric_energy": float(ee),
            "magnetic_energy": float(np.vdot(psi, Hpsi).real - ee),
            "link_energies": links.tolist(), "z": z.tolist(),
            "zz_connected_bond0": float(zz),
            "hexagon_x_string": ob.hexagon_x_string(lat, psi, hexagon)}


def run_dynamics(lat, g, abelian, excited, times, H=None, bounds=None):
    """Evolve and measure. H: sparse matrix (small n) or None (matrix-free Chebyshev)."""
    hexagon = lat.hexagons()[0] if lat.hexagons() else []
    psi0 = exact.product_state(lat, excited)
    if H is None:
        dim = 2 ** lat.n
        mv = lambda v: hm.apply_h(lat, g, v, abelian)
        H = LinearOperator((dim, dim), matvec=mv, rmatvec=mv, dtype=complex)
        H.fused = hm.cheb_recur(lat, g, abelian)        # parallel fused Chebyshev step
    out = {k: [] for k in ("electric_energy", "magnetic_energy", "link_energies", "z",
                           "zz_connected_bond0", "hexagon_x_string")}
    for k, psi in enumerate(exact.evolve_iter(H, psi0, times, bounds)):
        m = measure(lat, g, psi, abelian, hexagon, H @ psi if not isinstance(
            H, LinearOperator) else None)
        for key, v in m.items():
            out[key].append(v)
        if lat.n > 20:
            log(f"  t={times[k]:.3f} EE={m['electric_energy']:.6f} EB={m['magnetic_energy']:.6f}")
    out["times"] = list(map(float, times))
    return out


def base_record(N, g, state, model, lat):
    return {"N": N, "g": g, "state": state, "model": model, "L1": lat.L1, "L2": lat.L2,
            "infinite_temperature": ob.infinite_temperature(lat, g, model == "abelian")}


def part_n12():
    N = 12
    lat = HoneycombTorus(*TORI[N])
    runs = []
    for g in (1.1, 1.25, 1.4):
        for model in MODELS:
            ab = model == "abelian"
            H = hm.build_sparse(lat, g, ab)
            HE, _, _ = hm.split_sparse(lat, g, ab)
            for state in STATES:
                exc = initial_state(lat, state)
                rec = base_record(N, g, state, model, lat)
                rec.update(run_dynamics(lat, g, ab, exc, time_grid(g, 41), H=H))
                E0 = rec["electric_energy"][0] + rec["magnetic_energy"][0]
                beta = thermal.beta_for_energy_exact(H, E0)
                th = thermal.thermal_exact(H, beta, {"electric_energy": HE})
                rec["thermal"] = {"beta": beta, "electric_energy": th["electric_energy"],
                                  "magnetic_energy": th["energy"] - th["electric_energy"],
                                  "method": "full diagonalization"}
                runs.append(rec)
                log(f"N=12 g={g} {model} {state}: beta={beta:.4f}")
    return runs


def part_n24(state, model, beta0):
    N, g = 24, 1.25
    lat = HoneycombTorus(*TORI[N])
    ab = model == "abelian"
    dim = 2 ** lat.n
    mv = lambda v: hm.apply_h(lat, g, v, ab)
    op = LinearOperator((dim, dim), matvec=mv, rmatvec=mv, dtype=complex)
    t0 = time.time()
    bf = PARTS / f"bounds_n24_{model}.json"        # Lanczos bounds, shared by all states
    if bf.exists():
        bounds = tuple(json.loads(bf.read_text()))
    else:
        bounds = exact.spectral_bounds(op)
        bf.write_text(json.dumps(bounds))
    log(f"N=24 {model} {state}: spectral bounds {bounds} ({time.time() - t0:.0f}s)")
    rec = base_record(N, g, state, model, lat)
    rec.update(run_dynamics(lat, g, ab, initial_state(lat, state), time_grid(g, 25),
                            bounds=bounds))
    log(f"dynamics done ({time.time() - t0:.0f}s)")
    E0 = rec["electric_energy"][0] + rec["magnetic_energy"][0]
    diag = hm.electric_diagonal(lat, g)
    beta, th, hist = thermal.beta_for_energy_typicality(
        mv, dim, E0, {"electric_energy": lambda v: diag * v}, beta0=beta0, step=0.05,
        n_samples=4, seed=1234, max_evals=8, bounds=bounds,
        recur=hm.cheb_recur(lat, g, ab))
    rec["thermal"] = {"beta": beta, "electric_energy": th["electric_energy"],
                      "magnetic_energy": th["energy"] - th["electric_energy"],
                      "method": "typicality, 4 vectors, seed 1234, Chebyshev exp(-beta H/2)",
                      "search": [[b, r["energy"]] for b, r in hist]}
    log(f"thermal beta={beta:.4f} E={th['energy']:.5f} target={E0:.5f} "
        f"evals={len(hist)} ({time.time() - t0:.0f}s)")
    rec.pop("z")          # optional for N=24; keeps summary.json small
    return rec


def part_string():
    g = 1.4
    lat = HoneycombTorus(3, 4, periodic=False)
    exc = initial_state(lat, "string")
    rec = {"L1": 3, "L2": 4, "periodic": False, "g": g, "state": "string", "model": "su2",
           "excited": exc, "bonds": [list(b) for b in lat.bonds]}
    rec.update(run_dynamics(lat, g, False, exc, time_grid(g, 41)))
    return rec


# ------------------------------------------------------------------ assemble
def flags(rec):
    """Observables whose late-time mean lies within FLAG_REL of the decohered value.
    Decohered value 0 (magnetic energy, z, zz, hexagon string): relative to the largest
    |value| along the run."""
    t = np.array(rec["times"])
    late = t >= t[-1] * (1 - LATE_FRAC)
    inf = rec["infinite_temperature"]
    series = {"electric_energy": (np.array(rec["electric_energy"]), inf["electric_energy"]),
              "magnetic_energy": (np.array(rec["magnetic_energy"]), inf["magnetic_energy"]),
              "link_energy_mean": (np.array(rec["link_energies"]).mean(axis=1),
                                   inf["link_energy"]),
              "zz_connected_bond0": (np.array(rec["zz_connected_bond0"]), 0.0),
              "hexagon_x_string": (np.array(rec["hexagon_x_string"]), 0.0)}
    if "z" in rec:
        series["z_mean"] = (np.array(rec["z"]).mean(axis=1), 0.0)
    out = []
    for name, (x, ref) in series.items():
        scale = abs(ref) if ref != 0 else np.abs(x).max()
        if scale > 0 and abs(x[late].mean() - ref) <= FLAG_REL * scale:
            out.append(name)
    return out


def assemble():
    runs = json.loads((PARTS / "n12.json").read_text())
    for s in STATES:
        for m in MODELS:
            f = PARTS / f"n24_{s}_{m}.json"
            if f.exists():
                runs.append(json.loads(f.read_text()))
            else:
                log("missing", f)
    for r in runs:
        r["flags"] = flags(r)
    (RES / "summary.json").write_text(json.dumps({"runs": runs,
        "flag_rule": f"late-time (last {LATE_FRAC:.0%}) mean within {FLAG_REL:.0%} of the "
                     "decohered value (scale = max|O(t)| when that value is 0)"}))
    s = json.loads((PARTS / "string_open_3x4_g1.4.json").read_text())
    (RES / "string_open_3x4_g1.4.json").write_text(json.dumps(s))
    figures(runs, s)
    log("assembled", len(runs), "runs")


def figures(runs, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = RES / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    colors = {"vacuum": "C0", "neel_half": "C1", "stripe": "C2"}

    def get(N, g, st, m):
        return next(r for r in runs if r["N"] == N and abs(r["g"] - g) < 1e-9
                    and r["state"] == st and r["model"] == m)

    def energy_panel(ax, N, g, model, key):
        for st in STATES:
            try:
                r = get(N, g, st, model)
            except StopIteration:
                continue
            t = np.array(r["times"]) / g ** 2
            ax.plot(t, np.array(r[key]) / N, color=colors[st], label=st)
            ax.axhline(r["thermal"][key] / N, color=colors[st], ls="--", lw=0.8)
            ax.axhline(r["infinite_temperature"][key] / N, color="k", ls=":", lw=0.8)
        ax.set_xlabel("$t/g^2$")
        ax.set_ylabel(f"{key.replace('_', ' ')} / N")
        ax.set_title(f"N={N}, g={g}, {model}")

    k = 0
    for N, gs in ((12, (1.1, 1.25, 1.4)), (24, (1.25,))):
        for g in gs:
            fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
            for j, model in enumerate(MODELS):
                energy_panel(axs[0, j], N, g, model, "electric_energy")
                energy_panel(axs[1, j], N, g, model, "magnetic_energy")
            axs[0, 0].legend(fontsize=8)
            fig.suptitle("dashed: thermal anchor, dotted: decohered (infinite T)")
            fig.tight_layout()
            fig.savefig(figdir / f"energies_N{N}_g{g}.png", dpi=110)
            plt.close(fig)
            k += 1
    # finite-size comparison N=12 vs N=24 (per plaquette), g=1.25
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for j, model in enumerate(MODELS):
        for st in STATES:
            for N, ls in ((12, ":"), (24, "-")):
                try:
                    r = get(N, 1.25, st, model)
                except StopIteration:
                    continue
                axs[j].plot(np.array(r["times"]) / 1.25 ** 2,
                            np.array(r["electric_energy"]) / N, ls=ls, color=colors[st],
                            label=f"{st} N={N}")
        axs[j].set_title(f"finite size, g=1.25, {model}")
        axs[j].set_xlabel("$t/g^2$"); axs[j].set_ylabel("electric energy / N")
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(figdir / "finite_size_g1.25.png", dpi=110); plt.close(fig)
    # hexagon string and zz for N=12 g=1.25
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for model, ls in (("su2", "-"), ("abelian", "--")):
        for st in STATES:
            r = get(12, 1.25, st, model)
            t = np.array(r["times"]) / 1.25 ** 2
            axs[0].plot(t, r["hexagon_x_string"], ls=ls, color=colors[st], label=f"{st} {model}")
            axs[1].plot(t, r["zz_connected_bond0"], ls=ls, color=colors[st])
    axs[0].set_title("hexagon X string, N=12, g=1.25"); axs[1].set_title("connected ZZ, bond 0")
    for a in axs:
        a.set_xlabel("$t/g^2$"); a.axhline(0, color="k", ls=":", lw=0.8)
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(figdir / "strings_zz_N12_g1.25.png", dpi=110); plt.close(fig)
    # open string: z heat map and energies
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    t = np.array(s["times"])
    im = axs[0].imshow(np.array(s["z"]).T, aspect="auto", origin="lower", cmap="RdBu",
                       vmin=-1, vmax=1, extent=[0, t[-1], -0.5, len(s["z"][0]) - 0.5])
    fig.colorbar(im, ax=axs[0], label="<Z_p>")
    axs[0].set_xlabel("t"); axs[0].set_ylabel("plaquette p"); axs[0].set_title("open 3x4 string, g=1.4")
    axs[1].plot(t, s["electric_energy"], label="electric")
    axs[1].plot(t, s["magnetic_energy"], label="magnetic")
    axs[1].plot(t, np.array(s["electric_energy"]) + np.array(s["magnetic_energy"]), "k:", label="total")
    axs[1].set_xlabel("t"); axs[1].legend()
    fig.tight_layout(); fig.savefig(figdir / "string_open_3x4_g1.4.png", dpi=110); plt.close(fig)


def main():
    PARTS.mkdir(parents=True, exist_ok=True)
    what = sys.argv[1]
    if what == "n12":
        (PARTS / "n12.json").write_text(json.dumps(part_n12()))
    elif what == "n24":
        state, model = sys.argv[2], sys.argv[3]
        n12 = json.loads((PARTS / "n12.json").read_text())
        beta0 = next(r["thermal"]["beta"] for r in n12 if r["state"] == state
                     and r["model"] == model and abs(r["g"] - 1.25) < 1e-9)
        rec = part_n24(state, model, beta0)
        (PARTS / f"n24_{state}_{model}.json").write_text(json.dumps(rec))
    elif what == "string":
        (PARTS / "string_open_3x4_g1.4.json").write_text(json.dumps(part_string()))
    elif what == "assemble":
        assemble()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
