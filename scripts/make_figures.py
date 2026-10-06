"""Figures (i)-(v) of reports/collaborator_summary.md -> reports/figures/fig*.png.
Reads only committed JSON results. Usage: python3 scripts/make_figures.py"""
import json
import pathlib
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "figures"
COL = {"vacuum": "C0", "neel_half": "C1", "stripe": "C2"}


def load(p):
    return json.loads((ROOT / p).read_text())


def fig1():
    runs = load("results/m3/summary.json")["runs"]
    fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for col, N in enumerate((12, 24)):
        for row, model in enumerate(("su2", "abelian")):
            ax = axs[row, col]
            for r in runs:
                if r["N"] != N or abs(r["g"] - 1.25) > 1e-9 or r["model"] != model:
                    continue
                t = np.array(r["times"]) / r["g"] ** 2
                c = COL[r["state"]]
                ax.plot(t, np.array(r["electric_energy"]) / N, color=c, label=r["state"])
                ax.axhline(r["thermal"]["electric_energy"] / N, color=c, ls="--", lw=1)
                ax.axhline(r["infinite_temperature"]["electric_energy"] / N, color="k", ls=":")
            ax.set_title(f"N = {N}, g = 1.25, {'SU(2)' if model == 'su2' else 'abelian twin'}")
            ax.set_ylabel(r"$\langle H_E\rangle / N$")
            if row == 1:
                ax.set_xlabel(r"$t/g^2$")
    axs[0, 0].legend(title="initial state")
    fig.suptitle("Exact dynamics. Dashed: thermal anchor. Dotted: decohered (infinite temperature)")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_exact_dynamics.png", dpi=120)
    plt.close(fig)


def fig2():
    d = load("results/m4/noise_runs.json")
    cs = d["configs"]
    lab = [f"{c['steps']} steps\np2={c['p2']:g}" for c in cs]
    x = np.arange(len(cs))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for key, name, m, off in (("raw_electric_energy", "raw noisy", "o", -0.2),
                              ("rescaled_electric_energy", "echo-rescaled", "s", 0.0),
                              ("zne_electric_energy", "ZNE (folds 1,3,5)", "^", 0.2)):
        err = [abs(c[key] - c["ideal_electric_energy"]) for c in cs]
        axs[0].bar(x + off, err, width=0.2, label=name)
    axs[0].set_xticks(x, lab)
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"$|\langle H_E\rangle - \langle H_E\rangle_{ideal}|$")
    axs[0].set_title("Error of the electric energy (N = 12, stripe, g = 1.25)")
    axs[0].legend()
    axs[1].plot(x, [c["F_predicted"] for c in cs], "ko-", label=r"predicted $\prod(1-\epsilon)$")
    axs[1].plot(x, [c["F_echo"] for c in cs], "C3s--", label="Loschmidt echo")
    axs[1].set_xticks(x, lab)
    axs[1].set_ylim(0, 1)
    axs[1].set_ylabel("global fidelity F")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "fig2_noise_mitigation.png", dpi=120)
    plt.close(fig)


def fig3():
    d = load("results/m5/resources.json")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, gs in zip(axs, ("native_zz", "cz")):
        for k, t in enumerate(d["tori"]):
            for p, ls in zip(d["p2_values"], ("-", "--", ":")):
                rows = sorted((r for r in d["rows"] if r["N"] == t["N"] and r["gate_set"] == gs
                               and abs(r["p2"] - p) < 1e-15), key=lambda r: r["steps"])
                ax.plot([r["steps"] for r in rows], [r["raw_fidelity"] for r in rows], ls,
                        color=f"C{k}", label=f"N={t['N']}, p2={p:g}" if gs == "native_zz" else None)
        ax.set_yscale("log")
        ax.set_xlabel("first-order Trotter steps")
        ax.set_title(f"{'native ZZ (ions)' if gs == 'native_zz' else 'CZ hardware'}: "
                     r"raw fidelity $e^{-p_2 N_{2q}}$")
    axs[0].set_ylabel("estimated raw circuit fidelity")
    axs[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_resources.png", dpi=120)
    plt.close(fig)


def fig4():
    run = load("results/m6/n48_stripe_g1.25.json")
    conv = load("results/m6/convergence.json")
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    for r in sorted(run["runs"], key=lambda r: r["chi"]):
        axs[0].plot(r["times"], r["electric_energy"], "o-", ms=3, label=f"chi={r['chi']}")
        axs[2].semilogy(range(len(r["error_estimate_per_step"])),
                        np.maximum(r["error_estimate_per_step"], 1e-16), "o-", ms=3,
                        label=f"chi={r['chi']}")
    axs[0].axhline(conv["tol_electric_energy"] * 100, color="k", ls=":", lw=1,
                   label="decohered")
    axs[0].set_xlabel("t"); axs[0].set_ylabel(r"$\langle H_E\rangle$")
    axs[0].set_title("N = 48 torus, stripe, g = 1.25, MPS")
    axs[0].legend(fontsize=8)
    steps = [s["step"] for s in conv["steps"]]
    for i, (a, b) in enumerate(conv["pairs"]):
        axs[1].semilogy(steps, np.maximum([s["d_electric_energy"][i] for s in conv["steps"]],
                                          1e-16), "o-", ms=3, label=rf"$\Delta E$ {a}->{b}")
        axs[1].semilogy(steps, np.maximum([s["d_z_max"][i] for s in conv["steps"]], 1e-16),
                        "s--", ms=3, color=f"C{i}", label=rf"max $\Delta Z$ {a}->{b}")
    axs[1].axhline(conv["tol_electric_energy"], color="k", ls="-", lw=0.8)
    axs[1].axhline(conv["tol_z"], color="k", ls="--", lw=0.8)
    axs[1].set_xlabel("Trotter step"); axs[1].set_title("change between successive chi")
    axs[1].legend(fontsize=7)
    axs[2].set_xlabel("Trotter step"); axs[2].set_ylabel("discarded-weight error estimate")
    axs[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_mps_convergence.png", dpi=120)
    plt.close(fig)


def fig5():
    items = load("results/m6/roadmap.json")["items"]
    fig, ax = plt.subplots(figsize=(12, 1.2 + 1.3 * len(items)))
    ax.axis("off")
    for k, it in enumerate(items):
        y = 1 - (k + 0.5) / len(items)
        ax.text(0.0, y, it["title"], fontsize=11, weight="bold", va="center")
        ax.text(0.27, y, "N = " + ", ".join(map(str, it["N"])), fontsize=10, va="center")
        ax.text(0.40, y, "\n".join(textwrap.wrap(it["purpose"], 85)), fontsize=8.5, va="center")
    ax.set_title("What remains for the GPU cluster and for hardware", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_cluster_roadmap.png", dpi=120)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for f in (fig1, fig2, fig3, fig4, fig5):
        try:
            f()
            print("ok", f.__name__)
        except FileNotFoundError as exc:
            print("skipped", f.__name__, exc)


if __name__ == "__main__":
    main()
