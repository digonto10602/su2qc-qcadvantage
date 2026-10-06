# M2 — two-qubit gate counts of the Trotter circuits

**What was asked.** Count the two-qubit gates that the Trotter circuits of `su2torus/trotter.py` actually use, per plaquette and per time step, for first order (order 1) and second order (order 2), under the two hardware conventions in METHODS §3. Then check the claim of $\approx 9.5$ per plaquette per step.

**Definitions.** A *Trotter step* approximates $e^{-iH\delta t}$ by a product of exponentials of the parts $H_E$ (electric), $H_B^A$ and $H_B^B$ (magnetic, on the two sublattices). An *entangler* is a two-qubit gate. With *native ZZ* (trapped ions), `cx` and `rzz` each cost 1. On *CZ hardware*, `cx` costs 1 and `rzz` costs 2. $n$ is the number of plaquettes, which equals the number of qubits.

**What was done.** `scripts/gate_counts.py` builds each step circuit and counts the `cx` and `rzz` gates with Qiskit's `count_ops`, without transpilation. The raw numbers are in `results/m2/gate_counts.json`.

| lattice | model | order | native ZZ / ($n\cdot$step) | CZ / ($n\cdot$step) |
|---|---|---|---|---|
| torus (any size; checked at $n$ = 12, 24, 48) | SU(2) | 1 | **9.5** | **11.0** |
| torus | SU(2) | 2 (one step on its own) | 15.0 | 18.0 |
| torus | SU(2) | 2 (consecutive steps, with the $U_E$ half-steps merged) | 13.5 | 15.0 |
| torus | abelian | 1 | 1.5 | 3.0 |
| torus | abelian | 2 | 3.0 | 6.0 |
| open $3\times4$ | SU(2) | 1 | 7.04 | 8.25 |

**What it means.**
- The first-order claim holds exactly on the torus. There are $1.5n$ bonds, each with one `rzz`. Each of the $n$ plaquettes has a flip with 3 neighbours, built from $2^3=8$ `cx`. That gives $8n+1.5n=9.5n$ native and $8n+3n=11n$ on CZ hardware. The open lattice uses fewer gates because boundary plaquettes have fewer neighbours.
- Second order costs about 1.4–1.6 times as much per step. The $U_B^A$ flip layer appears twice per step, and only the $U_E$ half-steps of neighbouring steps can merge. Second order pays off only if it allows steps more than about 1.4 times larger at the same accuracy.
- The abelian twin, which is the transverse-field Ising model, needs only the `rzz` layer: $1.5n$.

**Problems / assumptions.**
- No transpiler optimisation was applied. A compiler could merge the two `h` gates and the neighbouring `rz` gates on a qubit, but the two-qubit count is already at the Gray-code minimum of $2^k$ `cx` for a uniformly controlled rotation with $k$ controls.
- The order-2 "merged" row is computed by subtracting one $U_E$ layer per step. The circuits themselves do not merge the half-steps.
