# PLAN — su2qc-qcadvantage, cloud phase (M1–M6), then stop for the GPU cluster

Goal: bring the code for the "SU(2) torus frontier run" (`docs/idea.md`) to the point where everything that can be built and validated without a GPU cluster is done, tested and documented, together with an evidence package for collaborators. Specification: `docs/METHODS.md`.

Acceptance tests for M1–M3 are pre-written and LOCKED (`tests/acceptance/m1..m3`). For M4–M6, the `reviewer` agent writes the acceptance tests from the criteria below at the start of the milestone; you then lock them with `python3 tools/lock.py mN` before writing any implementation.

## M1 — model correctness (locked tests: G1–G5)
Implement `su2torus/lattice.py` and `su2torus/hamiltonian.py` (stubs give the API).

## M2 — exact dynamics, Trotter circuits, gate counts (locked tests: G6–G8)
Implement `su2torus/exact.py` and `su2torus/trotter.py`. Record the true two-qubit counts per plaquette per step (orders 1 and 2, native-ZZ and CZ) in `reports/m2_gate_counts.md`.

## M3 — small-torus physics (locked tests)
Implement `su2torus/observables.py` and `su2torus/thermal.py`. Then write `scripts/run_m3.py`, which produces `results/m3/summary.json`, the open-string file and the figures (METHODS sections 6–7). Run the N=24 jobs as one background batch. Write `reports/m3_physics.md`: what the curves show, which observables are flagged (too close to decohered values), and which initial state and coupling look best for the frontier run.

## M4 — noisy emulation and error mitigation (reviewer writes tests)
Install `qiskit-aer` if missing. Use the N=12 torus (and N=16 = 2×4 if time allows) with Trotter circuits of 2–10 steps and 2-qubit depolarizing noise $p_2\in\{7.9\times10^{-4},10^{-3},1.5\times10^{-3}\}$ (Helios, H2, Heron), plus 1-qubit $3\times10^{-5}$ and readout $10^{-3}$.

Mitigation:
- (a) forward–backward Loschmidt echo to estimate the global fidelity $F$, and rescale $\langle O\rangle_\text{mit}=O_\infty+(\langle O\rangle_\text{noisy}-O_\infty)/F$
- (b) zero-noise extrapolation with gate folding at factors 1, 3, 5

Acceptance criteria:
- noiseless emulation equals the statevector to $10^{-10}$
- the echo-estimated $F$ agrees within 25% with the fidelity predicted from the noise model (state the formula used, e.g. $\prod_\text{gates}(1-\epsilon_\text{gate})$ with Aer's depolarizing-parameter-to-infidelity conversion)
- mitigated $\langle H_E\rangle$ error is $\le 0.5\times$ raw error in at least 3 of 4 tested (steps, $p_2$) configurations
- results are seeded and reproducible

Outputs: `results/m4/*.json`, `reports/m4_noise.md`.

## M5 — classical-campaign code ready for the cluster (reviewer writes tests)
Install `quimb` if missing.
- (a) `su2torus/mps.py`: simulate the *same* Trotter circuits as an MPS (quimb `CircuitMPS`, snake ordering of plaquettes, max bond $\chi$ configurable, optional GPU backend flag).
- (b) `scripts/resource_table.py`: compile circuits for tori N = 48 (4×6), 56 (4×7), 64 (4×8), 96 (6×8). Report two-qubit gates per step, totals for 4–10 steps, and estimated raw fidelity $e^{-p_2N_{2q}}$ for the three $p_2$ values, in `results/m5/resources.json`.
- (c) Config-driven entry point `scripts/run_classical.py --config configs/*.json` for MPS runs, plus Slurm templates `scripts/slurm/*.sbatch` with placeholders (`<PARTITION>`, `<ACCOUNT>`, `<ENV>`). Never hard-code real cluster names.

Acceptance criteria:
- MPS matches the exact statevector at N=12 and N=24 for 4 steps at $\chi=256$ (state fidelity $>1-10^{-6}$)
- the discarded-weight/error estimate is reported and decreases with $\chi$
- the resource table matches M2 per-plaquette counts × N × steps
- the Slurm templates contain only placeholders

## M6 — convergence evidence and collaborator package (reviewer writes tests)
Run MPS on the N=48 torus (4×6), stripe state, $g=1.25$, 1–8 first-order steps with $\delta t$ chosen so that $2\delta t/g^2\approx0.9$, at $\chi\in\{32,64,128,256\}$ (cloud limits: each run $\le2$ h, $\le12$ GB). Record per step where observables stop converging in $\chi$.

Write `reports/collaborator_summary.md` (plain English, LaTeX math, every term defined) with figures:
- (i) exact 12/24-plaquette dynamics with thermal and decohered lines
- (ii) noisy-emulation raw versus mitigated
- (iii) the gate/fidelity resource table
- (iv) the MPS convergence breakdown at N=48
- (v) what remains for the GPU cluster (large-$\chi$ MPS, PEPS-BP, Pauli propagation, sign-free QMC thermal anchor at N=48–64) and the hardware runs

Acceptance criteria:
- every figure is referenced and its data file exists
- the numbers quoted in the summary are reproduced from the JSON files by `scripts/check_summary.py`

Update `README.md` with how to reproduce everything. Tag `v0.1-pre-cluster` and push.

## STOP after M6. Do not start cluster work.
