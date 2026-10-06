# One SU(2) idea with a credible path to quantum advantage on today's hardware

*Prepared 2 Oct 2026 for the project "Funky Ideas about QC advantage". Every paper cited below was checked against its arXiv/INSPIRE/journal record during preparation; estimates that are mine are labelled "[my estimate]".*

## 1. What was asked

One idea related to SU(2) gauge theory that can show clear quantum advantage, with a concrete step-by-step plan, feasible on quantum hardware available today, with every step connected to real published work and no invented physics.

## 2. What was done

Three parallel literature sweeps (non-Abelian lattice-gauge-theory hardware runs 2019–2026; the 2023–2026 "beyond-classical" simulation claims and the classical methods that answered them; the technical details of the minimal SU(2) truncation), followed by direct verification of ~25 primary sources, and a resource analysis (gate counts, circuit fidelity, classical cost arguments) that is labelled as mine wherever it is not taken from a paper.

## 3. The idea (my label: "the SU(2) torus frontier run")

Take the cheapest published, exactly gauge-invariant formulation of 2+1-dimensional SU(2) Yang–Mills theory — the one-qubit-per-plaquette, minimally truncated Hamiltonian that Ciavarella, de Putter, Younis and Rrapaj ran on all 156 qubits of IBM's `ibm_boston` in August 2026 [1] — and run it **not** as the string-breaking experiment they did (which tensor networks reproduced), but as a **thermalizing quench on a periodic lattice (a torus) of 48–64 plaquettes on Quantinuum's all-to-all trapped-ion machines (Helios, 98 qubits; H2, 56 qubits), for 7–10 Trotter steps, at the coupling where electric and magnetic energies compete**, with the observables, classical-benchmarking campaign and verification ladder that the 2026 beyond-classical papers established [5–8]. The deliverable is the first 2+1D non-Abelian entry in the community's Quantum Advantage Tracker [9], together with new physics: thermalization and energy transport in the minimal SU(2) theory compared with its Abelian twin.

### Jargon used below (defined once)

- **Lattice gauge theory (LGT), Hamiltonian/Kogut–Susskind form:** space is a lattice, time is continuous; gauge fields live on links; the energy has an *electric* part (a cost per unit of flux on a link) and a *magnetic* part (an operator that flips flux around a *plaquette*, the smallest closed loop). **Gauss's law** is the local constraint that flux lines cannot end in empty space.
- **Truncation, $j_{\max}=1/2$:** each SU(2) link carries a "spin" label $j = 0, 1/2, 1, \dots$; keeping only $j\le 1/2$ is the minimal truncation. **Local Krylov truncation** [2] builds the kept states by acting with plaquette operators on the empty vacuum, which fixes Gauss's law automatically.
- **Trotter step:** one discrete time slice of the digitized evolution $e^{-iH\delta t}$.
- **Two-qubit gate:** the entangling operation whose error dominates on today's machines.
- **Tensor network (TN):** compressed classical representation of the quantum state; **MPS** (chain-shaped), **PEPS** (grid-shaped); accuracy set by the **bond dimension** $\chi$ (bigger $\chi$ = more entanglement captured). **Belief propagation (BP):** an approximate way to contract 2D TNs that is exact on trees and degrades on loops. **Pauli propagation / sparse Pauli dynamics:** classical methods that evolve the *observable* as a sum of Pauli strings and discard small terms; they thrive on near-Clifford circuits.
- **Torus / periodic boundary conditions:** the lattice wraps around in both directions, so there is no boundary. **Min-cut:** the minimum number of bonds a classical MPS has to cut to split the system in two; it sets the entanglement an MPS must carry.
- **Stoquastic:** a Hamiltonian whose off-diagonal elements are all $\le 0$ in some basis, so its *thermal* properties can be computed by classical quantum Monte Carlo (QMC) without a sign problem.
- **ETH:** eigenstate thermalization hypothesis — the statement that a closed chaotic quantum system relaxes to thermal values.
- **ZNE / ODR / TREX / PEC:** error-mitigation techniques (zero-noise extrapolation, operator decoherence renormalization, readout twirling, probabilistic error cancellation) that reduce but do not correct hardware noise.
- **Loschmidt echo:** evolve forward, then backward; a perfect device returns exactly to the start, so the return probability measures the circuit's global fidelity.

### 3.1 The Hamiltonian (verbatim from [1])

Pure SU(2) on a **triangular** lattice, minimal local-Krylov truncation, one qubit per triangular plaquette $p$, $\hat P^{0,1}_p$ projectors onto the plaquette being unexcited/excited, $\hat n$ running over the **three** plaquettes that share a link with $p$:

$$
\hat H_E=\frac{3g^2}{8}\sum_p\Big(3\hat P^1_p-\sum_{\hat n}\hat P^1_p\hat P^1_{p+\hat n}\Big),\qquad
\hat H_B=-\frac{1}{g^2}\sum_p\Big(\prod_{\hat n}\hat C_{p+\hat n}\Big)\hat X_p,\qquad
\hat C_p=\hat P^0_p+\tfrac12\hat P^1_p ,
$$

with link electric energy $\hat E^2_l=\tfrac38\,(1-\hat Z_{p_1}\hat Z_{p_2})$ for the two plaquettes sharing link $l$.

Rewritten in Pauli operators [my derivation, using $\hat P^1=(1-\hat Z)/2$ and $\hat C=(3+\hat Z)/4$]:

$$
\hat H=\frac{3g^2}{16}\sum_{\langle pq\rangle}\big(1-\hat Z_p\hat Z_q\big)\;-\;\frac{1}{g^2}\sum_p \hat X_p\prod_{\hat n}\frac{3+\hat Z_{p+\hat n}}{4}.
$$

Reading: the plaquette qubits live on the *dual* graph of the triangular lattice, which is a **honeycomb** graph (every plaquette has exactly three link-sharing neighbours). The electric term is a ferromagnetic Ising coupling (a flux line is a domain wall); the magnetic term flips a plaquette with amplitude $1/g^2$ **halved for every excited link-neighbour** — that factor of $\tfrac12$ per neighbour is the SU(2) recoupling and is the *only* place where the non-Abelian structure survives at this truncation. Setting $\hat C\to 1$ gives the transverse-field Ising model on the honeycomb graph, which is the Wegner dual of the $\mathbb Z_2$ gauge theory on the same triangular lattice [my derivation]. So the SU(2) model and its Abelian twin can be run on identical qubits, at identical depth, differing only in those factors of $\tfrac12$.

Two properties matter for the plan: (i) the model is **stoquastic** in the $Z$ basis (off-diagonal elements $-(1/g^2)(1/2)^m\le0$), so its *thermal* expectation values are classically computable by QMC at any size even though its *real-time* dynamics are not [my observation]; (ii) this exact model class is already known to be chaotic/thermalizing: the honeycomb $j_{\max}=1/2$ SU(2) theory obeys ETH numerically [11], plaquette chains show an emergent energy-diffusion mode [12], and local thermalization of 151-plaquette chains has been measured on IBM hardware [13].

### 3.2 Gate cost per Trotter step [my estimate, consistent with [1,2]]

The flip term is an $X$-rotation of plaquette $p$ whose angle depends on three control qubits (angles $\theta,\theta/2,\theta/4,\theta/8$). A rotation controlled by $k$ qubits compiles into $2^k$ two-qubit gates (the standard "uniformly controlled rotation" construction); for four neighbours on a square lattice Ciavarella–Burbano–Bauer report exactly $16=2^4$ CNOT per plaquette per step [2], so for three neighbours the count is $8$. The electric term costs one native arbitrary-angle $ZZ$ gate per bond, $1.5$ per plaquette. **Total $\approx 9.5$ two-qubit gates per plaquette per first-order step** on all-to-all hardware with native $ZZ$ gates; about $12$ on IBM's heavy-hex with routing (Ciavarella's dense encoding: 7,634 CZ for a 128-plaquette lattice [1]). For comparison, a transverse-field Ising step costs $\approx1.5$–$2$ gates per spin — the non-Abelian recoupling makes each step roughly $5\times$ more expensive, which is the quantitative reason non-Abelian advantage lags Abelian spin models.

### 3.3 Why this is the right target, and why the alternatives are not

- **1+1D SU(2) is lost to MPS.** The 120-qubit SU(2) loop-string-hadron experiment on `ibm_boston` [10] is listed as "active" on the Tracker, but Roland Farrell reproduced it with an MPS of bond dimension 800 in 3,003 s on a laptop and pointed out that its headline observable (a conserved charge) would be exact even for a fully decohered device [9].
- **Open-boundary lattices on 156-qubit IBM chips stay within MPS reach** [my estimate]: the shortest side is $\le 12$–$13$ plaquettes, so the MPS min-cut is $\le13$ bonds and $\chi\sim 10^4$ captures even a volume-law state; this is why the 16×8 run of [1] converged at $\chi=540$. A torus doubles the cut: Quantinuum's 56-qubit 7×8 periodic Ising lattice (min-cut 14) is where $\chi=4000$ MPS "is too low to even quantify how wrong the extrapolation might be" beyond 9 Trotter steps [5].
- **The exact honeycomb $j_{\max}=1/2$ formulation** (Müller–Yao [14]; Turro–Ciavarella–Yao [15]) is physically cleaner (its plaquette amplitude depends on six external legs, $(-\tfrac12)^{c}$) but costs up to 70 CNOT per plaquette per step [15] — only 1–2 steps at frontier size. It is the right *small-scale comparator*, not the main run.
- **Provable (BQP-complete) SU(2) quantities exist** — Wilson-loop expectation values in SU(2)$_k$ Chern–Simons theory are Jones polynomials — but Quantinuum's own end-to-end study on H2-2 concludes advantage needs $\sim100$ qubits at two-qubit error $\le10^{-4}$ and $\gtrsim2{,}800$ crossings [16]; not today.
- **Thermal / finite-density SU(2)** would attack the sign problem directly but needs thermal-state preparation at scale; not today.
- No published theorem states that real-time dynamics of a Kogut–Susskind gauge theory is BQP-complete, so "clear advantage" for *any* gauge-theory dynamics today means the empirical standard now used by Quantinuum [5], IBM/Qedma [6] and IBM/Algorithmiq [7]: a regime where every known classical method is shown to be either unconverged or mutually inconsistent, with the quantum result validated by an independent trust ladder and posted publicly for classical attack [9].

### 3.4 Hardware fit (vendor numbers)

| Machine | Qubits / connectivity | Two-qubit infidelity | Demonstrated circuit volume |
|---|---|---|---|
| Quantinuum Helios | 98, all-to-all [17] | $7.9(2)\times10^{-4}$ [17] | 3,439 two-qubit gates at up to 90 qubits [18] |
| Quantinuum H2 | 56, all-to-all | $\sim1\times10^{-3}$ (data sheet) | $>2{,}000$ two-qubit gates, 20 Trotter steps [5] |
| IBM Heron r3 (`ibm_boston`) | 156, heavy-hex | $\sim1.5\times10^{-3}$ median (third-party) | 7,634 CZ, depth 218, this Hamiltonian [1] |
| IBM Nighthawk r2 | 120, square lattice | comparable to Heron (IBM) | "accurate estimates on 7,500+ gate circuits" (IBM) |

All-to-all connectivity is what makes the torus free; IBM machines serve as the cross-platform check (open boundaries).

### 3.5 Resource table for the main run [my estimates]

Honeycomb torus with $L_1\times L_2$ unit cells, $N=2L_1L_2$ plaquettes, $3N/2$ bonds, $9.5N$ two-qubit gates per first-order step:

| $N$ (cells) | gates / step | steps at 3,650 gates | steps at 4,250 gates | raw circuit fidelity on Helios |
|---|---|---|---|---|
| 48 (4×6) | 456 | 8 | 9 | 5.6% / 3.5% |
| 56 (4×7) | 532 | 7 | 8 | 5.6% / 3.5% |
| 64 (4×8) | 608 | 6 | 7 | 5.6% / 3.5% |
| 96 (6×8) | 912 | 4 | 4–5 | (future: needs $\sim2\times$ lower error) |

Raw fidelity $=e^{-\epsilon_{2q}\times\text{gates}}$; local observables see only their light cone and are attenuated less. Choose $\delta t$ so the largest flip angle per step, $2\delta t/g^2$, is $\approx0.8$–$1.0$ rad (Quantinuum used a field angle of $1.0$ rad per step [5]); 7–9 steps then reach $t/g^2\approx3$–$4$, past the $tJ\approx2$ point beyond which the independent comparative study of Vovrosh et al. finds that all large-scale classical methods disagree for near-critical 2D quenches with $L\ge10$ [8].

## 4. Step-by-step plan (6–9 months)

1. **Build and validate the model (weeks 1–3).** Implement $\hat H$ on honeycomb tori; derive the Pauli form; compile one Trotter step with uniformly controlled rotations and native $ZZ$ gates; confirm the gate count. Validate against exact diagonalization on 12- and 24-plaquette tori, and reproduce two published anchors: Ciavarella et al.'s string-breaking curves on small open lattices at $g=1.4$ [1], and the transverse-field Ising limit $\hat C\to1$. Cross-check the truncation against the exact honeycomb $j_{\max}=1/2$ model [14,15] and against $j_{\max}=1$ on $\le 4$ plaquettes to quantify what the minimal truncation misses at the chosen coupling (the published leakage bound is $0.095$ at $g=1.4$ [1]; at smaller $g$ it must be measured).
2. **Choose coupling, initial state and observables classically (weeks 3–6).** Locate the model's order–disorder crossover with DMRG/QMC (sign-free), and pick $g$ near it (bare ratio of magnetic to electric scales is $16/(3g^4)$; $g\approx1.1$–$1.3$ is the expected window [my estimate]) plus a strong-coupling control at $g=1.4$ where the truncation is validated. Pick initial product states in the plaquette basis (always gauge-invariant) of *intermediate* energy density — e.g. stripes, or two halves of the torus at different flux densities for transport — such that the QMC thermal values differ clearly from the infinite-temperature (fully decohered) values; this is the Farrell criterion [9]. Observables: link electric energy $\langle\hat E^2_l\rangle$ and its spatial profile (energy transport), connected $\langle\hat Z_p\hat Z_q\rangle$ correlators versus distance (light-cone spreading), the magnetic energy $\langle\hat X_p\prod\hat C\rangle$ (measured in the $X$ basis), total energy (conserved — a built-in Trotter+noise monitor), closed $X$-strings around hexagons of the plaquette graph (high weight, sensitive to global fidelity), and Rényi-2 entropy of 2–3 plaquettes via randomized measurements as in [13]. Define one scalar with error bars as the Tracker instance.
3. **Calibrate on emulator and small hardware (weeks 6–10).** Run the full-depth circuits on the Quantinuum emulator with its noise model, then on hardware for the 24-plaquette torus, where exact diagonalization gives ground truth at the *same depth* as the frontier run. Fix the mitigation recipe here (ZNE/ODR-style rescaling, readout twirling) and measure the global fidelity with the forward–backward Loschmidt echo — the Lewis-group "self-mitigation" that was introduced precisely for SU(2) Trotter circuits [19].
4. **Frontier runs (weeks 10–16).** $N=48$ and $56$ tori, 6–10 steps, two couplings, 2–3 initial states; the Abelian twin at identical depth and also to 20–30 steps (it is $\sim5\times$ cheaper); the unperturbed echo at full depth for every configuration; shot budgets of order $10^3$ per circuit. Cross-platform: the same Hamiltonian on IBM Heron in the geometric encoding of [1] (open 8×8, ancilla post-selection) on shared observables at matched depth.
5. **Classical campaign, in parallel (weeks 6–20).** MPS with $\chi$ up to $\sim10^4$ on GPUs; PEPS with belief propagation and full update [20]; PEPO; Pauli propagation [21] and sparse Pauli dynamics [22]; optionally neural quantum states. Record, for each method, the step at which it stops converging in its own control parameter and the step at which methods disagree with each other — exactly the evidence structure of [5–7]. Use QMC thermal values as the late-time anchor and exact/TN results as the early-time anchor; the quantum data must interpolate between them ("sandwich" check — my framing, modelled on the opposite-limits argument of [7]).
6. **Trust ladder and claim (weeks 20–26).** Apply the validation criteria of [7]: consistency under controlled noise changes, independence of the attenuation from the observable/perturbation, cross-device reproducibility, convergence in mitigation parameters, agreement with exact results at 24 plaquettes, energy conservation. Publish circuits, data, error bars and all classical baselines; submit as an "observable estimation" instance to the Quantum Advantage Tracker [9]; write the paper.
7. **Scale-up (when hardware allows).** The same circuits on a 96-plaquette Helios torus need $\sim8{,}000$ two-qubit gates for 9 steps — i.e. a two-qubit error near $4$–$5\times10^{-4}$ or an IBM Nighthawk-class gate budget with a folded-torus layout. At that size the MPS min-cut is $\approx16$–$24$ bonds and the claim becomes uncontestable for MPS, with PEPS/Pauli methods already failing at $\ge9$ steps.

## 5. What the result would mean

- A 2+1D non-Abelian gauge theory simulated in the regime where, by the community's current standard, "no known classical methods are both efficient and trustworthy" [5] — the first such entry for a non-Abelian theory on the Tracker, and the first SU(2) run on a torus at this scale.
- Physics: whether and how the minimal SU(2) theory thermalizes and transports energy, and how the factor-$\tfrac12$ recoupling changes relaxation relative to the $\mathbb Z_2$ twin — extending the ETH/hydrodynamics programme of [11–13] into 2D at scale.
- Methodology: a reusable, published SU(2)-specific verification kit (exact-depth small-torus anchor, sign-free thermal anchor, echo fidelity, Abelian twin) that later $j_{\max}=1$ or matter-coupled runs can inherit.

## 6. Problems and assumptions (read before committing)

1. **"Clear" means frontier-standard, not proven.** With 48–64 plaquettes and 7–9 steps the run sits at the same scale as the Quantinuum 2026 claim; an MPS with $\chi\sim10^5$ on a national-lab machine might still reproduce it, in which case the Tracker entry would become "superseded" — the honest risk every current claim carries. The 96-plaquette version (step 7) is the one that escapes MPS by construction.
2. **The claim is about the truncated Hamiltonian.** $j_{\max}=1/2$ is far from the continuum; the model's order–disorder crossover is a truncation artifact (continuum 2+1D SU(2) is always confining), and the regime that is classically hardest (strong magnetic term) is where the truncation is least faithful. State the result as "beyond-classical dynamics of a non-Abelian gauge-invariant lattice Hamiltonian," and quantify the truncation effect at small volume as in step 1.
3. **Gate counts and fidelities are estimates** ($8+1.5$ per plaquette-step; raw fidelities 3–6%). They are consistent with [1,2] but must be confirmed by compilation on the target machine; a $1.5\times$ miss removes 2–3 steps.
4. **Large $\delta t$ means large Trotter error.** The classical comparison is to the *same circuit*, so this does not weaken the advantage claim, but the physics claim about continuous-time SU(2) needs the $\delta t$-dependence study at small size.
5. **Observables can be fooled by decoherence** if their thermal value equals the infinite-temperature value; step 2 exists to prevent this.
6. **Access and cost.** Helios is commercial; H2 and IBM machines are available through DOE's QCUP (IBM, Quantinuum $N\ge56$, IonQ, IQM) — the ORNL application already in progress is the natural vehicle. Deep all-to-all trapped-ion circuits run at low shot rates, so the shot plan must be budgeted.
7. **Hilbert-space fragmentation** is generic in Kogut–Susskind truncations with hard constraints [23]; the SU(2) factors here are soft ($\tfrac12$, not $0$), so no exact fragmentation is expected, but slow prethermal dynamics at large $g$ would make the system easier classically — another reason to work near the crossover.

## 7. References (all verified during preparation)

1. A. N. Ciavarella, R. de Putter, E. Younis, E. Rrapaj, "Quantum simulations of two-dimensional non-Abelian adjoint string breaking," arXiv:2608.28752 (Aug 2026).
2. A. N. Ciavarella, I. M. Burbano, C. W. Bauer, "Efficient truncations of SU($N_c$) lattice gauge theory for quantum simulation," Phys. Rev. D 112, 054514 (2025), arXiv:2503.11888.
3. A. N. Ciavarella, C. W. Bauer, "Quantum simulation of SU(3) lattice Yang–Mills theory at leading order in large-$N_c$ expansion," Phys. Rev. Lett. 133, 111901 (2024), arXiv:2402.10265.
4. T. A. Cochran et al. (Google Quantum AI), "Visualizing dynamics of charges and strings in (2+1)D lattice gauge theories," Nature 642, 315 (2025), arXiv:2409.17142.
5. R. Haghshenas et al. (Quantinuum), "Digital quantum magnetism on a trapped-ion quantum computer," Nature 653, 56 (2026), arXiv:2503.20870.
6. E. Leviatan et al. (Qedma/RIKEN/IBM), "Resolving structure in prethermal Floquet dynamics with precision quantum computation," arXiv:2607.24937 (Jul 2026).
7. S. V. Barron et al. (IBM/Algorithmiq/Flatiron et al.), "Observable estimation in the absence of classical verification," arXiv:2607.25998 (Jul 2026).
8. J. Vovrosh et al., "Simulating dynamics of the two-dimensional transverse-field Ising model: a comparative study of large-scale classical numerics," arXiv:2511.19340 (Nov 2025).
9. Quantum Advantage Tracker, quantum-advantage-tracker.github.io; issues #149 and #227 on `su2_hadron_dynamics_lsh_x_100_meson` (R. Farrell, MPS $\chi=800$, 3,003 s, Apple M3 Pro).
10. F. Ilčić, R. Majumdar, E. Mathew, M. O. Ali, N. Earnest-Noble, I. Raychowdhury, "Observation of robust and coherent non-Abelian hadron dynamics on noisy quantum processors," arXiv:2602.18080 (Feb 2026).
11. L. Ebner, B. Müller, A. Schäfer, C. Seidl, X. Yao, "Eigenstate thermalization in 2+1 dimensional SU(2) lattice gauge theory," Phys. Rev. D 109, 014504 (2024), arXiv:2308.16202.
12. F. Turro, X. Yao, "Emergent hydrodynamic mode on SU(2) plaquette chains and quantum simulation," Phys. Rev. D 111, 094502 (2025), arXiv:2502.17551.
13. J.-W. Chen, Y.-T. Chen, G. Meher, B. Müller, A. Schäfer, X. Yao, "Local thermalization of SU(2) lattice gauge fields on quantum computers," arXiv:2603.23948 (Mar 2026).
14. B. Müller, X. Yao, "Simple Hamiltonian for quantum simulation of strongly coupled (2+1)D SU(2) lattice gauge theory on a honeycomb lattice," Phys. Rev. D 108, 094505 (2023), arXiv:2307.00045.
15. F. Turro, A. Ciavarella, X. Yao, "Classical and quantum computing of shear viscosity for (2+1)D SU(2) gauge theory," Phys. Rev. D 109, 114511 (2024), arXiv:2402.04221.
16. T. Laakkonen et al. (Quantinuum), "End-to-end quantum algorithms for the Jones polynomial," PRX Quantum 7, 020355 (2026), arXiv:2503.05625.
17. A. Ransford et al. (Quantinuum), "A 98-qubit trapped-ion quantum computer with all-to-all connectivity," Nature 655, 81 (2026).
18. E. Granet et al. (Quantinuum), "Superconducting pairing correlations on a trapped-ion quantum computer," arXiv:2511.02125 (2025).
19. S. A Rahman, R. Lewis, E. Mendicelli, S. Powell, "Self-mitigating Trotter circuits for SU(2) lattice gauge theory on a quantum computer," Phys. Rev. D 106, 074502 (2022), arXiv:2205.09247.
20. J. Tindall, M. Fishman, E. M. Stoudenmire, D. Sels, "Efficient tensor network simulation of IBM's Eagle kicked Ising experiment," PRX Quantum 5, 010308 (2024), arXiv:2306.14887.
21. M. S. Rudolph et al., "Pauli propagation: a computational framework for simulating quantum systems," arXiv:2505.21606 (2025).
22. T. Begušić, G. K.-L. Chan, "Real-time operator evolution in two and three dimensions via sparse Pauli dynamics," PRX Quantum 6, 020302 (2025), arXiv:2409.03097.
23. A. N. Ciavarella, C. W. Bauer, J. C. Halimeh, "Generic Hilbert space fragmentation in Kogut–Susskind lattice gauge theories," Phys. Rev. D 112, L091501 (2025), arXiv:2502.03533.

Context items also checked: T. Hayata, Y. Hidaka, Y. Kikuchi, q-deformed SU(2)$_3$ on Quantinuum H2-1 (4 plaquettes, ~650 two-qubit gates), Phys. Rev. Research 8, 033137 (2026), arXiv:2601.13530; G. C. Santra et al., quantum resources in non-Abelian LGTs, arXiv:2510.07385; DOE QCUP machine list (IBM, Quantinuum, IonQ, IQM).
