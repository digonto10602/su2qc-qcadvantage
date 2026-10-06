# M6 gate audit (final cloud milestone, before tag `v0.1-pre-cluster`)

Reviewer: independent `reviewer` subagent. I did not modify any implementation file.

## Verdict: PASS-WITH-FIXES

The data, the scripts and the tests are sound. Three statements in `reports/collaborator_summary.md` are not supported by the data, or leave out a caveat that changes the team decision. All three need text edits only, with no reruns. Fix them before the tag.

## Checks run
| check | result |
|---|---|
| `python3 -m pytest -q` | 167 passed (exit 0) |
| `python3 tools/check_lock.py` | LOCK OK |
| `scripts/check_summary.py` | 20 markers, 0 failures |
| `scripts/analyze_convergence.py --out <scratch>` | byte-identical to `results/m6/convergence.json` |
| $\delta t$ | $2\cdot0.703125/1.25^2=0.900$ ✓ (PLAN M6) |
| decohered $H_E$ | $\tfrac{3g^2}{16}\cdot\tfrac{3N}{2}=21.094$; stripe step 0 = 21.094 ✓ ($\beta=0$) |
| gate count | 456 = 9.5·48 ✓ |
| ZNE improvement (M4 JSON) | 10.14, 15.48, 6.78, 5.80 → "5.8–15.5×" ✓ |
| per-run cap | the χ=256 run took 6935 s < 7200 s ✓ (only just) |
| prior audit B1–B3 (`reports/m3_physics.md`) | B1 fixed (:54–62, full table and coupling-dependent wording); B2 fixed (:66, 77/41/27%, N=24 42%); B3 fixed (:51, :63, ±0.04 in β, ±0.2 in $H_E$). Residual in `project_report.md`, see N10. |

### Convergence data (χ=128→256 pair; recomputed from `n48_stripe_g1.25.json`)
| step | t | max\|ΔZ\| | mean\|ΔZ\| | #plaq >0.01 | ΔE | χ=256 1−F (cumulative) | retained fraction this step F_k/F_{k−1} |
|---|---|---|---|---|---|---|---|
| 2 | 1.41 | 0.0007 | 0.0002 | 0 | 0.007 | 0.0002 | 0.9998 |
| 3 | 2.11 | 0.013 | 0.004 | 4/48 | 0.032 | 0.014 | 0.986 |
| 4 | 2.81 | 0.037 | 0.015 | 32/48 | 0.030 | 0.124 | 0.889 |
| 5 | 3.52 | 0.049 | 0.017 | 35/48 | 0.049 | 0.300 | 0.799 |
| 6 | 4.22 | 0.031 | 0.013 | 31/48 | 0.033 | 0.501 | 0.713 |
| 8 | 5.63 | 0.067 | 0.015 | 33/48 | 0.030 | 0.795 | 0.640 |

How max|ΔZ| at step 3 changes as χ doubles: 0.037 → 0.023 → 0.013, roughly halving each time. A χ=512 run would probably converge at step 3. At step 4 it goes 0.057 → 0.050 → 0.037, so convergence is slow. The ⟨H_E⟩ criterion is met at **all 8 steps for all pairs** (`last_converged_step_electric_energy = [8,8,8]`).

**Is "MPS fails after 2–3 steps at χ ≤ 256" justified?** Yes, with a qualification. The claim rests on the chi-to-chi comparison of local $Z_p$, not on the cumulative error estimate, so the cumulative nature of the estimate does not bias it. Requiring agreement at all steps up to K (the prefix rule) is the right choice, because truncation errors propagate forward. The per-step retained fraction gives an independent confirmation: it drops below 0.99 at step 3 and to about 0.64 per step by steps 7–8, so the onset is real. The precise statement is: converged through step 2 ($t=1.41$), marginal at step 3 (4 of 48 plaquettes exceed the tolerance, by at most 0.003), and clearly unconverged from step 4.

## BLOCKING findings (text-only fixes)

**B1 `reports/collaborator_summary.md:67`: the explanation of why the energy is insensitive is wrong.** The summary says "⟨H_E⟩ stays near the decohered value (≈21.09), because this state starts at β=0". In the data, ⟨H_E⟩ swings from 21.0 to 23.63, up to 12% away from 21.09, while the χ=128/256 difference never exceeds 0.049. Truncation does not pull ⟨H_E⟩ toward the β=0 value. The real reason is that ⟨H_E⟩ is a sum over 72 bonds, so its tolerance of 0.211 amounts to a *mean* per-bond tolerance of about 0.01, whereas the $Z_p$ test is a *maximum* over 48 sites. Mean |ΔZ| for 128→256 also stays ≤ 0.017. Fix: "⟨H_E⟩ agrees between successive χ to within the 1% tolerance at all 8 steps (max difference 0.049), even though it swings by up to 2.5 around 21.09. It is a sum over 72 bonds, so local errors average out. The maximum over the 48 local magnetizations is the sensitive probe." Also state the fact `last_converged_step_electric_energy = [8,8,8]` with a marker.

**B2 `reports/collaborator_summary.md:66`: the cumulative error estimate is over-read.** "error estimate 0.79, so most of the state's weight has been truncated away" presents the estimate as a direct measure of accuracy. It is actually $1-\prod_i(1-w_i)$ over every truncation since $t=0$ (`su2torus/mps.py:146`, quimb `fidelity_estimate`), i.e. a heuristic *global* infidelity. For 48 qubits, global fidelity falls roughly exponentially with N even when local observables are still accurate. Fix: "By step 8 the estimated overlap of the χ=256 state with the untruncated Trotter state is 0.21. This estimate is cumulative; per step, the fraction of weight kept falls from 0.986 (step 3) to about 0.64 (steps 7–8). Local observables degrade more slowly (mean |ΔZ| ≈ 0.015), but they are not converged." Also replace "fails after 2–3 Trotter steps" (:69) with "is converged through step 2, marginal at step 3, and unconverged from step 4. The step-3 differences halve with each doubling of χ, so χ≈512 may recover step 3."

**B3 `reports/collaborator_summary.md:79,82` (§7 item 4 and the "Decision"): the Trotter step size is not taken into account.** At $2\delta t/g^2=0.9$ the first-order circuit is a coarse periodic (Floquet) drive. It does not conserve $H$, so a thermal anchor of $H$ at $\beta\approx1.44$ (the QMC item) is not the late-time reference *for this circuit*. Independent check in this audit (`split_sparse` + `expm_multiply`, $g=1.25$, $\delta t=0.703$, vacuum, steps 1–10):
- N=12: Trotter ⟨H_E⟩/E_∞ against exact: step 2 0.92 vs 0.75, step 3 0.56 vs 0.38, step 8 0.54 vs 0.31. The drift of ⟨H⟩ is +0.7 to +3.1 (E_∞ = 5.27).
- N=16: the same pattern (step 8: 0.54 vs 0.32).
- The stripe ($\beta=0$) is barely affected.

The recommended vacuum target is therefore strongly affected by Trotter error at this δt. Fix: add to §8 and to roadmap item 4: "At $2\delta t/g^2=0.9$ the circuit deviates strongly from continuous-time dynamics from the vacuum (N=12: up to 23% of the decohered ⟨H_E⟩). A thermal anchor of $H$ applies only at smaller δt or second order. For the coarse circuit, the classical reference must simulate the same circuit (MPS/PEPS/Pauli propagation). Decide the δt for the physics comparison before QMC time is spent."

## NON-BLOCKING findings

- **N1 `collaborator_summary.md:68`**: "6935 s on 4 cloud vCPUs". The runs used `OPENBLAS_NUM_THREADS=1` and overlapped the M3 numba batch (the M6 log runs 04:22–06:43 and `logs/m3.log` ends at 05:28). Timings are therefore single-threaded with contention, not a clean benchmark. Say so, because it affects the χ=512+ extrapolation (the observed ratios 3.5/4.5/5.8× per doubling are approaching χ³ = 8×).
- **N2 `collaborator_summary.md:55`**: "With mitigation, 4–6 steps at N = 48–64 are realistic". The raw F at N=64, 6 steps, $p_2=10^{-3}$ is $e^{-3.65}=0.026$, and at N=48, 4 steps it is 0.161. ZNE was validated only down to F≈0.33 at N=12 (M4). State the validated range, or write "plausible, untested below F≈0.33".
- **N3 `collaborator_summary.md:28`**: "relaxes toward a thermal state" omits the fact (from `m3_physics.md:62`) that at N=24 the late-time mean is still 13%±4% below the anchor and fig. 1 shows undamped oscillations to $t=4g^2$. Add one clause.
- **N4 unmarked data numbers**: 1140 (:38, `noise_runs.json configs/steps=10/n_2q`), ≈21.09 (:67, `convergence.json steps/0/electric_energy/0`), 0.703 (:62, `config/dt`), and the derived values 5.8–15.5× (:48), ±0.04 (:28, from an audit re-run not stored in any JSON) and 2.3× (:86). Add markers where possible. Store the typicality spread in `results/m3/summary.json` or cite its source.
- **N5 fig 4 (`scripts/make_figures.py:110-118`)**: the tolerance lines (:114–115) have no legend entries. Add `label="tol E (1%)"` and `label="tol Z (0.01)"`. Step-0 zeros are floored to 1e-16, which squeezes the informative range 1e-4…1 into the top quarter of the plot: use `set_ylim(1e-5, 1)` or start at step 1. Change the right panel's ylabel to "cumulative 1 − fidelity estimate".
- **N6 fig 3**: the N=96/p₂=10⁻³ and N=64/p₂=1.5×10⁻³ curves coincide exactly (same $p_2N$), which hides the green dotted line. Mention this in the caption.
- **N7 fig 5**: the bold titles "PEPS with belief-propagation gauging" and "Pauli propagation of local observables" overlap the "N = …" column. Move the column right (x ≈ 0.31) or wrap the titles.
- **N8 README reproducibility**:
  - (a) `logs/` is git-ignored, so on a fresh clone `nohup scripts/run_m3_batch.sh > logs/m3.log` (README:27) fails. Add `mkdir -p logs`.
  - (b) `numba` is missing from `requirements.txt`, although the M3 N=24 timing (~4 h) depends on it (the fallback in `hamiltonian.py:19-22` is silent).
  - (c) No versions are pinned. The tested set is qiskit 2.5.2, qiskit-aer 0.17.2, quimb 1.15.0, numba 0.68.0. Pin them or list them.
  - (d) README:31 should note that the ~2.3 h is the total over four χ runs and that χ=256 alone is ~1.9 h.
- **N9 fig 1**: correct and readable. At N=24 the stripe and Néel-half thermal lines are hidden under the decohered line ($\beta=0$). Say so in the caption.
- **N10 `reports/project_report.md`** (linked from README:5 as the "Full report") is stale:
  - :249 says "runs finished so far … pending";
  - :688 says the typicality uncertainty is "not quantified";
  - :700 says "Running now";
  - :253 and :478 quote β = 1.44008 / 1.440075 with no error bar (residual of the previous B3).

  Either update it or mark it as a historical build log in the README.

## Data-to-claim check (summary §3–§7)
All 20 marked numbers match the JSON. Fig. 2 and the ZNE ranges match. The claim "global-fidelity rescaling over-corrects" is supported (echo-rescaled error > raw error in all four configurations). The stoquasticity claim is correct, since off-diagonal elements are $-\tfrac1{g^2}\prod c\le0$. Figures (i)–(v) are referenced, and their data files exist.
