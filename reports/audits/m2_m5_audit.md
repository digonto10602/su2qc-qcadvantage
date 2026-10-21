# Combined gate audit — M2, M3, M4, M5 (independent reviewer)

Date: 2026-10-06. Scope: su2torus/{exact,trotter,hamiltonian,observables,thermal,noise,mps}.py,
scripts/{run_m3,run_m4,run_m5_check,resource_table,run_classical}.py, reports/m2_gate_counts.md,
reports/m3_physics.md, reports/m4_noise.md, results/m3/summary.json, results/m4/noise_runs.json,
results/m5/*.json.

Suite: `python3 -m pytest -q` -> 167 passed, 0 failed. `python3 tools/check_lock.py` -> `LOCK OK`.

## Verdicts

| milestone | verdict |
|---|---|
| M2 | PASS |
| M3 | BLOCKING (report text only; code and data are correct) |
| M4 | PASS |
| M5 | PASS |

## Independent checks run by the reviewer (scratch scripts, not committed)

1. **Real-time Chebyshev** (`exact.chebyshev_step` + fused `hamiltonian.cheb_recur`, bounds from
   `exact.spectral_bounds`), 3x3 torus (n = 18), g = 1.25, stripe, both models, 7 times in
   [0, 4g^2] plus a negative time: max state error vs `scipy expm_multiply` on the sparse H is
   1.1e-14.
2. **Bound padding.** Same lattice: the padded Lanczos bounds lie 1.05 (SU(2)) and 1.25 (abelian)
   energy units outside the true extremal eigenvalues (computed with eigsh tol 1e-12). Lanczos Ritz
   values always lie inside the spectrum and the tol = 1e-6 error is < 1e-4 here, so the 5% + 1e-3
   padding is safe by four orders of magnitude. The N = 24 bounds in `results/m3/parts/bounds_n24_*.json`
   are consistent with this (and energy is conserved to <= 6e-12 in all N = 24 runs, which would not
   happen if the spectrum left [-1, 1]).
3. **Imaginary-time Chebyshev** (`thermal.chebyshev_imag`): direction of e^{-beta H/2} v agrees with
   `expm_multiply` to 1.3e-14 at beta = 1.5 and 1.5e-15 at beta = -0.8. The coefficients
   ive(k, a)(-1)^k (2 - delta_k0) are the correct modified-Bessel expansion for both signs of a.
4. **Fast numba observables** (n > 16 path): Z, ZZ on all bonds, hexagon X-string and <H_B> agree
   with brute-force numpy to <= 4e-15 at n = 18.
5. **Symmetry claim (m3_physics.md point 3)** on the even-L1 4x2 torus (n = 16): abelian stripe and
   neel_half <H_E>(t) constant at 7.03125 (= <H_E>_inf) to 1e-14 over t in [0, 6]; SU(2) not
   constant (7.01..7.65). The N = 24 abelian data in summary.json have std 1.3e-13. Claim confirmed.
6. **MPS environment sweep** (`mps._observables`): at N = 12, chi = 256 and a deliberately truncated
   chi = 6 (error estimate 0.14), the sweep's <H_E> and <Z_p> equal those of the dense contraction of
   the same MPS to <= 9e-16 (so normalisation of a truncated MPS is handled). Vs the exact Trotter
   state at chi = 256: |dH_E| = 5e-12, |dZ| = 2e-9.
7. **Typicality statistics at N = 24** (SU(2) vacuum, the committed beta = 1.4401, committed bounds):
   8 fresh random vectors (seeds 5000-5007) give per-vector <H> from -0.83 to +1.05 and
   <H> = 0.19 +- 0.24, <H_E> = 5.27 +- 0.16 (jackknife). Implied std of a 4-vector estimate:
   ~0.35 in <H>, ~0.22 in <H_E> (~4%), i.e. ~0.03-0.04 in beta (d<H>/dbeta ~ -9).
   At N = 12 the same 4-vector estimator has std 0.40 in <H> (Boltzmann participation ratio only 5.7).

## BLOCKING findings

**B1 (M3) `reports/m3_physics.md:48` — claim "SU(2) thermalizes faster than its abelian twin" is not
supported as stated.** From summary.json (late = last 25%):
- SU(2) N = 12 deviations from the thermal anchor are 1.2-12.5%, not "2-9%" (vacuum g = 1.25: 2.178
  vs 2.488 = -12.5%).
- At g = 1.4 the abelian twin is *closer* to its anchor than SU(2): 0.3%, 0.9%, 0.7% (abelian) vs
  3.7%, 3.4%, 1.2% (SU(2)). At g = 1.25 abelian neel_half/stripe are 6.8%/5.8% off, comparable to SU(2).
  These abelian rows are omitted from the table ("in the JSON"), which hides the counter-evidence.
- Fix: add the abelian g = 1.1 and 1.4 rows to the table and rewrite point 1, e.g. "From the vacuum
  at g <= 1.25, SU(2) gets closer to its anchor than the abelian twin (12% vs 22% at g = 1.25, 9% vs
  44% at g = 1.1); at g = 1.4 both are within 4% and the abelian twin is closer." Adjust the "clear
  SU(2)/abelian difference" phrase at :52 accordingly.

**B2 (M3) `reports/m3_physics.md:52` — "E_E rises from 0 to about 40-45% of the decohered value" for
g ~ 1.25-1.4 is wrong at g = 1.4** (N = 12: 1.773/6.615 = 27%). Fix: "27% (g = 1.4) to 42% (g = 1.25)".

**B3 (M3) `reports/m3_physics.md:49` (and `:19`, `:38`) — N = 24 thermal anchors are quoted without
their statistical error, and the finite-size claim relies on a difference below it.** "1.448 against
1.440, so the anchor is already close to its large-N value" compares numbers whose 4-vector
uncertainty is ~+-0.03-0.04 in beta and ~+-0.2 in <H_E> (check 7). `reports/project_report.md:253,478`
quotes beta = 1.44008/1.440075 (6 digits). Fix: state the anchor as beta = 1.44 +- 0.04,
<H_E>_th = 5.15 +- 0.2 (or recompute with a jackknife over the 4 stored-per-vector values, see N2) and
soften point 2 to "consistent within the typicality error". The SU(2)-vs-thermal gap at N = 24
(4.47 vs 5.15, 13% +- 4%) survives, so the recommendation does not change.

## NON-BLOCKING findings

**N1 (M3) `su2torus/thermal.py:85-107`, `scripts/run_m3.py:127-134` — no error bar is produced.**
Fix: accumulate per-vector (den, num) in `thermal_typicality`, return a jackknife error, store it in
`thermal` as `electric_energy_err`/`beta_err`.

**N2 (M3) `su2torus/thermal.py:124,143` — 2% criterion with E0 = 0.** The tolerance is
0.25 * 0.02 * max(1, |E0|), i.e. an absolute 0.005 for the vacuum. All 24 runs meet 2%: N = 12 to
<= 4e-13 relative; N = 24 worst |dE| = 0.0021 (vacuum; 0.03% of E0 - E_gs ~ 6). The search is a
correct Illinois false position (monotone E(beta) for fixed vectors guarantees the bracket), uses
<= 6 of the 8 allowed evaluations, and the returned beta is the best evaluated point. But nothing
asserts the 2% match after `max_evals`. Fix: in run_m3.py raise/flag if |E_th - E0| > 0.02 *
max(1, |E0|), and document in METHODS/report that for E0 = 0 the scale is 1 energy unit.

**N3 (M3) `reports/m3_physics.md:50` — symmetry argument is missing one step.** U H U^dag = K - H and
U H_E U^dag = K - H_E (K = 2<H_E>_inf) give <H_E>_psi(t) = K - <H_E>_{U psi}(-t). Constancy also needs
time reversal (H, H_E and the product state are real, so <H_E>(-t) = <H_E>(t)) and translation
invariance (U psi is the stripe/half shifted by one/two rows). Fix: add "and, because H and the
initial state are real, <H_E>(t) is even in t". Conclusion is correct (check 5). Also say why it
fails for SU(2): Pi_A X changes the neighbour-dependent amplitudes c(n_q) of the B flips.

**N4 (M3) METHODS §5 deviation.** Typicality uses a Chebyshev e^{-beta H/2} instead of `expm_multiply`.
Verified identical to 1e-14 (check 3); documented in STATUS and the report. Keep it listed for
Digonto's acknowledgement.

**N5 (M4) `su2torus/noise.py:88-91`, `scripts/run_m4.py:192` — F_pred and F_echo do not see the
same noise.** F_pred is computed on `evo` (no X preparation gates, no readout); F_echo's return
probability includes readout ((1-r)^12 = 0.988, lowering F_echo by ~0.6%) and the X gates. Harmless
against the 25% criterion (observed 2.8-13.4%), but state it in reports/m4_noise.md, or include
readout in F_pred as an option.

**N6 (M4) `su2torus/noise.py:40` — `_readout_rate` reads the private attribute
`NoiseModel._default_readout_error`.** May break silently (returning 0) on an Aer upgrade. Fix: pass
the readout rate explicitly to `run_noisy`/`echo_fidelity`, or assert the attribute exists when the
model has readout errors.

**N7 (M4) conventions verified.** Aer `depolarizing_error(p, k)` is (1-p) rho + p I/d; its probability
of a non-identity Pauli is p(1-1/d^2) = 15p/16 (2q), 3p/4 (1q): correct. Folding G -> G (G^-1 G)^k
scales each gate's error by 1, 3, 5; Aer `run` does not transpile, so folded pairs are not cancelled
(confirmed by the ZNE values: deviations 0.036, 0.097, 0.154 at factors 1, 3, 5). Richardson weights
15/8, -5/4, 3/8 are correct. ZNE error / raw error = 0.065-0.172 (4/4 configs). Readout error is not
folded, so it remains in the ZNE residual (bias of order 2r per bond); mention in the report.

**N8 (M5) `scripts/resource_table.py:37` — raw fidelity e^{-p2 N_2q}** follows PLAN but differs from
the M4 per-gate infidelity 15p2/16; add a one-line note to resources.json/collaborator summary so the
two numbers are not compared directly.

**N9 (M5) `su2torus/mps.py:78,109-115` — performance only.** `_site_arrays` copies every tensor to
host memory after every step (GPU backend), and each bond's ZZ re-contracts the transfer chain from
scratch (O(n_bonds x span x chi^3)). Fine at N <= 48, chi <= 256; for N = 96 on GPU, keep the sweep on
the device and reuse left environments.

**N10 (M2)** `reports/m2_gate_counts.md:24`: "Gray-code minimum of 2^k cx" is the minimum for a
*generic* uniformly controlled rotation; add "generic" (the SU(2) angles have product structure,
which has not been exploited). Counts 9.5n/11n (order 1), 15n/18n and 13.5n/15n (order 2) verified by
hand: 1.5n bonds, 8 cx per flip with the closing toggle.

## What was found correct (no action)
- Chebyshev series `exact.cheb_series` (recurrence, truncation kmax = x + 10x^{1/3} + 40 with
  1e-15 early stop), sign handling for t < 0, global phase e^{-ict}.
- numba kernels: flip amplitude depends only on neighbours of p (unchanged by the flip), so H is
  symmetric; per-chunk reductions have no write races.
- Abelian/SU(2) beta = 0 for stripe and neel_half on even-L1 tori: exactly half (18/36) of the bonds
  are cut and <H_B>(0) = 0, so E0 = Tr H / D. Also true for the 4x6 stripe (36/72).
- M4 echo F within 25% (max 13.4%); seeded 4000-shot repro block present.
- M5: snake order along the shorter side; gate matrices via `Operator(op).reverse_qargs()` match
  quimb's big-endian convention (N = 24 fidelity 1 - 1.1e-9; error estimate 3.4e-9 is conservative).
