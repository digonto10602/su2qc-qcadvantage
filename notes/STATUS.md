# STATUS (coordinator-only; keep < 60 lines)

## Position
- Milestone: ALL DONE (M1–M6). Tagged v0.1-pre-cluster. STOPPED per CLAUDE.md (no cluster work).
- Locked acceptance milestones: m1–m6 (full suite: 171 passed after post-M6 fixes; merged to main in PR #1)
- Tags on GitHub: m1-done … m6-done, v0.1-pre-cluster (pushed from laptop via scripts/push_tags.sh); v0.1.1 = main after post-M6 fixes (b856dd3); all 8 tags on GitHub

## Attempt log (failing acceptance test -> attempts)
- (none; no test needed escalation)

## Milestone results
- M1: lattice + H (sparse, Pauli, numba matrix-free). Audit PASS.
- M2: order-1 step = 9.5n (native ZZ) / 11n (CZ), exactly as claimed; order 2: 15n/18n (13.5n/15n fused). Audit PASS.
- M3: 18 N=12 + 6 N=24 runs + open string; E conserved < 1e-7. Audit: report wording fixed (B1–B3). Key finding: stripe/neel_half have β=0 on even-L1 tori (incl. N=48 frontier torus); abelian <H_E>(t) exactly constant there (symmetry).
- M4: N=12 exact density matrix; ZNE cuts <H_E> error 5.8–15.5x (4/4); echo F within 3–13% of prod(1-15p2/16); global rescaling makes errors worse. Audit PASS.
- M5: MPS (quimb) N=24 chi=256 1-F=1.1e-9; resource table N=48–96; config entry point + Slurm placeholders. Audit PASS.
- M6: N=48 stripe chi 32–256 (chi=256: 6935 s, under the 2 h cap). Converged through step 2, marginal at 3, unconverged from 4; chi=256 cumulative error 0.79 at step 8. Audit PASS after text fixes. Summary: 21 verified numbers.

## Notes
- Cloud git proxy accepts branch pushes only: tags are pushed from a laptop with scripts/push_tags.sh.
- Thread oversubscription (numba + OpenBLAS) slowed SVDs ~400x; MPS runs use OPENBLAS_NUM_THREADS=1.
- N=24 thermal anchors: Chebyshev e^{-βH/2} instead of expm_multiply (METHODS §5), equal to 1e-9; 4-vector typicality error ≈ ±0.2 in <H_E>, ±0.04 in β.
- M5 check uses dt=0.2 (reviewer choice; at dt≈0.7 chi=256 cannot reach 1-1e-6 at N=24).
- Post-M6 cleanups (all audit non-blocking items + code review, 10 findings) done: jackknife typicality errors (N=24 SU(2) vacuum anchor 5.15±0.27); public-API readout rate; readout-corrected echo (excess vs F_pred is not readout); resource table has both fidelity conventions; MPS observables on-device, resumable run_classical, batch set -e. See reports/audits/ and git log.

## Decisions waiting for Digonto
1. Initial state for the physics comparison at N=48–64: PLAN prescribes stripe, but it is at β=0 (thermal = decohered) on even-L1 tori. Recommend the SU(2) vacuum (β≈1.4) or an odd-row state. Evidence: reports/m3_physics.md §3.
2. Open-lattice H_E form (audit N1): code uses the bond form; METHODS' 3n_p form differs on open lattices (string E_E(0)=4.41 vs 5.88). Evidence: reports/audits/m1_audit.md.
3. Trotter step 2δt/g²=0.9 (PLAN M6) departs strongly from e^{-iHt} (audit: N=12/16 Trotter <H_E> 0.54 vs exact 0.31 of E_inf at step 8). Decide δt before cluster QMC thermal-anchor work. Evidence: reports/audits/m6_audit.md.
4. ~~Tags~~ DONE: all 8 tags (m1-done … m6-done, v0.1-pre-cluster, v0.1.1) pushed from the laptop with scripts/push_tags.sh.
