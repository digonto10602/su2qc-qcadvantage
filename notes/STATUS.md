# STATUS (coordinator-only; keep < 60 lines)

## Position
- Milestone: M6 gate audit in progress; then tag v0.1-pre-cluster and STOP
- Locked acceptance milestones: m1–m6 (full suite: 167 passed)
- Tags (local only, see Notes): m1-done … m5-done

## Attempt log (failing acceptance test -> attempts)
- (none; no test needed escalation)

## Milestone results
- M1: lattice + H (sparse, Pauli, numba matrix-free). Audit PASS.
- M2: order-1 step = 9.5n (native ZZ) / 11n (CZ), exactly as claimed; order 2: 15n/18n (13.5n/15n fused). Audit PASS.
- M3: 18 N=12 + 6 N=24 runs + open string; E conserved < 1e-7. Audit: report wording fixed (B1–B3). Key finding: stripe/neel_half have β=0 on even-L1 tori (incl. N=48 frontier torus); abelian <H_E>(t) exactly constant there (symmetry).
- M4: N=12 exact density matrix; ZNE cuts <H_E> error 5.8–15.5x (4/4); echo F within 3–13% of prod(1-15p2/16); global rescaling makes errors worse. Audit PASS.
- M5: MPS (quimb) N=24 chi=256 1-F=1.1e-9; resource table N=48–96; config entry point + Slurm placeholders. Audit PASS.
- M6: N=48 stripe chi 32–256 (chi=256 run 6935 s ≈ 1.9 h, under the 2 h cap). All chi pairs converge to step 2, fail from step 3. Summary: 20 verified numbers.

## Notes
- Git tags do not reach the remote (proxy accepts branch pushes only). Tags exist locally.
- Thread oversubscription (numba + OpenBLAS) slowed SVDs ~400x; MPS runs use OPENBLAS_NUM_THREADS=1.
- N=24 thermal anchors: Chebyshev e^{-βH/2} instead of expm_multiply (METHODS §5), equal to 1e-9; 4-vector typicality error ≈ ±0.2 in <H_E>, ±0.04 in β.
- M5 check uses dt=0.2 (reviewer choice; at dt≈0.7 chi=256 cannot reach 1-1e-6 at N=24).
- Non-blocking audit items open: jackknife errors in typicality; noise.py reads a private Aer attribute; F_echo includes readout but F_pred does not; resource table uses e^{-p2 N} vs M4's 15p2/16; MPS observable sweep needs optimisation before N=96.

## Decisions waiting for Digonto
1. Initial state for the physics comparison at N=48–64: PLAN prescribes stripe, but it is at β=0 (thermal = decohered) on even-L1 tori. Recommend the SU(2) vacuum (β≈1.4) or an odd-row state. Evidence: reports/m3_physics.md §3.
2. Open-lattice H_E form (audit N1): code uses the bond form; METHODS' 3n_p form differs on open lattices (string E_E(0)=4.41 vs 5.88). Evidence: reports/audits/m1_audit.md.
3. Tags must be pushed from a machine with tag-push rights (git push origin --tags).
