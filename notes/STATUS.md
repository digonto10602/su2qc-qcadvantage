# STATUS (coordinator-only; keep < 60 lines)

## Position
- Milestone: M3 (N=24 batch running: logs/m3.log); M4 tests being written by reviewer
- Locked acceptance milestones: m1, m2, m3
- Last tag: m1-done (local only, see note)

## Attempt log (failing acceptance test -> attempts)
- (none)

## Milestone results
- M1: lattice + H (sparse, Pauli, numba matrix-free 0.8 s/matvec at n=24). Audit PASS (reports/audits/m1_audit.md).
- M2: Chebyshev/expm_multiply evolution; Trotter circuits. Gate counts: order 1 = 9.5n (ZZ) / 11n (CZ) exactly (claim verified); order 2 = 15n/18n per isolated step, 13.5n/15n fused. reports/m2_gate_counts.md

## Notes
- Git tags do not reach the remote (proxy accepts branch pushes only; `git push origin <tag>` silently no-ops). Tags exist locally.
- Audit N1 (for Digonto, non-blocking): METHODS' 3n_p form of H_E equals the bond form only with 3 neighbours per plaquette; on open lattices the code (and locked tests) use the bond form (no boundary-link energy). Open-string E_E(0)=4.41 (bond) vs 5.88 (3n_p) at g=1.4.
- METHODS §5 asks for expm_multiply in typicality; at n=24 it took ~186 s per vector per beta. N=24 uses an imaginary-time Chebyshev expansion instead (same result to 1e-9, unit test); root search = bracketed false position (Illinois), <= 8 evaluations.

## Decisions waiting for Digonto
- (none)
