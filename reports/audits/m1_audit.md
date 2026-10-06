# M1 gate audit (reviewer, independent)

Scope: `su2torus/lattice.py`, `su2torus/hamiltonian.py`, `tests/unit/test_hamiltonian_unit.py`, against METHODS sections 1-2 and PLAN M1.

**Verdict: PASS.** No blocking findings. There are four non-blocking items. Item N1 is a physics and spec question for Digonto.

## What was checked
- `python3 tools/check_lock.py` printed `LOCK OK`.
- `pytest tests/acceptance/m1 tests/unit`: 51 passed.
- Full suite: 63 passed, 9 failed, 24 errors. Every failure and error is in `tests/acceptance/m3`, which is not implemented yet. M1 and M2 are green.
- **Bit order.** The code uses $b=\sum_p n_p2^p$, with $|1\rangle$ meaning excited and $Z=-1$. `build_pauli` uses `from_sparse_list`, which is little-endian, and matches the brute-force reference to about $4\times10^{-15}$ on open (3,3) and (2,3) lattices, in addition to the locked torus checks.
- **Neighbour rule.** It follows METHODS exactly: up(i,j) is linked to down(i,j), down(i-1,j) and down(i,j-1), with wrapping on the torus and dropping on open lattices. The graph is bipartite by construction.
- **Hexagon cyclic order.** I checked it by hand. In the order up(i,j), down(i,j-1), up(i,j-1), down(i-1,j-1), up(i-1,j), down(i-1,j), each consecutive pair shares a link, and so do the last and first entries.
- **$H_E$ diagonal.** $(3g^2/16)\sum(1-ZZ) = (3g^2/8)\times$(number of cut bonds). This is correct. On a torus it is identical to the $3n_p-\sum n_pn_q$ form (checked numerically on (3,2): difference 0).
- **Flip amplitude.** $-(1/g^2)\,0.5^m$, where $m$ is the number of excited neighbours of $p$. In the abelian case it is $-(1/g^2)$. Correct.
- **Pauli coefficients.** $-(1/g^2)\,3^{k-|S|}/4^k$ on $X_p\prod_{q\in S}Z_q$. Correct. $X_p$ commutes with $Z_S$ because $p\notin nb(p)$, so the operator order does not matter.
- **Gather vs scatter in `apply_h`.** The gather form is valid. $H_{b,\,b\oplus 2^p}$ depends only on the occupations of the neighbours of $p$, and flipping $p$ does not change them. So $H$ is real symmetric, and gather equals scatter. Independent check at $n=24$ on (4,3):
  - numba vs numpy fallback: maximum difference 0.0, for both SU(2) and abelian;
  - $|\langle\phi|H\psi\rangle-\langle H\phi|\psi\rangle| < 6\times10^{-17}$;
  - time per warm call: 0.55 s with numba (the first call takes 6.5 s because of JIT compilation and cache loading) and 6-10 s with the numpy fallback.

## BLOCKING
None.

## NON-BLOCKING

**N1. The electric energy on open lattices is a spec ambiguity with physics consequences.**
- Where: `su2torus/hamiltonian.py:33-39` and METHODS section 2, line 22.
- The problem: METHODS writes $\frac{3g^2}{8}\sum_p(3n_p-\sum_q n_pn_q)$ "=" $\frac{3g^2}{16}\sum_{bonds}(1-Z_pZ_q)$. These two expressions are equal only when every plaquette has degree 3.
- Size of the effect: on open lattices, the $3n_p$ form adds $\frac{3g^2}{8}\sum_p(3-\deg p)\,n_p$. This term is the energy of the boundary links, assuming vacuum outside, which is the same assumption the code makes in $H_B$ by setting $c=1$ for dropped neighbours. For the M3 open `string` state on the 3x4 lattice at $g=1.4$:
  - bond form: $E_E(0)=4.41$;
  - $3n_p$ form: $E_E(0)=5.88$.
- What is locked: the code uses the bond form. That matches the locked `_ref.py` and `tests/acceptance/m3/test_m3_outputs.py:24-27,76`, so it must not change unilaterally.
- Fix:
  - Record this in STATUS under "Decisions waiting for Digonto". The comparison with Ciavarella's open-lattice string-breaking curves in M3/M6 may need the $3n_p$ form, which would be an optional `boundary="vacuum"` flag.
  - Change the METHODS wording to "= (torus)".

**N2. Degenerate tori are accepted silently.**
- Where: `su2torus/lattice.py:31-41`.
- The problem: with $L_1=1$ or $L_2=1$ on a torus, the neighbour sets merge multiple shared links into one, so the degree drops below 3 and the Hamiltonian is wrong.
- Fix: add `if self.periodic and min(L1, L2) < 2: raise ValueError`.

**N3. Memory guard for the sparse builders.**
- Where: `su2torus/hamiltonian.py:51-70`.
- The problem: `split_sparse` and `build_sparse` build int64 COO arrays of size about $n\cdot 2^n$. At $n=24$ that is more than 10 GB, close to the 12 GB stop rule.
- Fix: add a guard such as `if lat.n > 20: raise ValueError("use apply_h")`, and state in the docstring that `apply_h` is the $n=24$ path.

**N4. Hexagons on open lattices are not specified.**
- Where: `su2torus/lattice.py:56-66`.
- The problem: the code returns only complete hexagons on open lattices. That is sensible, but METHODS does not specify it.
- Fix: add one sentence to METHODS section 1, or keep the docstring as the convention.
