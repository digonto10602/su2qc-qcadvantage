# METHODS — conventions and specification (binding for code and tests)

Source physics: Ciavarella, de Putter, Younis, Rrapaj, arXiv:2608.28752 (minimal local-Krylov truncation of pure SU(2) on a triangular lattice, one qubit per plaquette); background and references in `docs/idea.md`.

## 1. Lattice

Triangular-lattice sites are $(i,j)$ with $0\le i<L_1$ and $0\le j<L_2$. Each cell owns two plaquettes:

- up triangle $\{(i,j),(i+1,j),(i,j+1)\}$: kind 0, sublattice A
- down triangle $\{(i+1,j),(i,j+1),(i+1,j+1)\}$: kind 1, sublattice B

The plaquette index is $p = 2(iL_2+j)+\text{kind}$, so $n=2L_1L_2$.

Neighbours: $\mathrm{up}(i,j)$ shares a link with $\mathrm{down}(i,j)$, $\mathrm{down}(i-1,j)$ and $\mathrm{down}(i,j-1)$. Indices wrap mod $L_1, L_2$ on the torus; on open lattices, out-of-range neighbours are dropped. The plaquette graph is the honeycomb graph and is bipartite (A–B).

The hexagon around site $(i,j)$ is $\{\mathrm{up}(i,j),\mathrm{up}(i-1,j),\mathrm{up}(i,j-1),\mathrm{down}(i-1,j),\mathrm{down}(i,j-1),\mathrm{down}(i-1,j-1)\}$, listed in cyclic order.

## 2. Hamiltonian and basis

Qubit $p$ represents plaquette $p$. $|1\rangle$ means excited ($n_p=1$, $Z_p=-1$). The basis index is $b=\sum_p n_p 2^p$, the Qiskit little-endian convention.

$$H_E=\frac{3g^2}{8}\sum_p\Big(3n_p-\sum_{q\in nb(p)}n_pn_q\Big)=\frac{3g^2}{16}\sum_{\langle pq\rangle}(1-Z_pZ_q)=\frac{g^2}{2}\sum_{\langle pq\rangle}E^2_{pq},\qquad E^2_{pq}=\tfrac38(1-Z_pZ_q)$$

$$H_B=\sum_p a_p X_p,\qquad a_p=-\frac{1}{g^2}\prod_{q\in nb(p)}c(n_q),\quad c(0)=1,\ c(1)=\tfrac12 .$$

Abelian twin: $c\equiv1$, which is the transverse-field Ising model on the honeycomb graph. $H_B^A$ and $H_B^B$ are the sublattice parts of $H_B$. Both parts are sums of mutually commuting terms, because each term's controls lie on the other sublattice.

Pauli form, using $c(n_q)=(3+Z_q)/4$: $\prod_{q}(3+Z_q)/4=4^{-k}\sum_{S\subseteq nb(p)}3^{k-|S|}\prod_{q\in S}Z_q$, where $k=|nb(p)|$.

## 3. Trotter circuits and gate counting

Writing $U_X(t)=e^{-iH_Xt}$:

- order 1: $U(\delta t)=U_{B}(\delta t)\,U_{A}(\delta t)\,U_E(\delta t)$, with $U_E$ applied first
- order 2: $U(\delta t)=U_E(\tfrac{\delta t}{2})U_A(\tfrac{\delta t}{2})U_B(\delta t)U_A(\tfrac{\delta t}{2})U_E(\tfrac{\delta t}{2})$

Equality is required up to a global phase.

- $U_E$: one `rzz` per bond, angle $-2\cdot\tfrac{3g^2}{16}\,t$ (the constant term only adds a global phase).
- Flip of $p$: $e^{-i t a_p X_p}=H_p\,e^{-it a_p Z_p}\,H_p$, with $a_p=\sum_S w_S Z_S$ and $w_S=-\frac{1}{g^2}3^{k-|S|}/4^k$. Implement it as a Gray-code uniformly controlled $R_z$: `rz(2 t w_S)` on $p$ while the parity register on $p$ holds $S$, with a `cx`(neighbour → p) for every toggle. That is $2^k$ `cx` in total (8 on the torus).
- Abelian flip: `rx(2 t (-1/g^2))`.

Gate-count convention:

- native ZZ (trapped ions): `cx` = 1 entangler and `rzz` = 1
- CZ hardware: `cx` = 1 and `rzz` = 2

The expected first-order torus counts are $9.5\,n$ (native ZZ) and $11\,n$ (CZ). Acceptance requires $\le$ these values.

## 4. Observables and initial states

Observables: $\langle Z_p\rangle$, link energies $\langle E^2_{pq}\rangle$, connected $\langle Z_pZ_q\rangle_c$, $\langle H_E\rangle$, $\langle H_B\rangle$, and the hexagon string $\langle\prod_{p\in\text{hex}}X_p\rangle$.

Infinite-temperature (fully decohered) values: link energy $3/8$, $\langle H_E\rangle=\frac{g^2}{2}\cdot\frac38\cdot n_\text{bonds}$, and $0$ for the others.

Initial product states (all gauge-invariant):

- `vacuum`: nothing excited
- `neel_half`: $\mathrm{up}(i,j)$ for $i<\lfloor L_1/2\rfloor$
- `stripe`: $\mathrm{up}(i,j)$ for even $i$
- `string`: $\{\mathrm{index}(k,i,\lfloor L_2/2\rfloor)\}$ over all $i$ and both kinds (open lattices)

Tori used: $N=12\to(L_1,L_2)=(3,2)$ and $N=24\to(4,3)$. Open string lattice: $3\times4$.

## 5. Thermal anchors

The late-time value of a thermalizing run is the canonical average at the $\beta$ for which $\langle H\rangle_\beta$ equals the initial energy ($\beta$ may be negative). For $n\le14$ use full diagonalization. For $n=24$ use quantum typicality with matrix-free $e^{-\beta H/2}$ (`expm_multiply` on a `LinearOperator`; pass both `matvec` and `rmatvec`, since $H$ is Hermitian). Use $\le 8$ bisection steps and 4–8 random vectors. The required energy match is 2%.

## 6. M3 output schema — `results/m3/summary.json`

```json
{"runs": [{"N": 12, "g": 1.25, "state": "stripe", "model": "su2" | "abelian",
           "times": [...], "electric_energy": [...], "magnetic_energy": [...],
           "link_energies": [[...per bond...] per time], "z": [[...] per time] (optional for N=24),
           "thermal": {"beta": b, "electric_energy": x, "magnetic_energy": y},
           "infinite_temperature": {...section 4...},
           "flags": ["observable names whose late-time value is within 10% of the decohered value"]}]}
```

Required grid:

- $N=12$: $g\in\{1.1,1.25,1.4\}$ × 3 states × 2 models
- $N=24$: $g=1.25$ × 3 states × 2 models

Time grid: $t\in[0,4g^2]$ with $\ge21$ points. Write $\ge6$ figures to `results/m3/figures/`.

Open string: `results/m3/string_open_3x4_g1.4.json` with `times`, `electric_energy`, `magnetic_energy`, `z`, `link_energies`.

## 7. Performance hints (cloud machine: 4 vCPU, 16 GB)

- **Matrix-free $H$ at $n=24$.** Reshape $\psi$ to $(2,)^n$, where axis $n-1-p$ is qubit $p$. A flip is `np.flip` along that axis times a broadcast amplitude over the neighbour axes. Precompute the diagonal once (134 MB). Never build $2^{24}\times2^{24}$ matrices.
- **Batch and background long runs.** Put all long runs in one script started with `nohup ... &` and a log file. Poll with `sleep 590; tail -5 log`, not with frequent short checks.
- **Memory.** A complex vector at $n=24$ is 268 MB. Keep at most about 20 of them alive.
