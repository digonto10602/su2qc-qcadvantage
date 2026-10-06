# M3 — small-torus physics (exact dynamics, N = 12 and N = 24)

## What was asked
Compute the exact real-time dynamics of the minimal SU(2) model and of its abelian twin from three initial product states:
- on the 12-plaquette torus ($3\times2$), at couplings $g\in\{1.1,1.25,1.4\}$;
- on the 24-plaquette torus ($4\times3$), at $g=1.25$;
- as an open $3\times4$ "string" at $g=1.4$.

Compare the late-time values with the **thermal anchor** and the **decohered values**, flag observables that cannot be told apart from decoherence, and recommend an initial state and coupling for the frontier run.

## Definitions
- **Thermal anchor.** A closed system that thermalizes approaches the canonical average $\mathrm{Tr}(e^{-\beta H}O)/Z$ at the inverse temperature $\beta$ fixed by the conserved energy, $\langle H\rangle_\beta=E(0)$. $\beta<0$ means a negative temperature (energy above the infinite-temperature value).
- **Decohered (infinite-temperature) value.** The value in the maximally mixed state, which is what a fully depolarized quantum computer outputs. $\langle E^2_{pq}\rangle=3/8$, $\langle H_E\rangle_\infty=\tfrac{g^2}{2}\tfrac38 n_\text{bonds}$, and $0$ for $H_B$, $Z_p$, connected $ZZ$ and the hexagon string.
- **Flag.** An observable whose late-time mean (last 25% of the time window) lies within 10% of its decohered value. For observables whose decohered value is 0, "10%" means 10% of the largest $|O(t)|$ along the run. A flagged observable cannot distinguish a working quantum computer from a fully decohered one.

## What was done
- **N = 12.** Sparse $H$, `expm_multiply` evolution on 41 times in $t\in[0,4g^2]$, and $\beta$ from full diagonalization.
- **N = 24.** Matrix-free $H$ (numba). Chebyshev propagation on 25 times. Thermal values from quantum typicality: 4 random vectors with a Chebyshev $e^{-\beta H/2}$, and $\beta$ found by bracketed false position in at most 8 evaluations.
- **Energy conservation.** Total energy is conserved to $<10^{-7}$ in every run; $\sim10^{-14}$ at N = 12. The thermal energies match $E(0)$ within the 2% tolerance; at N = 24 the worst mismatch is 0.002.
- **Outputs.** Data in `results/m3/summary.json` and `results/m3/string_open_3x4_g1.4.json`. Figures in `results/m3/figures/`: energies at each $(N,g)$, finite size, hexagon string and $ZZ$, and the open string. The batch is reproduced by `scripts/run_m3_batch.sh`.

### Electric energy $\langle H_E\rangle$: start, late-time mean, thermal anchor, decohered value

| N | g | model | state | $E_E(0)$ | late | thermal | decohered | $\beta$ |
|---|---|---|---|---|---|---|---|---|
| 12 | 1.1 | SU(2) | vacuum | 0 | 3.146 | 2.883 | 4.084 | 1.100 |
| 12 | 1.1 | SU(2) | neel_half | 2.723 | 3.857 | 3.662 | 4.084 | 0.452 |
| 12 | 1.1 | SU(2) | stripe | 5.445 | 4.646 | 4.508 | 4.084 | −0.460 |
| 12 | 1.25 | SU(2) | vacuum | 0 | 2.178 | 2.488 | 5.273 | 1.448 |
| 12 | 1.25 | SU(2) | neel_half | 3.516 | 4.385 | 4.298 | 5.273 | 0.618 |
| 12 | 1.25 | SU(2) | stripe | 7.031 | 6.447 | 6.263 | 5.273 | −0.632 |
| 12 | 1.4 | SU(2) | vacuum | 0 | 1.773 | 1.710 | 6.615 | 1.624 |
| 12 | 1.4 | SU(2) | neel_half | 4.410 | 5.106 | 4.939 | 6.615 | 0.662 |
| 12 | 1.4 | SU(2) | stripe | 8.820 | 8.405 | 8.307 | 6.615 | −0.670 |
| 12 | 1.25 | abelian | vacuum | 0 | 3.147 | 4.009 | 5.273 | 0.981 |
| 12 | 1.25 | abelian | neel_half | 3.516 | 4.524 | 4.854 | 5.273 | 0.277 |
| 12 | 1.25 | abelian | stripe | 7.031 | 6.023 | 5.693 | 5.273 | −0.277 |
| 24 | 1.25 | SU(2) | vacuum | 0 | 4.466 | 5.154 | 10.547 | 1.440 |
| 24 | 1.25 | SU(2) | neel_half | 10.547 | 10.818 | 10.547 | 10.547 | 0.000 |
| 24 | 1.25 | SU(2) | stripe | 10.547 | 10.662 | 10.548 | 10.547 | 0.000 |
| 24 | 1.25 | abelian | vacuum | 0 | 6.292 | 8.129 | 10.547 | 0.998 |
| 24 | 1.25 | abelian | neel_half | 10.547 | 10.547 | 10.546 | 10.547 | 0.000 |
| 24 | 1.25 | abelian | stripe | 10.547 | 10.547 | 10.547 | 10.547 | 0.000 |

(The abelian rows at $g=1.1$ and $1.4$ are in the JSON.) Open string at $g=1.4$: $\langle Z\rangle$ on the six string plaquettes goes $-1\to-0.78\ (t=1.96)\to-0.46\ (3.92)\to-0.26\ (7.84)$; off the string it goes $+1\to+0.35$. The string melts but has not dissolved by $t=4g^2$.

## What it means
1. **SU(2) thermalizes faster than its abelian twin.** At N = 12 the SU(2) late-time electric energies sit within 2–9% of the thermal anchor. The abelian twin is often 10–25% away, for example from the vacuum: 3.15 against 4.01 at N = 12 and 6.29 against 8.13 at N = 24. The neighbour-dependent flip amplitudes $c(n)=1,\tfrac12$ make the SU(2) dynamics more strongly interacting. This is the physics contrast the frontier run is designed to see.
2. **Finite-size check.** The vacuum $\beta$ hardly changes from N = 12 to N = 24 (1.448 against 1.440), so the thermal anchor is already close to its large-$N$ value.
3. **Important: stripe and Néel-half sit at infinite temperature on even-$L_1$ tori.** Any initial state in which exactly half the bonds are cut has $E_E(0)=\langle H_E\rangle_\infty$ and $E_B(0)=0$, so $\beta=0$. This holds for both states on the $4\times3$ torus and for the stripe on the $4\times6$ frontier torus (N = 48). Their thermal values *are* the decohered values, and the electric energy is flagged in all four of those N = 24 runs. For the abelian twin it is even exact: $U=\prod_{p\in A}X_p\prod_pZ_p$ maps $H\to2\langle H_E\rangle_\infty-H$ and maps these states to translates of themselves, so $\langle H_E\rangle(t)$ is exactly constant (confirmed numerically).
4. **Recommendation for the frontier run.**
   - The SU(2) vacuum at $g\approx1.25$–$1.4$ is the best candidate. It has a finite positive temperature ($\beta\approx1.4$–$1.6$), a large dynamical signal (its $E_E$ rises from 0 to about 40–45% of the decohered value), no flags, and a clear SU(2)/abelian difference.
   - If an imbalanced state is preferred, the N = 12 stripe works ($\beta\approx-0.6$), but on tori with even $L_1$ (N = 24, 48) it must be replaced by a state with an odd number of excited rows.
   - PLAN M6 prescribes the stripe at N = 48. That run is still useful for the MPS convergence study, because entanglement growth is what limits MPS. It is not a good *physics* target, because its late-time value cannot be told apart from decoherence. This is listed as a decision for Digonto.

## Problems / assumptions
- At N = 24, $\langle Z_p\rangle$ is not stored (the schema makes it optional) to keep `summary.json` small. Link energies, the $ZZ$ on bond 0 and the hexagon string are stored.
- METHODS §5 specifies `expm_multiply` for typicality. At N = 24 it cost about 186 s per vector per $\beta$, so a Chebyshev expansion of $e^{-\beta H/2}$ was used instead. It is identical to $10^{-9}$ in a unit test.
- The flag rule (last 25% of the window, 10%) is a choice. With $t_\text{max}=4g^2$ the slow abelian runs have not fully relaxed.
