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
| 12 | 1.1 | abelian | vacuum | 0 | 2.071 | 3.682 | 4.084 | 0.477 |
| 12 | 1.1 | abelian | neel_half | 2.723 | 3.418 | 3.946 | 4.084 | 0.150 |
| 12 | 1.1 | abelian | stripe | 5.445 | 4.749 | 4.222 | 4.084 | −0.150 |
| 12 | 1.25 | abelian | vacuum | 0 | 3.147 | 4.009 | 5.273 | 0.981 |
| 12 | 1.25 | abelian | neel_half | 3.516 | 4.524 | 4.854 | 5.273 | 0.277 |
| 12 | 1.25 | abelian | stripe | 7.031 | 6.023 | 5.693 | 5.273 | −0.277 |
| 12 | 1.4 | abelian | vacuum | 0 | 3.076 | 3.067 | 6.615 | 1.567 |
| 12 | 1.4 | abelian | neel_half | 4.410 | 5.588 | 5.640 | 6.615 | 0.408 |
| 12 | 1.4 | abelian | stripe | 8.820 | 7.642 | 7.590 | 6.615 | −0.408 |
| 24 | 1.25 | SU(2) | vacuum | 0 | 4.466 | 5.154 | 10.547 | 1.440 |
| 24 | 1.25 | SU(2) | neel_half | 10.547 | 10.818 | 10.547 | 10.547 | 0.000 |
| 24 | 1.25 | SU(2) | stripe | 10.547 | 10.662 | 10.548 | 10.547 | 0.000 |
| 24 | 1.25 | abelian | vacuum | 0 | 6.292 | 8.129 | 10.547 | 0.998 |
| 24 | 1.25 | abelian | neel_half | 10.547 | 10.547 | 10.546 | 10.547 | 0.000 |
| 24 | 1.25 | abelian | stripe | 10.547 | 10.547 | 10.547 | 10.547 | 0.000 |

N = 24 thermal values are 4-vector typicality estimates with a statistical uncertainty of about $\pm0.2$ in $\langle H_E\rangle$ and $\pm0.04$ in $\beta$ (audit re-run with 8 fresh vectors: $\langle H_E\rangle=5.27\pm0.16$ for the SU(2) vacuum). Open string at $g=1.4$: $\langle Z\rangle$ on the six string plaquettes goes $-1\to-0.78\ (t=1.96)\to-0.46\ (3.92)\to-0.26\ (7.84)$; off the string it goes $+1\to+0.35$. The string melts but has not dissolved by $t=4g^2$.

## What it means
1. **SU(2) versus abelian relaxation depends on the coupling.** The table below gives the late-time distance from the thermal anchor, $|E_\text{late}-E_\text{th}|/E_\text{th}$, at N = 12:

   | g | SU(2) vacuum / neel_half / stripe | abelian vacuum / neel_half / stripe |
   |---|---|---|
   | 1.1 | 9.1% / 5.3% / 3.1% | 43.8% / 13.4% / 12.5% |
   | 1.25 | 12.5% / 2.0% / 2.9% | 21.5% / 6.8% / 5.8% |
   | 1.4 | 3.7% / 3.4% / 1.2% | 0.3% / 0.9% / 0.7% |

   At $g\le1.25$ the abelian twin is much further from its anchor than SU(2) by $t=4g^2$, most clearly from the vacuum (at N = 24: 6.29 against 8.13, while SU(2) is 13%±4% below its anchor). At $g=1.4$ both are within about 4%, and the abelian twin is closer. "SU(2) relaxes faster" therefore holds only for $g\lesssim1.25$ and a time window of $4g^2$.
2. **Finite-size check.** The vacuum $\beta$ is $1.448$ (N = 12, exact) and $1.44\pm0.04$ (N = 24, typicality). These agree within the statistical error. The data cannot resolve a finite-size shift smaller than that.
3. **Important: stripe and Néel-half sit at infinite temperature on even-$L_1$ tori.** Any initial state in which exactly half the bonds are cut has $E_E(0)=\langle H_E\rangle_\infty$ and $E_B(0)=0$, so $\beta=0$. This holds for both states on the $4\times3$ torus and for the stripe on the $4\times6$ frontier torus (N = 48). Their thermal values *are* the decohered values, and the electric energy is flagged in all four of those N = 24 runs. For the abelian twin it is even exact: $U=\prod_{p\in A}X_p\prod_pZ_p$ maps $H\to2\langle H_E\rangle_\infty-H$ and maps these states to translates of themselves. Translations commute with $H$, and $H$ is real (time-reversal symmetric) while the initial states are real, so $\langle H_E\rangle(t)=2\langle H_E\rangle_\infty-\langle H_E\rangle(t)$, i.e. $\langle H_E\rangle(t)$ is exactly constant (confirmed numerically).
4. **Recommendation for the frontier run.**
   - The SU(2) vacuum at $g\approx1.25$–$1.4$ is the best candidate. It has a finite positive temperature ($\beta\approx1.4$–$1.6$), a large dynamical signal (its late-time $E_E$ is 77%, 41% and 27% of the decohered value at $g=1.1$, 1.25 and 1.4 for N = 12, and 42% at N = 24, $g=1.25$), no flags, and a clear SU(2)/abelian difference.
   - If an imbalanced state is preferred, the N = 12 stripe works ($\beta\approx-0.6$), but on tori with even $L_1$ (N = 24, 48) it must be replaced by a state with an odd number of excited rows.
   - PLAN M6 prescribes the stripe at N = 48. That run is still useful for the MPS convergence study, because entanglement growth is what limits MPS. It is not a good *physics* target, because its late-time value cannot be told apart from decoherence. This is listed as a decision for Digonto.

## Problems / assumptions
- At N = 24, $\langle Z_p\rangle$ is not stored (the schema makes it optional) to keep `summary.json` small. Link energies, the $ZZ$ on bond 0 and the hexagon string are stored.
- METHODS §5 specifies `expm_multiply` for typicality. At N = 24 it cost about 186 s per vector per $\beta$, so a Chebyshev expansion of $e^{-\beta H/2}$ was used instead. It is identical to $10^{-9}$ in a unit test. No error bar is computed in the code; the audit's 8-vector re-run gives the uncertainty quoted above.
- The flag rule (last 25% of the window, 10%) is a choice. With $t_\text{max}=4g^2$ the slow abelian runs have not fully relaxed.
