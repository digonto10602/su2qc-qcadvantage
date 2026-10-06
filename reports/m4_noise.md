# M4 — noisy emulation and error mitigation (N = 12 torus)

## What was asked
Emulate the first-order Trotter circuits of the 12-plaquette torus under realistic hardware noise, and test two error-mitigation methods: Loschmidt-echo rescaling and zero-noise extrapolation. The test checks that the echo estimates the circuit fidelity, and that mitigation at least halves the error of the electric energy $\langle H_E\rangle$ in at least 3 of 4 configurations.

## Definitions
- **Depolarizing noise.** After every gate on $k$ qubits ($d=2^k$), the state is replaced by $(1-p)\rho+p\,\mathrm{Tr}_k(\rho)\otimes I/d$. Here $p=p_2$ for two-qubit gates (`cx`, `rzz`) and $p_1=3\times10^{-5}$ for one-qubit gates. In addition, each measured bit is flipped with probability $10^{-3}$ (readout error).
- **Global fidelity $F$.** The probability that the whole circuit runs without any error. It is predicted from the noise model as
  $$F_\text{pred}=\prod_\text{gates}(1-\epsilon_\text{gate}),\qquad \epsilon=p\,(1-1/d^2):\quad \epsilon_{2q}=\tfrac{15}{16}p_2,\ \epsilon_{1q}=\tfrac34 p_1 .$$
  $\epsilon$ is the process infidelity of the Aer depolarizing channel: the probability that it applies a non-identity Pauli operator. Readout error is not included in $F_\text{pred}$.
- **Loschmidt echo.** Run the circuit forward and then exactly backward ($U^\dagger U$), and measure the probability $P$ of returning to the initial bit string. Under a global-depolarizing model, $P=F^2+(1-F^2)/D$ with $D=2^{12}$, so $F_\text{echo}=\sqrt{(P-1/D)/(1-1/D)}$.
- **Echo rescaling.** $\langle O\rangle_\text{mit}=O_\infty+(\langle O\rangle_\text{noisy}-O_\infty)/F$, where $O_\infty$ is the decohered value ($\langle H_E\rangle_\infty=5.2734$).
- **Zero-noise extrapolation (ZNE).** Every gate $G$ is replaced by $G(G^{-1}G)^k$ (*gate folding*), which multiplies the noise by a factor 1, 3 or 5. A quadratic is then fitted through the three noisy values and extrapolated to zero noise: $\langle O\rangle_0=\tfrac{15}{8}O_1-\tfrac54O_3+\tfrac38O_5$ (Richardson extrapolation).

## What was done
- `su2torus/noise.py` implements the noise model, exact noisy emulation (Aer density matrix, with readout applied analytically), seeded shot sampling, the echo, folding and ZNE.
- `scripts/run_m4.py` runs the stripe state at $g=1.25$ with $\delta t=0.45g^2$ ($2\delta t/g^2=0.9$, the frontier-run choice).
- The reference ("ideal") is the noiseless Trotter circuit, because mitigation targets hardware noise, not Trotter error.
- All configurations are exact density-matrix calculations, so they have no shot noise. A seeded 4000-shot run (`repro` block) shows bit-for-bit reproducibility.
- Data: `results/m4/noise_runs.json`. Figure: `results/m4/figures/noise_mitigation.png`.

| steps | $p_2$ | 2q gates | ideal | raw (error) | echo-rescaled (error) | ZNE (error) | $F_\text{pred}$ | $F_\text{echo}$ |
|---|---|---|---|---|---|---|---|---|
| 2 | $7.9\times10^{-4}$ | 228 | 6.1378 | 6.1019 (0.036) | 6.2332 (0.095) | 6.1343 (0.0035) | 0.840 | 0.863 |
| 4 | $1.0\times10^{-3}$ | 456 | 6.5661 | 6.4246 (0.141) | 6.9514 (0.385) | 6.5569 (0.0091) | 0.645 | 0.686 |
| 6 | $1.5\times10^{-3}$ | 684 | 6.3286 | 6.0836 (0.245) | 7.1750 (0.846) | 6.2925 (0.036) | 0.376 | 0.426 |
| 10 | $1.0\times10^{-3}$ | 1140 | 6.2695 | 6.0087 (0.261) | 7.2613 (0.992) | 6.2245 (0.045) | 0.334 | 0.370 |

## What it means
- **The echo measures the fidelity.** $F_\text{echo}$ is within 3–13% of $F_\text{pred}$; the criterion allows 25%. The echo is always slightly higher than the prediction, because some errors (e.g. $Z$ errors on qubits in a $Z$ eigenstate) do not change the returned bit string.
- **ZNE meets the criterion in 4 of 4 configurations.** It reduces the error by a factor of 5.8–15.5, well beyond the required factor of 2. ZNE is the method recorded in the results file (`mitigation_method: "zne"`).
- **Global echo rescaling fails.** It makes the error 2.7–3.8 times *larger*. The noisy $\langle H_E\rangle$ moves toward its decohered value much more slowly than $F$ decays. A single depolarizing error damages only a few local link energies; it does not scramble the whole state. Dividing by the global $F$ therefore overcorrects. Global-fidelity rescaling is not suitable for local observables such as link energies; observable-specific decay factors or ZNE should be used instead.
- **Cost of one exact emulation.** At 10 steps (1140 two-qubit gates, $F\approx0.33$), raw noise biases $\langle H_E\rangle$ by about $0.26$, roughly 25% of the distance between the ideal value and the decohered value. Folding by 5 still leaves a signal that can be extrapolated.

## Problems / assumptions
- Only the four required (steps, $p_2$) points were run. Each density-matrix configuration costs 6–45 minutes on 4 vCPUs (ZNE needs folds 1, 3, 5 plus the echo), and the full grid would have taken about 4.5 h. N = 16 was not run.
- Exact (infinite-shot) values are used for mitigation. With hardware shot budgets, ZNE amplifies statistical noise by $\sqrt{(15/8)^2+(5/4)^2+(3/8)^2}\approx2.3$.
- The noise is purely depolarizing and Markovian. Coherent errors, crosstalk and leakage are not modelled.
