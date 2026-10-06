# su2qc-qcadvantage

Code toward a beyond-classical quantum simulation of minimally truncated (2+1)D SU(2) lattice gauge theory on tori: the "SU(2) torus frontier run". The idea, the hardware case and the references are in [`docs/idea.md`](docs/idea.md). The binding specification is [`docs/METHODS.md`](docs/METHODS.md), the milestones are in [`plan/PLAN.md`](plan/PLAN.md), and progress is tracked in [`notes/STATUS.md`](notes/STATUS.md).

**Status:** cloud phase complete (M1–M6, tag `v0.1-pre-cluster`). Main findings: [`reports/collaborator_summary.md`](reports/collaborator_summary.md). Full report with plotting data: [`reports/project_report.md`](reports/project_report.md).

```bash
python3 -m pip install -r requirements.txt
python3 -m pytest -q          # acceptance tests (m1-m6 locked) + unit tests
```

The physics model is the one-qubit-per-plaquette SU(2) Hamiltonian of Ciavarella, de Putter, Younis and Rrapaj, arXiv:2608.28752, studied here on periodic honeycomb plaquette graphs.

## Layout
- `su2torus/`: the package. `lattice` (honeycomb plaquette graph), `hamiltonian` (sparse, Pauli and matrix-free numba H), `exact` (expm_multiply/Chebyshev evolution), `trotter` (circuits, gate counts), `observables`, `thermal` (exact and typicality anchors), `noise` (Aer noise, echo, ZNE), `mps` (quimb CircuitMPS).
- `scripts/`: the run scripts below. `scripts/slurm/`: cluster templates with placeholders only.
- `configs/`: JSON configs for `scripts/run_classical.py`.
- `results/m2`–`results/m6`: JSON data and figures.
- `reports/`: milestone reports, audits and the collaborator summary.

## Reproducing everything
Tested on 4 vCPUs and 16 GB. Times are wall-clock.

```bash
python3 -m pip install -r requirements.txt
python3 scripts/gate_counts.py                 # M2 -> results/m2/gate_counts.json (seconds)
nohup scripts/run_m3_batch.sh > logs/m3.log 2>&1 &   # M3: run_m3 n12 / string / n24 x6 / assemble -> results/m3 (~4 h)
python3 scripts/run_m4.py                      # M4 -> results/m4/noise_runs.json + figure (~1.5 h)
OPENBLAS_NUM_THREADS=1 python3 scripts/run_m5_check.py      # M5 -> results/m5/mps_check.json (~1 min)
python3 scripts/resource_table.py              # M5 -> results/m5/resources.json
OPENBLAS_NUM_THREADS=1 python3 scripts/run_classical.py --config configs/n48_stripe_g1.25.json   # M6 -> results/m6/n48_stripe_g1.25.json (~2.3 h)
python3 scripts/analyze_convergence.py         # M6 -> results/m6/convergence.json
python3 scripts/make_figures.py                # reports/figures/fig1..fig5
python3 scripts/check_summary.py               # verifies every number quoted in reports/collaborator_summary.md
python3 -m pytest -q                           # all acceptance + unit tests
```

Notes:
- **Thread contention.** Set `OPENBLAS_NUM_THREADS=1` for MPS runs whenever another multithreaded job (numba, Aer) shares the CPUs. Oversubscribed BLAS made SVDs about 400× slower here.
- **Cluster runs.** Fill the placeholders in `scripts/slurm/mps_cpu.sbatch` or `scripts/slurm/mps_gpu.sbatch` and pass a config such as `configs/n96_stripe_g1.25_gpu.json`.
