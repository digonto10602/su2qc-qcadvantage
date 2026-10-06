# su2qc-qcadvantage

Code toward a beyond-classical quantum simulation of minimally truncated (2+1)D SU(2) lattice gauge theory on tori: the "SU(2) torus frontier run". The idea, the hardware case and the references are in [`docs/idea.md`](docs/idea.md). The binding specification is [`docs/METHODS.md`](docs/METHODS.md), the milestones are in [`plan/PLAN.md`](plan/PLAN.md), and progress is tracked in [`notes/STATUS.md`](notes/STATUS.md).

**Status:** cloud build phase (M1–M6). This phase stops before the GPU-cluster phase.

```bash
python3 -m pip install -r requirements.txt
python3 -m pytest -q          # acceptance tests (m1-m3 locked) + unit tests
```

The physics model is the one-qubit-per-plaquette SU(2) Hamiltonian of Ciavarella, de Putter, Younis and Rrapaj, arXiv:2608.28752, studied here on periodic honeycomb plaquette graphs.
