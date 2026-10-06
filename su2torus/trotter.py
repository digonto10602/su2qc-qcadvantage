"""Trotter circuits (docs/METHODS.md, section 3).

Step ordering, with U_X(t) = exp(-i H_X t):
    order 1:  U_step(dt) = U_BB(dt) U_BA(dt) U_E(dt)      (U_E acts first)
    order 2:  U_step(dt) = U_E(dt/2) U_BA(dt/2) U_BB(dt) U_BA(dt/2) U_E(dt/2)
Equality is required up to a global phase.

Allowed gates: rx ry rz h x y z s sdg sx sxdg p u cx rzz (and barrier).
Each flip exp(-i a_p(Z_nb) X_p dt) is a uniformly controlled rotation on qubit p
controlled by its neighbours (Gray-code construction: 2**(#neighbours) cx gates).
"""
from __future__ import annotations

ALLOWED_GATES = {"rx", "ry", "rz", "h", "x", "y", "z", "s", "sdg", "sx", "sxdg",
                 "p", "u", "cx", "rzz", "barrier"}


def trotter_step(lat, g: float, dt: float, order: int = 1, abelian: bool = False):
    """One Trotter step as a qiskit QuantumCircuit on lat.n qubits."""
    raise NotImplementedError


def trotter_circuit(lat, g: float, dt: float, steps: int, order: int = 1,
                    abelian: bool = False, initial=None):
    """X gates on `initial` (excited plaquettes), then `steps` Trotter steps."""
    raise NotImplementedError


def two_qubit_counts(circ) -> dict:
    """{'native_zz': n_cx + n_rzz, 'cz': n_cx + 2*n_rzz}  (METHODS section 3)."""
    raise NotImplementedError
