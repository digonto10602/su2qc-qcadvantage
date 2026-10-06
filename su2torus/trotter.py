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

from qiskit import QuantumCircuit

ALLOWED_GATES = {"rx", "ry", "rz", "h", "x", "y", "z", "s", "sdg", "sx", "sxdg",
                 "p", "u", "cx", "rzz", "barrier"}


def _electric(qc, lat, g, t):
    """U_E(t): one rzz per bond (the constant in H_E is a global phase)."""
    theta = -2.0 * (3.0 * g ** 2 / 16.0) * t
    for p, q in lat.bonds:
        qc.rzz(theta, p, q)


def _flip(qc, lat, g, p, t, abelian):
    """exp(-i t a_p X_p) = H exp(-i t a_p Z_p) H, as a Gray-code uniformly controlled Rz:
    rz(2 t w_S) while the parity on p holds S, one cx(neighbour -> p) per toggle."""
    if abelian:
        qc.rx(2.0 * t * (-1.0 / g ** 2), p)
        return
    nb = lat.neighbors[p]
    k = len(nb)
    w0 = -(1.0 / g ** 2) / 4.0 ** k
    qc.h(p)
    prev = 0
    for i in range(2 ** k):
        gray = i ^ (i >> 1)
        if i > 0:
            qc.cx(nb[(gray ^ prev).bit_length() - 1], p)
        qc.rz(2.0 * t * w0 * 3.0 ** (k - bin(gray).count("1")), p)
        prev = gray
    if k > 0:                                   # close the Gray cycle: parity back to Z_p
        qc.cx(nb[prev.bit_length() - 1], p)
    qc.h(p)


def _magnetic(qc, lat, g, t, sub, abelian):
    for p in range(lat.n):
        if lat.sublattice[p] == sub:
            _flip(qc, lat, g, p, t, abelian)


def trotter_step(lat, g: float, dt: float, order: int = 1, abelian: bool = False):
    """One Trotter step as a qiskit QuantumCircuit on lat.n qubits."""
    qc = QuantumCircuit(lat.n)
    if order == 1:
        _electric(qc, lat, g, dt)
        _magnetic(qc, lat, g, dt, 0, abelian)
        _magnetic(qc, lat, g, dt, 1, abelian)
    elif order == 2:
        _electric(qc, lat, g, dt / 2)
        _magnetic(qc, lat, g, dt / 2, 0, abelian)
        _magnetic(qc, lat, g, dt, 1, abelian)
        _magnetic(qc, lat, g, dt / 2, 0, abelian)
        _electric(qc, lat, g, dt / 2)
    else:
        raise ValueError("order must be 1 or 2")
    return qc


def trotter_circuit(lat, g: float, dt: float, steps: int, order: int = 1,
                    abelian: bool = False, initial=None):
    """X gates on `initial` (excited plaquettes), then `steps` Trotter steps."""
    qc = QuantumCircuit(lat.n)
    for p in (initial or []):
        qc.x(p)
    step = trotter_step(lat, g, dt, order, abelian)
    for _ in range(steps):
        qc.compose(step, inplace=True)
    return qc


def two_qubit_counts(circ) -> dict:
    """{'native_zz': n_cx + n_rzz, 'cz': n_cx + 2*n_rzz}  (METHODS section 3)."""
    ops = circ.count_ops()
    cx, rzz = ops.get("cx", 0), ops.get("rzz", 0)
    return {"native_zz": cx + rzz, "cz": cx + 2 * rzz}
