"""Noisy emulation and error mitigation (PLAN M4; API fixed by tests/acceptance/m4).

Noise model (Aer convention, k qubits, d = 2**k): after each gate
    E(rho) = (1 - p) rho + p Tr_k(rho) (x) I/d,
with p = p2 on every cx / rzz and p = p1 on every 1-qubit gate, plus a symmetric readout
error r on every qubit. Predicted global fidelity: the probability that no gate applies a
non-identity Pauli, F = prod_gates (1 - eps), eps = p (1 - 1/d^2) (process infidelity):
eps_2q = 15 p2 / 16, eps_1q = 3 p1 / 4.
"""
from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error

from .trotter import ALLOWED_GATES

ONE_QUBIT_GATES = ALLOWED_GATES - {"cx", "rzz", "barrier"}
TWO_QUBIT_GATES = {"cx", "rzz"}


def noise_model(p2, p1=3e-5, readout=1e-3) -> NoiseModel:
    """Aer NoiseModel: depolarizing(p2) on cx/rzz, depolarizing(p1) on 1-qubit gates,
    symmetric readout error. Zero rates add nothing."""
    nm = NoiseModel(basis_gates=sorted(ALLOWED_GATES - {"barrier"}))
    if p2 > 0:
        nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), sorted(TWO_QUBIT_GATES))
    if p1 > 0:
        nm.add_all_qubit_quantum_error(depolarizing_error(p1, 1), sorted(ONE_QUBIT_GATES))
    if readout > 0:
        nm.add_all_qubit_readout_error(ReadoutError([[1 - readout, readout],
                                                     [readout, 1 - readout]]))
    return nm


def _readout_rate(noise) -> float:
    """Symmetric readout flip rate of a noise model built by noise_model (0 if none)."""
    err = getattr(noise, "_default_readout_error", None)
    if err is None:
        return 0.0
    return float(np.asarray(err.probabilities)[0, 1])


def _apply_readout(probs, n, r):
    if r == 0:
        return probs
    M = np.array([[1 - r, r], [r, 1 - r]])         # M[measured, true]
    t = probs.reshape((2,) * n)
    for ax in range(n):
        t = np.moveaxis(np.tensordot(M, t, axes=([1], [ax])), 0, ax)
    return t.reshape(-1)


def run_noisy(circ, noise=None, shots=None, seed=None) -> np.ndarray:
    """Distribution of b = sum_p n_p 2^p after measuring all qubits (readout included).
    shots=None: exact (statevector if ideal, else Aer density matrix + analytic readout).
    shots=int: frequencies from Aer sampling with seed_simulator=seed."""
    n = circ.num_qubits
    ideal = noise is None or noise.is_ideal()
    if shots is None:
        if ideal:
            return np.abs(Statevector.from_instruction(circ).data) ** 2
        qc = circ.copy()
        qc.save_probabilities()
        sim = AerSimulator(method="density_matrix", noise_model=noise)
        res = sim.run(qc, shots=1, seed_simulator=seed).result()
        probs = np.asarray(res.data(0)["probabilities"], dtype=float)
        return _apply_readout(probs, n, _readout_rate(noise))
    qc = circ.copy()
    qc.measure_all()
    sim = AerSimulator(noise_model=None if ideal else noise)
    counts = sim.run(qc, shots=shots, seed_simulator=seed).result().get_counts()
    out = np.zeros(2 ** n)
    for key, c in counts.items():
        out[int(key.replace(" ", ""), 2)] += c
    return out / shots


def _gate_counts(circ):
    ops = circ.count_ops()
    n2 = sum(v for k, v in ops.items() if k in TWO_QUBIT_GATES)
    n1 = sum(v for k, v in ops.items() if k not in TWO_QUBIT_GATES and k != "barrier")
    return n2, n1


def predicted_fidelity(circ, p2, p1=3e-5) -> float:
    """prod_gates (1 - eps): eps = 15 p2/16 per 2-qubit gate, 3 p1/4 per 1-qubit gate."""
    n2, n1 = _gate_counts(circ)
    return float((1 - 15 * p2 / 16) ** n2 * (1 - 3 * p1 / 4) ** n1)


def echo_circuit(circ, initial):
    qc = QuantumCircuit(circ.num_qubits)
    for p in initial or []:
        qc.x(p)
    qc.compose(circ, inplace=True)
    qc.compose(circ.inverse(), inplace=True)
    return qc


def echo_fidelity(circ, initial, noise, shots=None, seed=None) -> float:
    """Loschmidt echo x(initial) circ circ^-1: return probability P of the initial bitstring;
    global-depolarizing model P = F^2 + (1 - F^2)/D gives F = sqrt((P - 1/D)/(1 - 1/D))."""
    D = 2 ** circ.num_qubits
    P = run_noisy(echo_circuit(circ, initial), noise, shots, seed)[sum(1 << p for p in initial or [])]
    return float(np.sqrt(max(0.0, (P - 1 / D) / (1 - 1 / D))))


def mitigate_rescale(noisy, F, O_inf):
    """O_inf + (noisy - O_inf) / F."""
    return O_inf + (np.asarray(noisy) - O_inf) / F if np.ndim(noisy) else \
        O_inf + (noisy - O_inf) / F


def zne(values, scales=(1, 3, 5)) -> float:
    """Richardson extrapolation to scale 0 (polynomial through all points)."""
    x = np.asarray(scales, dtype=float)
    y = np.asarray(values, dtype=float)
    w = [np.prod([-x[j] / (x[i] - x[j]) for j in range(len(x)) if j != i]) for i in range(len(x))]
    return float(np.dot(w, y))


def fold_gates(circ, factor) -> QuantumCircuit:
    """Per-gate folding G -> G (G^-1 G)^k, factor = 2k + 1 (odd >= 1)."""
    if factor < 1 or factor % 2 != 1:
        raise ValueError("factor must be an odd integer >= 1")
    k = (factor - 1) // 2
    out = QuantumCircuit(*circ.qregs, *circ.cregs)
    for inst in circ.data:
        out.append(inst.operation, inst.qubits, inst.clbits)
        if inst.operation.name == "barrier":
            continue
        inv = inst.operation.inverse()
        for _ in range(k):
            out.append(inv, inst.qubits, inst.clbits)
            out.append(inst.operation, inst.qubits, inst.clbits)
    return out
