"""M4: noisy emulation and error mitigation -- API-level acceptance tests.

New public API: module `su2torus/noise.py`
------------------------------------------
ONE_QUBIT_GATES = ALLOWED_GATES - {"cx", "rzz", "barrier"}      (trotter.ALLOWED_GATES)

noise_model(p2, p1=3e-5, readout=1e-3) -> qiskit_aer.noise.NoiseModel
    Aer `depolarizing_error(p2, 2)` on every `cx` and `rzz`, `depolarizing_error(p1, 1)` on every
    gate in ONE_QUBIT_GATES, and the symmetric `ReadoutError([[1-r, r], [r, 1-r]])` (r = readout)
    on every qubit. Aer convention for k qubits, d = 2**k:
        E(rho) = (1 - p) rho + p Tr_k(rho) (x) I/d,   applied after the gate.
    Zero rates add nothing, so noise_model(0, 0, 0).is_ideal() is True.

run_noisy(circ, noise=None, shots=None, seed=None) -> np.ndarray, shape (2**n,), float
    Distribution of the bitstring b = sum_p n_p 2**p (METHODS section 2) when all qubits of `circ`
    (which has no measurements) are measured at the end, readout error included.
    shots=None: exact distribution (statevector if `noise` is None or ideal, otherwise the
                Aer density-matrix method; readout error applied analytically as the per-bit
                confusion matrix).
    shots=int : empirical frequencies counts/shots from Aer sampling with seed_simulator=seed;
                the same seed must give identical output.

predicted_fidelity(circ, p2, p1=3e-5) -> float
    Global fidelity predicted from the noise model = probability that no gate suffers a
    non-identity Pauli error:  F = prod_gates (1 - eps_gate),  with the process (entanglement)
    infidelity of the Aer depolarizing channel eps = p (1 - 1/d**2):
        eps_2q = 15 p2 / 16  for each cx / rzz,   eps_1q = 3 p1 / 4  for each 1-qubit gate.
    (Readout is not included. The average gate infidelity p (d-1)/d is NOT the convention used.)

echo_fidelity(circ, initial, noise, shots=None, seed=None) -> float
    Forward-backward Loschmidt echo. `circ` is the evolution only (no preparation). The echo
    circuit is: x on every plaquette in `initial`, then circ, then circ.inverse(). With
    P = run_noisy(echo, noise, shots, seed)[b_init], b_init = sum_{p in initial} 2**p, D = 2**n:
        F = sqrt(max(0, (P - 1/D) / (1 - 1/D)))        (global-depolarizing model, F_echo = F**2).

mitigate_rescale(noisy, F, O_inf) -> O_inf + (noisy - O_inf) / F        (numpy-broadcasting)

zne(values, scales=(1, 3, 5)) -> float
    Richardson extrapolation: the polynomial of degree len(values)-1 through (scales, values),
    evaluated at scale 0.

fold_gates(circ, factor) -> QuantumCircuit
    Gate folding for odd integer factor >= 1 (global circ (circ^-1 circ)^k or per-gate G (G^-1 G)^k).
    Same unitary up to a global phase, and the numbers of 2-qubit and of 1-qubit gates are each
    exactly `factor` times those of `circ` (so the noise model scales the error by `factor`).
"""
import functools
import itertools
import numpy as np
import pytest
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from scipy.sparse.linalg import expm_multiply
from su2torus.lattice import HoneycombTorus
from su2torus.trotter import trotter_circuit
from su2torus.noise import (noise_model, run_noisy, predicted_fidelity, echo_fidelity,
                            mitigate_rescale, zne, fold_gates)
from tests.acceptance._ref import ref_neighbors, ref_hamiltonian, basis_vec

G = 1.25
DT = 0.45 * G ** 2          # 2 dt / g^2 = 0.9, as in the frontier run (PLAN M6)
TWO_Q = {"cx", "rzz"}


def stripe(L1, L2):
    return [2 * (i * L2 + j) for i in range(0, L1, 2) for j in range(L2)]


@functools.lru_cache(maxsize=None)
def ref_parts(L1, L2, g):
    return ref_hamiltonian(ref_neighbors(L1, L2, True), g, parts=True)


def ref_trotter_state(L1, L2, g, dt, steps, init):
    """Order-1 Trotter state from brute-force H parts (METHODS section 3); U_E acts first."""
    HE, HA, HB = ref_parts(L1, L2, g)
    dE = HE.diagonal()
    v = basis_vec(2 * L1 * L2, init)
    for _ in range(steps):
        v = np.exp(-1j * dE * dt) * v
        v = expm_multiply(-1j * dt * HA, v)    # commuting terms: exact sublattice exponential
        v = expm_multiply(-1j * dt * HB, v)
    return v, dE


def n_gates(circ):
    ops = circ.count_ops()
    n2 = sum(v for k, v in ops.items() if k in TWO_Q)
    n1 = sum(v for k, v in ops.items() if k not in TWO_Q and k != "barrier")
    return n2, n1


# ---------------------------------------------------------------- noiseless emulation
def test_noiseless_equals_statevector():
    L1, L2 = 3, 2                                   # N = 12 torus
    init = stripe(L1, L2)
    lat = HoneycombTorus(L1, L2)
    circ = trotter_circuit(lat, G, DT, 2, 1, False, init)
    psi, dE = ref_trotter_state(L1, L2, G, DT, 2, init)
    ref = np.abs(psi) ** 2
    assert noise_model(0.0, 0.0, 0.0).is_ideal()
    for noise in (None, noise_model(0.0, 0.0, 0.0)):
        p = np.asarray(run_noisy(circ, noise))
        assert p.shape == (2 ** 12,)
        # 1e-10: PLAN criterion; both sides are double-precision exact evolutions
        assert np.abs(p - ref).max() < 1e-10
        assert abs(p @ dE - ref @ dE) < 1e-10


# ---------------------------------------------------------------- noise model = Aer channel
def _pauli_full(n, qubits, labels):
    P = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]),
         "Z": np.diag([1.0, -1.0])}
    ops = ["I"] * n
    for q, l in zip(qubits, labels):
        ops[q] = l
    M = np.array([[1.0]])
    for q in reversed(range(n)):                    # little endian: qubit 0 is rightmost factor
        M = np.kron(M, P[ops[q]])
    return M


def _depolarize(rho, n, qubits, p):
    k = len(qubits)
    twirl = sum(_pauli_full(n, qubits, ls) @ rho @ _pauli_full(n, qubits, ls).conj().T
                for ls in itertools.product("IXYZ", repeat=k)) / 4 ** k
    return (1 - p) * rho + p * twirl                # = (1-p) rho + p Tr_k(rho) (x) I/2^k


def _brute_force_probs(qc, p2, p1, r):
    n = qc.num_qubits
    rho = np.zeros((2 ** n, 2 ** n), complex)
    rho[0, 0] = 1
    for inst in qc.data:
        qs = [qc.find_bit(q).index for q in inst.qubits]
        one = QuantumCircuit(n)
        one.append(inst.operation, qs)
        U = Operator(one).data
        rho = U @ rho @ U.conj().T
        rho = _depolarize(rho, n, qs, p2 if len(qs) == 2 else p1)
    probs = np.real(np.diag(rho))
    out = np.zeros_like(probs)
    for b in range(2 ** n):
        for m in range(2 ** n):
            flips = bin(b ^ m).count("1")
            out[m] += probs[b] * r ** flips * (1 - r) ** (n - flips)
    return out


def _test_circuit():
    qc = QuantumCircuit(3)
    qc.x(2); qc.h(0); qc.rzz(0.7, 0, 1); qc.cx(1, 2); qc.rz(0.4, 1); qc.h(1)
    qc.cx(0, 1); qc.rx(0.9, 2); qc.rzz(-1.1, 1, 2); qc.h(0); qc.h(2)
    return qc


@pytest.mark.parametrize("p2,p1,r", [(0.04, 0.0, 0.0), (0.0, 0.03, 0.0), (0.0, 0.0, 0.1),
                                     (0.04, 0.02, 0.03)])
def test_noise_model_channel(p2, p1, r):
    qc = _test_circuit()
    ref = _brute_force_probs(qc, p2, p1, r)
    p = np.asarray(run_noisy(qc, noise_model(p2, p1, r)))
    # 1e-10: exact density-matrix calculation on 3 qubits, only round-off remains
    assert np.abs(p - ref).max() < 1e-10


def test_sampling_seeded():
    qc = _test_circuit()
    nm = noise_model(0.04, 0.02, 0.03)
    shots = 20000
    a = np.asarray(run_noisy(qc, nm, shots=shots, seed=11))
    b = np.asarray(run_noisy(qc, nm, shots=shots, seed=11))
    c = np.asarray(run_noisy(qc, nm, shots=shots, seed=12))
    assert np.array_equal(a, b), "same seed must give identical results"
    assert not np.array_equal(a, c), "different seeds should give different samples"
    assert abs(a.sum() - 1) < 1e-12 and np.allclose(a * shots, np.round(a * shots), atol=1e-6)
    # 0.02 > 5 sigma of a frequency at 20000 shots (sigma <= 0.5/sqrt(20000) = 0.0035)
    assert np.abs(a - _brute_force_probs(qc, 0.04, 0.02, 0.03)).max() < 0.02


# ---------------------------------------------------------------- fidelity prediction / echo
def test_predicted_fidelity_formula():
    lat = HoneycombTorus(3, 2)
    for steps, p2 in [(2, 7.9e-4), (10, 1.5e-3)]:
        circ = trotter_circuit(lat, G, DT, steps, 1, False, None)
        n2, n1 = n_gates(circ)
        ref = (1 - 15 * p2 / 16) ** n2 * (1 - 3 * 3e-5 / 4) ** n1
        assert abs(predicted_fidelity(circ, p2, 3e-5) / ref - 1) < 1e-12   # same closed form


def test_echo_fidelity():
    L1, L2 = 2, 2                                   # n = 8: density matrix is cheap
    lat = HoneycombTorus(L1, L2)
    init = stripe(L1, L2)
    circ = trotter_circuit(lat, G, DT, 1, 1, False, None)
    assert abs(echo_fidelity(circ, init, None) - 1) < 1e-10   # noiseless echo returns exactly
    nm = noise_model(1.5e-3)
    echo = QuantumCircuit(lat.n)
    for p in init:
        echo.x(p)
    echo.compose(circ, inplace=True)
    echo.compose(circ.inverse(), inplace=True)
    P = np.asarray(run_noisy(echo, nm))[sum(1 << p for p in init)]
    D = 2 ** lat.n
    ref = np.sqrt(max(0.0, (P - 1 / D) / (1 - 1 / D)))
    assert abs(echo_fidelity(circ, init, nm) - ref) < 1e-10   # same exact distribution
    assert 0 < ref < 1


# ---------------------------------------------------------------- mitigation primitives
def test_mitigate_rescale():
    assert abs(mitigate_rescale(5.9, 0.5, 5.0) - 6.8) < 1e-12
    x = mitigate_rescale(np.array([1.0, 2.0]), 0.25, 3.0)
    assert np.allclose(x, [3.0 - 8.0, 3.0 - 4.0], atol=1e-12)


def test_zne_richardson():
    # Lagrange weights at 0 for nodes 1, 3, 5: 15/8, -5/4, 3/8
    v = [2.0, 1.3, 0.7]
    assert abs(zne(v) - (15 / 8 * v[0] - 5 / 4 * v[1] + 3 / 8 * v[2])) < 1e-12
    a, b, c = 0.4, -0.7, 0.05
    f = lambda s: a + b * s + c * s * s
    assert abs(zne([f(1), f(3), f(5)], (1, 3, 5)) - a) < 1e-12    # exact for a quadratic
    assert abs(zne([f(1), f(2)], (1, 2)) - (2 * f(1) - f(2))) < 1e-12


@pytest.mark.parametrize("factor", [1, 3, 5])
def test_fold_gates(factor):
    lat = HoneycombTorus(2, 2)
    circ = trotter_circuit(lat, G, DT, 1, 1, False, stripe(2, 2))
    folded = fold_gates(circ, factor)
    n2, n1 = n_gates(circ)
    assert n_gates(folded) == (factor * n2, factor * n1)
    A, B = Operator(folded).data, Operator(circ).data
    k = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    phase = A[k] / B[k]
    # 1e-9: products of ~10^3 exact gate matrices, round-off only
    assert abs(abs(phase) - 1) < 1e-9 and np.abs(A - phase * B).max() < 1e-9
