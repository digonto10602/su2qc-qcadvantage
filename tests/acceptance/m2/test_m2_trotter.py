"""M2 / G6-G8: Trotter circuits and gate counts (LOCKED)."""
import numpy as np
import pytest
import scipy.linalg as sla
from qiskit.quantum_info import Operator, Statevector
from su2torus.lattice import HoneycombTorus
from su2torus.trotter import trotter_step, trotter_circuit, two_qubit_counts, ALLOWED_GATES
from tests.acceptance._ref import ref_hamiltonian, basis_vec


def ref_step(lat, g, dt, order, abelian):
    HE, HA, HB = (m.toarray() for m in ref_hamiltonian(lat.neighbors, g, abelian, parts=True))
    U = lambda H, t: sla.expm(-1j * H * t)
    if order == 1:
        return U(HB, dt) @ U(HA, dt) @ U(HE, dt)
    return U(HE, dt / 2) @ U(HA, dt / 2) @ U(HB, dt) @ U(HA, dt / 2) @ U(HE, dt / 2)


def equal_up_to_phase(A, B, tol=1e-9):
    k = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    phase = A[k] / B[k]
    return abs(abs(phase) - 1) < tol and np.abs(A - phase * B).max() < tol


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("abelian", [False, True])
def test_step_unitary(order, abelian):                          # G7
    lat = HoneycombTorus(2, 2)
    g, dt = 1.25, 0.37
    qc = trotter_step(lat, g, dt, order, abelian)
    assert qc.num_qubits == lat.n
    assert set(qc.count_ops()) <= ALLOWED_GATES
    assert equal_up_to_phase(Operator(qc).data, ref_step(lat, g, dt, order, abelian))


def test_circuit_with_initial_state():                          # G7
    lat = HoneycombTorus(2, 2)
    g, dt, steps = 1.1, 0.3, 3
    init = [0, 3, 6]
    sv = Statevector.from_instruction(trotter_circuit(lat, g, dt, steps, 1, False, init)).data
    U = ref_step(lat, g, dt, 1, False)
    ref = np.linalg.matrix_power(U, steps) @ basis_vec(lat.n, init)
    assert abs(abs(np.vdot(ref, sv)) - 1) < 1e-9


@pytest.mark.parametrize("order,lo,hi", [(1, 0.8, 1.3), (2, 1.7, 2.4)])
def test_trotter_error_scaling(order, lo, hi):                  # G6
    lat = HoneycombTorus(2, 2)
    g, T = 1.25, 1.0
    H = ref_hamiltonian(lat.neighbors, g).toarray()
    init = [0, 2]
    exact_state = sla.expm(-1j * H * T) @ basis_vec(lat.n, init)
    dts, errs = [0.2, 0.1, 0.05], []
    for dt in dts:
        sv = Statevector.from_instruction(
            trotter_circuit(lat, g, dt, int(round(T / dt)), order, False, init)).data
        errs.append(np.sqrt(max(0.0, 2 - 2 * abs(np.vdot(exact_state, sv)))))
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    assert lo < slope < hi, f"slope {slope:.3f}, errors {errs}"


def test_gate_counts():                                         # G8
    lat = HoneycombTorus(3, 2)
    c = two_qubit_counts(trotter_step(lat, 1.25, 0.3, 1, False))
    ops = trotter_step(lat, 1.25, 0.3, 1, False).count_ops()
    assert c["native_zz"] == ops.get("cx", 0) + ops.get("rzz", 0)
    assert c["cz"] == ops.get("cx", 0) + 2 * ops.get("rzz", 0)
    assert c["native_zz"] / lat.n <= 9.5 + 1e-12, c
    assert c["cz"] / lat.n <= 11 + 1e-12, c
    ca = two_qubit_counts(trotter_step(lat, 1.25, 0.3, 1, True))
    assert ca["native_zz"] / lat.n <= 1.5 + 1e-12, ca
