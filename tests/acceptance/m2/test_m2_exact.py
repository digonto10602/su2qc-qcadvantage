"""M2 / G6: exact evolution (LOCKED)."""
import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator
from su2torus.lattice import HoneycombTorus
from su2torus import exact, hamiltonian as hm
from tests.acceptance._ref import ref_hamiltonian, basis_vec


def test_product_state():
    lat = HoneycombTorus(2, 2)
    v = exact.product_state(lat, [0, 3])
    assert v.dtype == np.complex128 and v.shape == (256,)
    assert v[(1 << 0) | (1 << 3)] == 1 and np.abs(v).sum() == 1


def test_ground_state():
    lat = HoneycombTorus(2, 2)
    H = ref_hamiltonian(lat.neighbors, 1.25)
    E0, psi = exact.ground_state(H)
    w = np.linalg.eigvalsh(H.toarray())
    assert abs(E0 - w[0]) < 1e-9
    assert abs(np.vdot(psi, H @ psi).real - w[0]) < 1e-8


def test_evolve_matches_dense_expm_and_conserves_energy():
    lat = HoneycombTorus(2, 2)
    H = ref_hamiltonian(lat.neighbors, 1.25)
    psi0 = basis_vec(lat.n, [0, 2, 4])
    times = np.linspace(0, 3, 7)
    out = exact.evolve(H, psi0, times)
    Hd = H.toarray()
    E0 = np.vdot(psi0, Hd @ psi0).real
    for t, psi in zip(times, out):
        assert np.abs(psi - sla.expm(-1j * Hd * t) @ psi0).max() < 1e-9
        assert abs(np.vdot(psi, Hd @ psi).real - E0) < 1e-10
        assert abs(np.linalg.norm(psi) - 1) < 1e-12


def test_evolve_accepts_linear_operator():
    lat = HoneycombTorus(3, 2)
    g = 1.1
    dim = 2 ** lat.n
    mv = lambda v: hm.apply_h(lat, g, v)
    op = LinearOperator((dim, dim), matvec=mv, rmatvec=mv, dtype=complex)  # H Hermitian
    psi0 = basis_vec(lat.n, [0, 6])
    times = np.array([0.0, 0.5, 2.0])
    a = exact.evolve(op, psi0, times)
    b = exact.evolve(ref_hamiltonian(lat.neighbors, g), psi0, times)
    assert np.abs(a - b).max() < 1e-8
