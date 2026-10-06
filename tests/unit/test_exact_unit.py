import numpy as np
from scipy.sparse.linalg import aslinearoperator
from su2torus.lattice import HoneycombTorus
from su2torus import exact, hamiltonian as hm


def test_chebyshev_matches_expm_multiply_both_directions():
    lat = HoneycombTorus(3, 2)
    H = hm.build_sparse(lat, 1.4)
    psi0 = exact.product_state(lat, [0, 4, 8])
    times = np.array([0.0, 0.7, 3.1, 7.84])
    a = exact.evolve(aslinearoperator(H), psi0, times)
    b = exact.evolve(H, psi0, times)
    assert np.abs(a - b).max() < 1e-10
    back = exact.chebyshev_step(H.dot, a[-1], -7.84, exact.spectral_bounds(H))
    assert np.abs(back - psi0).max() < 1e-10


def test_fused_recurrence_matches_matvec():
    lat = HoneycombTorus(3, 2)
    for ab in (False, True):
        H = hm.build_sparse(lat, 1.25, ab)
        bd = exact.spectral_bounds(H)
        psi0 = exact.product_state(lat, [1, 2, 9])
        a = exact.chebyshev_step(hm.cheb_recur(lat, 1.25, ab), psi0, 0.9, bd)
        b = exact.chebyshev_step(H.dot, psi0, 0.9, bd)
        assert np.abs(a - b).max() < 1e-12
        from su2torus import thermal
        x = thermal.chebyshev_imag(hm.cheb_recur(lat, 1.25, ab), psi0, -0.4, bd)
        y = thermal.chebyshev_imag(H.dot, psi0, -0.4, bd)
        assert np.abs(x - y).max() < 1e-12 * np.abs(y).max()
