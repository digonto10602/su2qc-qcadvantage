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
