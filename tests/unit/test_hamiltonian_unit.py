import numpy as np
import pytest
from su2torus.lattice import HoneycombTorus
from su2torus import hamiltonian as hm


@pytest.mark.parametrize("abelian", [False, True])
@pytest.mark.parametrize("dims,periodic", [((3, 2), True), ((3, 3), False)])
def test_numpy_fallback_matches_sparse(dims, periodic, abelian):
    lat = HoneycombTorus(*dims, periodic=periodic)
    rng = np.random.default_rng(0)
    psi = rng.normal(size=2 ** lat.n) + 1j * rng.normal(size=2 ** lat.n)
    ref = hm.build_sparse(lat, 1.3, abelian) @ psi
    assert np.abs(hm._apply_numpy(lat, 1.3, psi, abelian) - ref).max() < 1e-11
    assert np.abs(hm.apply_h(lat, 1.3, psi, abelian) - ref).max() < 1e-11
