import numpy as np
from su2torus.lattice import HoneycombTorus
from su2torus import observables as ob


def test_fast_kernels_match_numpy(monkeypatch):
    lat = HoneycombTorus(3, 3, periodic=False)          # n = 18
    rng = np.random.default_rng(4)
    psi = rng.normal(size=2 ** lat.n) + 1j * rng.normal(size=2 ** lat.n)
    psi /= np.linalg.norm(psi)
    fast = (ob.z_expectations(lat, psi), ob.link_energies(lat, psi),
            ob.hexagon_x_string(lat, psi, lat.hexagons()[0]))
    monkeypatch.setattr(ob, "_FAST_N", 99)
    slow = (ob.z_expectations(lat, psi), ob.link_energies(lat, psi),
            ob.hexagon_x_string(lat, psi, lat.hexagons()[0]))
    for a, b in zip(fast, slow):
        assert np.allclose(a, b, atol=1e-12)
