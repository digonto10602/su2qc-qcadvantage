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


def test_mps_observables_match_dense_with_truncation():
    """Grouped environment sweep == dense expectation values of the same (truncated) MPS."""
    from su2torus.lattice import initial_state
    from su2torus import mps
    lat = HoneycombTorus(3, 2)
    for chi in (6, 64):
        res = mps.simulate_mps(lat, 1.25, 0.4, 2, initial=initial_state(lat, "stripe"), chi=chi)
        psi = mps.mps_statevector(res)
        assert np.allclose(res["z"][-1], ob.z_expectations(lat, psi), atol=1e-12)
        assert abs(res["electric_energy"][-1] - ob.electric_energy(lat, 1.25, psi)) < 1e-11
