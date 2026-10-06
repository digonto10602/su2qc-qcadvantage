"""M3: observables and thermal anchors (LOCKED)."""
import numpy as np
import pytest
from su2torus.lattice import HoneycombTorus
from su2torus import observables as ob, thermal, hamiltonian as hm
from tests.acceptance._ref import ref_hamiltonian, z_op, x_string_op


def rand_state(n, seed=3):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    return v / np.linalg.norm(v)


def ev(op, psi):
    return np.vdot(psi, op @ psi).real


def test_z_link_zz():
    lat = HoneycombTorus(2, 2)
    psi = rand_state(lat.n)
    Z = [z_op(lat.n, p) for p in range(lat.n)]
    assert np.allclose(ob.z_expectations(lat, psi), [ev(z, psi) for z in Z], atol=1e-12)
    ref_links = [0.375 * (1 - ev(Z[p] @ Z[q], psi)) for p, q in lat.bonds]
    assert np.allclose(ob.link_energies(lat, psi), ref_links, atol=1e-12)
    p, q = 0, 5
    ref_c = ev(Z[p] @ Z[q], psi) - ev(Z[p], psi) * ev(Z[q], psi)
    assert abs(ob.zz_connected(lat, psi, p, q) - ref_c) < 1e-12


@pytest.mark.parametrize("abelian", [False, True])
def test_energies(abelian):
    lat = HoneycombTorus(2, 2)
    g = 1.25
    psi = rand_state(lat.n, 5)
    HE, HA, HB = ref_hamiltonian(lat.neighbors, g, abelian, parts=True)
    assert abs(ob.electric_energy(lat, g, psi) - ev(HE, psi)) < 1e-11
    assert abs(ob.magnetic_energy(lat, g, psi, abelian) - ev(HA + HB, psi)) < 1e-11


def test_hexagon_string():
    lat = HoneycombTorus(3, 3)
    psi = rand_state(lat.n, 7)
    h = lat.hexagons()[4]
    assert abs(ob.hexagon_x_string(lat, psi, h) - ev(x_string_op(lat.n, h), psi)) < 1e-12


def test_infinite_temperature():
    lat = HoneycombTorus(3, 2)
    g = 1.1
    d = ob.infinite_temperature(lat, g)
    H = ref_hamiltonian(lat.neighbors, g)
    HE, _, _ = ref_hamiltonian(lat.neighbors, g, parts=True)
    dim = 2 ** lat.n
    assert abs(d["link_energy"] - 0.375) < 1e-14
    assert abs(d["electric_energy"] - HE.diagonal().sum() / dim) < 1e-12
    assert abs(d["magnetic_energy"] - (H - HE).diagonal().sum() / dim) < 1e-14
    assert d["zz_connected"] == 0 and d["hexagon_x_string"] == 0


def test_thermal_exact_and_beta():
    lat = HoneycombTorus(2, 2)
    g = 1.25
    H = ref_hamiltonian(lat.neighbors, g)
    w, V = np.linalg.eigh(H.toarray())
    beta = 0.7
    rho = (V * np.exp(-beta * (w - w[0]))) @ V.conj().T
    rho /= np.trace(rho)
    Z0 = z_op(lat.n, 0).toarray()
    res = thermal.thermal_exact(H, beta, {"z0": Z0})
    assert abs(res["energy"] - np.trace(rho @ H.toarray()).real) < 1e-10
    assert abs(res["z0"] - np.trace(rho @ Z0).real) < 1e-10
    b = thermal.beta_for_energy_exact(H, res["energy"])
    assert abs(b - beta) < 1e-6
    bneg = thermal.beta_for_energy_exact(H, 0.5 * (w.mean() + w[-1]))
    assert bneg < 0


def test_typicality_matches_exact():
    lat = HoneycombTorus(3, 2)            # n = 12, dim = 4096
    g = 1.25
    H = ref_hamiltonian(lat.neighbors, g)
    beta = 0.5
    Z0 = z_op(lat.n, 0)
    ex = thermal.thermal_exact(H, beta, {"z0": Z0.toarray()})
    ty = thermal.thermal_typicality(lambda v: hm.apply_h(lat, g, v), 2 ** lat.n, beta,
                                    {"z0": lambda v: Z0 @ v}, n_samples=12, seed=11)
    assert abs(ty["energy"] - ex["energy"]) < 0.02 * abs(ex["energy"]) + 0.05
    assert abs(ty["z0"] - ex["z0"]) < 0.05
