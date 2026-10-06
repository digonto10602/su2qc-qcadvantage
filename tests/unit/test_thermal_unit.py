from su2torus.lattice import HoneycombTorus
from su2torus import hamiltonian as hm, thermal


def test_beta_search_typicality_matches_exact():
    lat = HoneycombTorus(3, 2)
    g = 1.25
    H = hm.build_sparse(lat, g)
    E = thermal.thermal_exact(H, 0.4, {})["energy"]
    b, res, hist = thermal.beta_for_energy_typicality(
        lambda v: hm.apply_h(lat, g, v), 2 ** lat.n, E, {}, beta0=0.0, n_samples=8, seed=2)
    assert len(hist) <= 8
    assert abs(res["energy"] - E) < 0.005 * max(1, abs(E))
    assert abs(b - 0.4) < 0.1


def test_chebyshev_typicality_equals_expm_multiply():
    lat = HoneycombTorus(3, 2)
    g = 1.25
    H = hm.build_sparse(lat, g)
    from su2torus.exact import spectral_bounds
    bd = spectral_bounds(H)
    mv = lambda v: hm.apply_h(lat, g, v)
    for beta in (0.6, -0.3):
        a = thermal.thermal_typicality(mv, 2 ** lat.n, beta, {}, 3, 5)
        b = thermal.thermal_typicality(mv, 2 ** lat.n, beta, {}, 3, 5, bounds=bd)
        assert abs(a["energy"] - b["energy"]) < 1e-9
