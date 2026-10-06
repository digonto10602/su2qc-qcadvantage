"""M1 / G2-G5: Hamiltonian correctness (LOCKED)."""
import numpy as np
import pytest
import scipy.sparse as sp
from su2torus.lattice import HoneycombTorus
from su2torus import hamiltonian as hm
from tests.acceptance._ref import ref_hamiltonian

G_VALUES = [0.9, 1.25, 1.4]


def norm_diff(a, b):
    d = (a - b)
    d = d.toarray() if sp.issparse(d) else np.asarray(d)
    return np.abs(d).max()


@pytest.mark.parametrize("dims", [(2, 2), (3, 2)])
@pytest.mark.parametrize("g", G_VALUES)
@pytest.mark.parametrize("abelian", [False, True])
def test_sparse_matches_bruteforce(dims, g, abelian):          # G2, G4, G5
    lat = HoneycombTorus(*dims)
    H = hm.build_sparse(lat, g, abelian)
    assert H.shape == (2 ** lat.n, 2 ** lat.n)
    assert norm_diff(H, ref_hamiltonian(lat.neighbors, g, abelian)) < 1e-12


@pytest.mark.parametrize("dims", [(2, 2), (3, 2)])
@pytest.mark.parametrize("abelian", [False, True])
def test_pauli_form_matches(dims, abelian):                     # G2
    lat = HoneycombTorus(*dims)
    g = 1.25
    P = hm.build_pauli(lat, g, abelian).to_matrix(sparse=True)
    assert norm_diff(P, ref_hamiltonian(lat.neighbors, g, abelian)) < 1e-12


@pytest.mark.parametrize("abelian", [False, True])
def test_split_parts(abelian):
    lat = HoneycombTorus(2, 2)
    g = 1.1
    HE, HA, HB = hm.split_sparse(lat, g, abelian)
    rE, rA, rB = ref_hamiltonian(lat.neighbors, g, abelian, parts=True)
    for a, b in ((HE, rE), (HA, rA), (HB, rB)):
        assert norm_diff(a, b) < 1e-12


def test_hermitian_and_stoquastic():                            # G3
    lat = HoneycombTorus(3, 2)
    H = hm.build_sparse(lat, 1.25)
    assert norm_diff(H, H.conj().T) < 1e-14
    off = (H - sp.diags(H.diagonal())).tocoo()
    assert off.data.max() <= 1e-15, "off-diagonal elements must be <= 0"


@pytest.mark.parametrize("m", [0, 1, 2, 3])
def test_flip_amplitude(m):                                     # G5
    g = 1.3
    assert abs(hm.flip_amplitude(g, m) - (-(0.5 ** m) / g ** 2)) < 1e-15
    assert abs(hm.flip_amplitude(g, m, abelian=True) - (-1 / g ** 2)) < 1e-15


def test_explicit_matrix_elements():                            # G5
    lat = HoneycombTorus(2, 2)
    g = 1.2
    H = hm.build_sparse(lat, g).tolil()
    p = 0
    nb = lat.neighbors[p]
    for m in range(4):
        b = sum(1 << q for q in nb[:m])
        assert abs(H[b ^ 1, b] - (-(0.5 ** m) / g ** 2)) < 1e-14
    assert abs(H[1, 0] - (-1 / g ** 2)) < 1e-14


@pytest.mark.parametrize("abelian", [False, True])
def test_matrix_free(abelian):
    lat = HoneycombTorus(3, 2)
    g = 1.25
    rng = np.random.default_rng(1)
    psi = rng.normal(size=2 ** lat.n) + 1j * rng.normal(size=2 ** lat.n)
    ref = ref_hamiltonian(lat.neighbors, g, abelian) @ psi
    assert np.abs(hm.apply_h(lat, g, psi, abelian) - ref).max() < 1e-11
