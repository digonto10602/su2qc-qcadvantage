"""M1 / G1: lattice geometry (LOCKED)."""
import pytest
from su2torus.lattice import HoneycombTorus, initial_state
from tests.acceptance._ref import ref_neighbors, ref_bonds

SIZES = [(2, 2), (3, 2), (2, 3), (3, 3), (4, 3), (3, 4)]


@pytest.mark.parametrize("L1,L2", SIZES)
def test_torus_counts_and_neighbors(L1, L2):
    lat = HoneycombTorus(L1, L2, periodic=True)
    assert lat.n == 2 * L1 * L2
    assert all(len(nb) == 3 for nb in lat.neighbors)
    assert [list(x) for x in lat.neighbors] == ref_neighbors(L1, L2, True)
    assert [tuple(b) for b in lat.bonds] == ref_bonds(ref_neighbors(L1, L2, True))
    assert len(lat.bonds) == 3 * lat.n // 2


@pytest.mark.parametrize("L1,L2", [(3, 3), (3, 4), (4, 4)])
def test_open_neighbors(L1, L2):
    lat = HoneycombTorus(L1, L2, periodic=False)
    assert [list(x) for x in lat.neighbors] == ref_neighbors(L1, L2, False)


@pytest.mark.parametrize("L1,L2", SIZES)
def test_bipartite_sublattice(L1, L2):
    lat = HoneycombTorus(L1, L2)
    assert lat.sublattice == [p % 2 for p in range(lat.n)]
    for p, q in lat.bonds:
        assert lat.sublattice[p] != lat.sublattice[q]


def test_index_convention():
    lat = HoneycombTorus(3, 4)
    assert lat.index(0, 1, 2) == 2 * (1 * 4 + 2)
    assert lat.index(1, 1, 2) == 2 * (1 * 4 + 2) + 1
    assert lat.index(1, -1, 5) == lat.index(1, 2, 1)


@pytest.mark.parametrize("L1,L2", [(3, 3), (3, 4), (4, 4)])
def test_hexagons(L1, L2):
    lat = HoneycombTorus(L1, L2)
    hexes = lat.hexagons()
    assert len(hexes) == L1 * L2
    nb = [set(x) for x in lat.neighbors]
    for h in hexes:
        assert len(h) == 6 and len(set(h)) == 6
        for a, b in zip(h, h[1:] + h[:1]):
            assert b in nb[a], "hexagon not in cyclic order"
    i, j = 1, 2
    expect = {lat.index(0, i, j), lat.index(0, i - 1, j), lat.index(0, i, j - 1),
              lat.index(1, i - 1, j), lat.index(1, i, j - 1), lat.index(1, i - 1, j - 1)}
    assert expect in [set(h) for h in hexes]


def test_initial_states():
    lat = HoneycombTorus(4, 3)
    assert initial_state(lat, "vacuum") == []
    assert initial_state(lat, "neel_half") == sorted(lat.index(0, i, j) for i in range(2) for j in range(3))
    assert initial_state(lat, "stripe") == sorted(lat.index(0, i, j) for i in (0, 2) for j in range(3))
    o = HoneycombTorus(3, 4, periodic=False)
    assert initial_state(o, "string") == sorted(o.index(k, i, 2) for i in range(3) for k in (0, 1))
