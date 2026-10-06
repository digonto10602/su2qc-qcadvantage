"""Plaquette graph of the triangular lattice (= honeycomb graph).

Conventions (see docs/METHODS.md, section 1):
- Triangular-lattice sites (i, j), 0 <= i < L1, 0 <= j < L2.
- Each cell (i, j) owns two triangular plaquettes:
    kind 0 = "up"   triangle with vertices (i,j), (i+1,j), (i,j+1)      -> sublattice A
    kind 1 = "down" triangle with vertices (i+1,j), (i,j+1), (i+1,j+1)  -> sublattice B
- Plaquette index p = 2*(i*L2 + j) + kind, so n = 2*L1*L2 plaquettes.
- up(i,j) shares a link with down(i,j), down(i-1,j), down(i,j-1).
  Periodic: indices wrap mod L1, L2. Open: neighbours outside the patch are dropped.
- A hexagon is the set of 6 plaquettes touching one triangular-lattice site (i,j):
  up(i,j), up(i-1,j), up(i,j-1), down(i-1,j), down(i,j-1), down(i-1,j-1),
  returned in cyclic order (consecutive entries share a link, last-first too).
"""
from __future__ import annotations


class HoneycombTorus:
    """Plaquette graph on an L1 x L2 patch of the triangular lattice.

    Attributes (all must exist after __init__):
        L1, L2, periodic
        n          : int, number of plaquettes (= 2*L1*L2)
        neighbors  : list[list[int]], sorted link-sharing neighbours of each plaquette
        bonds      : list[tuple[int,int]], each pair (p, q) with p < q exactly once, sorted
        sublattice : list[int], 0 for up (A), 1 for down (B)
    """

    def __init__(self, L1: int, L2: int, periodic: bool = True):
        raise NotImplementedError

    def index(self, kind: int, i: int, j: int) -> int:
        """Plaquette index; i, j are taken mod L1, L2 when periodic."""
        raise NotImplementedError

    def hexagons(self) -> list[list[int]]:
        """All complete hexagons, each in cyclic order. Periodic: exactly L1*L2 of them."""
        raise NotImplementedError


def initial_state(lat: HoneycombTorus, name: str) -> list[int]:
    """Excited plaquettes of the named initial product state (docs/METHODS.md, section 4):
    'vacuum'   : []
    'neel_half': up plaquettes index(0,i,j) with i <  L1//2
    'stripe'   : up plaquettes index(0,i,j) with i % 2 == 0
    'string'   : index(k, i, L2//2) for all i and k in {0, 1}   (used on open lattices)
    Returned sorted.
    """
    raise NotImplementedError
