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
        self.L1, self.L2, self.periodic = int(L1), int(L2), bool(periodic)
        if self.periodic and min(self.L1, self.L2) < 2:
            raise ValueError("torus needs L1, L2 >= 2 (L = 1 merges distinct links)")
        self.n = 2 * self.L1 * self.L2
        nb = [set() for _ in range(self.n)]
        for i in range(self.L1):
            for j in range(self.L2):
                u = self.index(0, i, j)
                for (a, b) in ((i, j), (i - 1, j), (i, j - 1)):
                    if self._inside(a, b):
                        d = self.index(1, a, b)
                        nb[u].add(d)
                        nb[d].add(u)
        self.neighbors = [sorted(s) for s in nb]
        self.bonds = sorted({(min(p, q), max(p, q)) for p in range(self.n) for q in self.neighbors[p]})
        self.sublattice = [p % 2 for p in range(self.n)]

    def _inside(self, i: int, j: int) -> bool:
        return self.periodic or (0 <= i < self.L1 and 0 <= j < self.L2)

    def index(self, kind: int, i: int, j: int) -> int:
        """Plaquette index; i, j are taken mod L1, L2 when periodic."""
        if self.periodic:
            i, j = i % self.L1, j % self.L2
        elif not (0 <= i < self.L1 and 0 <= j < self.L2):
            raise IndexError((kind, i, j))
        return 2 * (i * self.L2 + j) + kind

    def hexagons(self) -> list[list[int]]:
        """All complete hexagons, each in cyclic order. Periodic: exactly L1*L2 of them."""
        out = []
        for i in range(self.L1):
            for j in range(self.L2):
                cells = [(0, i, j), (1, i, j - 1), (0, i, j - 1),
                         (1, i - 1, j - 1), (0, i - 1, j), (1, i - 1, j)]
                if all(self._inside(a, b) for _, a, b in cells):
                    out.append([self.index(k, a, b) for k, a, b in cells])
        return out


def initial_state(lat: HoneycombTorus, name: str) -> list[int]:
    """Excited plaquettes of the named initial product state (docs/METHODS.md, section 4):
    'vacuum'   : []
    'neel_half': up plaquettes index(0,i,j) with i <  L1//2
    'stripe'   : up plaquettes index(0,i,j) with i % 2 == 0
    'string'   : index(k, i, L2//2) for all i and k in {0, 1}   (used on open lattices)
    Returned sorted.
    """
    L1, L2 = lat.L1, lat.L2
    if name == "vacuum":
        ex = []
    elif name == "neel_half":
        ex = [lat.index(0, i, j) for i in range(L1 // 2) for j in range(L2)]
    elif name == "stripe":
        ex = [lat.index(0, i, j) for i in range(0, L1, 2) for j in range(L2)]
    elif name == "string":
        ex = [lat.index(k, i, L2 // 2) for i in range(L1) for k in (0, 1)]
    else:
        raise ValueError(f"unknown initial state {name!r}")
    return sorted(ex)
