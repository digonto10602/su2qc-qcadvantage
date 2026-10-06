"""Independent brute-force references for acceptance tests (LOCKED; do not edit).

Built from triangle geometry and the Hamiltonian formulas in docs/METHODS.md,
without using any su2torus code."""
import itertools
import numpy as np
import scipy.sparse as sp


def triangles(L1, L2):
    """Vertex sets of all plaquettes, indexed p = 2*(i*L2+j)+kind (unwrapped coords)."""
    tri = {}
    for i in range(L1):
        for j in range(L2):
            tri[2 * (i * L2 + j) + 0] = [(i, j), (i + 1, j), (i, j + 1)]
            tri[2 * (i * L2 + j) + 1] = [(i + 1, j), (i, j + 1), (i + 1, j + 1)]
    return tri


def _edge(a, b, L1, L2, periodic):
    """Link label = (type, base vertex); types (1,0), (0,1), (1,-1). Distinguishes the two
    different links that join the same pair of sites on a 2-wide torus."""
    d = (b[0] - a[0], b[1] - a[1])
    if d not in ((1, 0), (0, 1), (1, -1)):
        a, d = b, (-d[0], -d[1])
    base = (a[0] % L1, a[1] % L2) if periodic else a
    return (d, base)


def ref_neighbors(L1, L2, periodic):
    tri = triangles(L1, L2)
    edges = {p: {_edge(a, b, L1, L2, periodic) for a, b in itertools.combinations(vs, 2)}
             for p, vs in tri.items()}
    n = 2 * L1 * L2
    nb = [sorted(q for q in range(n) if q != p and edges[p] & edges[q]) for p in range(n)]
    return nb


def ref_bonds(nb):
    return sorted({(min(p, q), max(p, q)) for p in range(len(nb)) for q in nb[p]})


def occ(b, n):
    return [(b >> p) & 1 for p in range(n)]


def ref_hamiltonian(nb, g, abelian=False, parts=False):
    """Brute force: H_E diagonal and X flips with amplitude -(1/g^2) prod c(n_q)."""
    n = len(nb)
    dim = 2 ** n
    bonds = ref_bonds(nb)
    diag = np.zeros(dim)
    rows, cols, vals, subl = [], [], [], []
    for b in range(dim):
        o = occ(b, n)
        diag[b] = 0.5 * g ** 2 * sum(0.375 * (1 - (1 - 2 * o[p]) * (1 - 2 * o[q])) for p, q in bonds)
        for p in range(n):
            m = sum(o[q] for q in nb[p])
            amp = -(1.0 / g ** 2) * (1.0 if abelian else 0.5 ** m)
            rows.append(b ^ (1 << p)); cols.append(b); vals.append(amp); subl.append(p % 2)
    HE = sp.diags(diag).tocsr()
    rows, cols, vals, subl = map(np.array, (rows, cols, vals, subl))
    HA = sp.csr_matrix((vals[subl == 0], (rows[subl == 0], cols[subl == 0])), shape=(dim, dim))
    HB = sp.csr_matrix((vals[subl == 1], (rows[subl == 1], cols[subl == 1])), shape=(dim, dim))
    if parts:
        return HE, HA, HB
    return (HE + HA + HB).tocsr()


def z_op(n, p):
    d = np.array([1 - 2 * ((b >> p) & 1) for b in range(2 ** n)], dtype=float)
    return sp.diags(d).tocsr()


def x_string_op(n, plaqs):
    mask = sum(1 << p for p in plaqs)
    dim = 2 ** n
    idx = np.arange(dim)
    return sp.csr_matrix((np.ones(dim), (idx ^ mask, idx)), shape=(dim, dim))


def basis_vec(n, excited):
    v = np.zeros(2 ** n, dtype=complex)
    v[sum(1 << p for p in excited)] = 1.0
    return v
