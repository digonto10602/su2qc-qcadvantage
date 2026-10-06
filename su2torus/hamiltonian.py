"""Minimal (one qubit per plaquette) SU(2) Hamiltonian, arXiv:2608.28752.

Basis convention (docs/METHODS.md, section 2): qubit p <-> plaquette p; |1> = excited
(n_p = 1, Z_p = -1). Basis index b = sum_p n_p 2**p (Qiskit little-endian).

    H_E = (3 g^2 / 8) sum_p ( 3 n_p - sum_{q in nb(p)} n_p n_q )
        = (3 g^2 / 16) sum_{bonds} (1 - Z_p Z_q)
    H_B = sum_p a_p X_p,   a_p = -(1/g^2) prod_{q in nb(p)} c(n_q),  c(0)=1, c(1)=1/2
    Abelian twin (abelian=True): c == 1, i.e. a_p = -(1/g^2).
    H_B^A / H_B^B : the X terms on sublattice A (up) / B (down) only.
"""
from __future__ import annotations

import itertools

import numpy as np
import scipy.sparse as sp

try:  # optional JIT for the n = 24 matrix-free product
    import numba
except ImportError:  # pragma: no cover
    numba = None


def flip_amplitude(g: float, m: int, abelian: bool = False) -> float:
    """Off-diagonal element <s'|H|s> for a flip of a plaquette with m excited neighbours."""
    return -1.0 / g ** 2 if abelian else -(0.5 ** m) / g ** 2


def electric_diagonal(lat, g: float) -> np.ndarray:
    """Diagonal of H_E over all 2**n basis states (float64)."""
    b = np.arange(2 ** lat.n, dtype=np.uint32 if lat.n < 32 else np.uint64)
    cut = np.zeros(b.shape, dtype=np.uint8)
    for p, q in lat.bonds:
        cut += (((b >> p) ^ (b >> q)) & 1).astype(np.uint8)
    return (3.0 * g ** 2 / 8.0) * cut.astype(np.float64)


def _flip_amplitudes(lat, g: float, p: int, b: np.ndarray, abelian: bool) -> np.ndarray:
    if abelian:
        return np.full(b.shape, -1.0 / g ** 2)
    m = np.zeros(b.shape, dtype=np.int64)
    for q in lat.neighbors[p]:
        m += (b >> q) & 1
    return -(0.5 ** m) / g ** 2


def split_sparse(lat, g: float, abelian: bool = False):
    """(H_E, H_B^A, H_B^B) as csr matrices; their sum equals build_sparse."""
    if lat.n > 20:
        raise MemoryError("n > 20: sparse H would need > 10 GB; use apply_h (matrix-free)")
    dim = 2 ** lat.n
    b = np.arange(dim, dtype=np.int64)
    HE = sp.diags(electric_diagonal(lat, g)).tocsr()
    parts = []
    for s in (0, 1):
        rows, cols, vals = [], [], []
        for p in range(lat.n):
            if lat.sublattice[p] != s:
                continue
            rows.append(b ^ (1 << p)); cols.append(b)
            vals.append(_flip_amplitudes(lat, g, p, b, abelian))
        if rows:
            M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(dim, dim))
        else:
            M = sp.csr_matrix((dim, dim))
        parts.append(M)
    return HE, parts[0], parts[1]


def build_sparse(lat, g: float, abelian: bool = False):
    """Full H as scipy.sparse.csr_matrix of shape (2**n, 2**n), float64."""
    HE, HA, HB = split_sparse(lat, g, abelian)
    return (HE + HA + HB).tocsr()


def flip_pauli_terms(lat, g: float, p: int, abelian: bool = False):
    """[(S, w_S)]: a_p = sum_S w_S prod_{q in S} Z_q (METHODS section 2)."""
    if abelian:
        return [((), -1.0 / g ** 2)]
    nb = lat.neighbors[p]
    k = len(nb)
    out = []
    for r in range(k + 1):
        for S in itertools.combinations(nb, r):
            out.append((S, -(1.0 / g ** 2) * 3.0 ** (k - r) / 4.0 ** k))
    return out


def build_pauli(lat, g: float, abelian: bool = False):
    """Full H as qiskit.quantum_info.SparsePauliOp (simplified)."""
    from qiskit.quantum_info import SparsePauliOp
    n = lat.n
    c = 3.0 * g ** 2 / 16.0
    terms = [("", [], c * len(lat.bonds))]
    terms += [("ZZ", [p, q], -c) for p, q in lat.bonds]
    for p in range(n):
        for S, w in flip_pauli_terms(lat, g, p, abelian):
            terms.append(("X" + "Z" * len(S), [p, *S], w))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n).simplify(atol=1e-14)


# ---------------------------------------------------------------- matrix-free H|psi>
_DIAG_CACHE: dict = {}


def _lat_key(lat, g):
    return (lat.L1, lat.L2, lat.periodic, float(g))


def _diag_cached(lat, g):
    key = _lat_key(lat, g)
    if key not in _DIAG_CACHE:
        if len(_DIAG_CACHE) > 2:
            _DIAG_CACHE.clear()
        _DIAG_CACHE[key] = electric_diagonal(lat, g)
    return _DIAG_CACHE[key]


def _nb_array(lat):
    k = max(len(x) for x in lat.neighbors)
    nb = -np.ones((lat.n, k), dtype=np.int64)
    for p, x in enumerate(lat.neighbors):
        nb[p, :len(x)] = x
    return nb


if numba is not None:
    @numba.njit(parallel=True, cache=True)
    def _apply_numba(psi, diag, nb, inv_g2, abelian, out):
        n = nb.shape[0]
        k = nb.shape[1]
        tab = np.empty(k + 1)
        for m in range(k + 1):
            tab[m] = -inv_g2 * 0.5 ** m
        for b in numba.prange(psi.shape[0]):
            acc = diag[b] * psi[b]
            for p in range(n):
                if abelian:
                    amp = -inv_g2
                else:
                    m = 0
                    for t in range(k):
                        q = nb[p, t]
                        if q >= 0:
                            m += (b >> q) & 1
                    amp = tab[m]
                acc += amp * psi[b ^ (1 << p)]
            out[b] = acc
        return out


def _apply_numpy(lat, g, psi, abelian):
    n = lat.n
    diag = _diag_cached(lat, g)
    out = diag * psi
    t = psi.reshape((2,) * n)
    o = out.reshape((2,) * n)
    for p in range(n):
        ax = n - 1 - p
        shape = [1] * n
        nb = lat.neighbors[p]
        for q in nb:
            shape[n - 1 - q] = 2
        if abelian:
            amp = np.full([1] * n, -1.0 / g ** 2)
        else:
            m = np.zeros([2] * len(nb))
            for a, _ in enumerate(nb):
                sh = [1] * len(nb); sh[a] = 2
                m = m + np.arange(2).reshape(sh)
            # order of broadcast axes follows increasing axis number = decreasing q
            order = np.argsort([n - 1 - q for q in nb])
            amp = (-(0.5 ** m) / g ** 2).transpose(order).reshape(shape)
        o += np.flip(t, axis=ax) * amp
    return out


def apply_h(lat, g: float, psi, abelian: bool = False):
    """Matrix-free H|psi> (numpy, no 2^n x 2^n matrix); needed for n = 24."""
    psi = np.ascontiguousarray(psi)
    if numba is not None:
        dt = np.result_type(psi.dtype, np.float64)
        psi = psi.astype(dt, copy=False)
        out = np.empty_like(psi)
        return _apply_numba(psi, _diag_cached(lat, g), _nb_array(lat), 1.0 / g ** 2,
                            bool(abelian), out)
    return _apply_numpy(lat, g, psi, abelian)
