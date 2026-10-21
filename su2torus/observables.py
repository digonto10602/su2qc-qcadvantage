"""Observables on state vectors (docs/METHODS.md, section 4). psi is a normalized numpy vector."""
from __future__ import annotations

import numpy as np

from . import hamiltonian as hm

try:
    import numba
except ImportError:  # pragma: no cover
    numba = None

_FAST_N = 16          # use the parallel kernels above this many qubits

if numba is not None:
    @numba.njit(parallel=True, cache=True)
    def _z_zz_numba(psi, n, P, Q, nchunk):
        """Per-chunk sums of |psi|^2 (1-2n_p) and |psi|^2 (1-2 parity(p,q))."""
        dim = psi.shape[0]
        zs = np.zeros((nchunk, n))
        zzs = np.zeros((nchunk, P.shape[0]))
        size = (dim + nchunk - 1) // nchunk
        for c in numba.prange(nchunk):
            for b in range(c * size, min(dim, (c + 1) * size)):
                w = psi[b].real ** 2 + psi[b].imag ** 2
                for p in range(n):
                    zs[c, p] += w * (1 - 2 * ((b >> p) & 1))
                for k in range(P.shape[0]):
                    zzs[c, k] += w * (1 - 2 * (((b >> P[k]) ^ (b >> Q[k])) & 1))
        return zs.sum(axis=0), zzs.sum(axis=0)

    @numba.njit(parallel=True, cache=True)
    def _flip_overlap(psi, mask, nchunk):
        dim = psi.shape[0]
        acc = np.zeros(nchunk)
        size = (dim + nchunk - 1) // nchunk
        for c in numba.prange(nchunk):
            for b in range(c * size, min(dim, (c + 1) * size)):
                v = np.conj(psi[b]) * psi[b ^ mask]
                acc[c] += v.real
        return acc.sum()


def _fast(lat, psi):
    return numba is not None and lat.n > _FAST_N


def _z_zz(lat, psi, pairs):
    P = np.array([p for p, _ in pairs] or [0], dtype=np.int64)[:len(pairs) or 0]
    Q = np.array([q for _, q in pairs] or [0], dtype=np.int64)[:len(pairs) or 0]
    return _z_zz_numba(np.ascontiguousarray(psi, dtype=np.complex128), lat.n, P, Q, 256)


def _probs(psi):
    return np.abs(psi) ** 2


def _bit_sign(n, p, b=None):
    """(1 - 2 n_p) over the basis, as int8."""
    if b is None:
        b = np.arange(2 ** n, dtype=np.uint32)
    return (1 - 2 * ((b >> p) & 1)).astype(np.int8)


def z_expectations(lat, psi):
    """Array of <Z_p>, length lat.n."""
    if _fast(lat, psi):
        return _z_zz(lat, psi, [])[0]
    n = lat.n
    pr = _probs(psi)
    out = np.empty(n)
    for p in range(n):
        m = pr.reshape(2 ** (n - 1 - p), 2, 2 ** p).sum(axis=(0, 2))
        out[p] = m[0] - m[1]
    return out


def zz_expectations(lat, psi, pairs):
    """Array of <Z_p Z_q> for the listed pairs."""
    if _fast(lat, psi):
        return _z_zz(lat, psi, list(pairs))[1]
    pr = _probs(psi)
    b = np.arange(pr.size, dtype=np.uint32)
    out = np.empty(len(pairs))
    for k, (p, q) in enumerate(pairs):
        par = ((b >> p) ^ (b >> q)) & 1
        out[k] = pr.sum() - 2.0 * pr[par == 1].sum()
    return out


def link_energies(lat, psi):
    """Array over lat.bonds of <E^2_pq> = (3/8) <1 - Z_p Z_q>."""
    return 0.375 * (1.0 - zz_expectations(lat, psi, lat.bonds))


def zz_connected(lat, psi, p: int, q: int) -> float:
    """<Z_p Z_q> - <Z_p><Z_q>."""
    z = z_expectations(lat, psi)
    return float(zz_expectations(lat, psi, [(p, q)])[0] - z[p] * z[q])


def electric_energy(lat, g: float, psi) -> float:
    """<H_E> = (g^2/2) * sum of link energies."""
    return float(0.5 * g ** 2 * link_energies(lat, psi).sum())


def magnetic_energy(lat, g: float, psi, abelian: bool = False) -> float:
    """<H_B>."""
    Hpsi = hm.apply_h(lat, g, psi, abelian)
    tot = np.vdot(psi, Hpsi).real
    return float(tot - electric_energy(lat, g, psi))


def hexagon_x_string(lat, psi, hexagon) -> float:
    """<prod_{p in hexagon} X_p>."""
    mask = sum(1 << p for p in hexagon)
    if _fast(lat, psi):
        return float(_flip_overlap(np.ascontiguousarray(psi, dtype=np.complex128), mask, 256))
    idx = np.arange(psi.size, dtype=np.uint32 if lat.n < 32 else np.uint64)
    return float(np.vdot(psi, psi[idx ^ mask]).real)


def infinite_temperature(lat, g: float, abelian: bool = False) -> dict:
    """Fully decohered (maximally mixed) values:
    {'link_energy', 'zz_connected', 'electric_energy', 'magnetic_energy', 'hexagon_x_string'}."""
    return {"link_energy": 0.375, "zz_connected": 0.0,
            "electric_energy": 0.5 * g ** 2 * 0.375 * len(lat.bonds),
            "magnetic_energy": 0.0, "hexagon_x_string": 0.0, "z": 0.0}
