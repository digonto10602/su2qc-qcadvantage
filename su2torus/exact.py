"""Exact (state-vector) tools."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh, expm_multiply
from scipy.special import jv


def product_state(lat, excited):
    """Complex128 basis vector with the listed plaquettes excited."""
    v = np.zeros(2 ** lat.n, dtype=np.complex128)
    v[sum(1 << int(p) for p in set(excited))] = 1.0
    return v


def ground_state(H):
    """(E0, psi0) of a sparse Hermitian matrix (Lanczos)."""
    w, V = eigsh(H, k=1, which="SA", tol=1e-12)
    return float(w[0]), V[:, 0]


def spectral_bounds(H, pad: float = 0.05):
    """Padded [Emin, Emax] from Lanczos (for Chebyshev propagation)."""
    lo = eigsh(H, k=1, which="SA", tol=1e-6, return_eigenvectors=False)[0]
    hi = eigsh(H, k=1, which="LA", tol=1e-6, return_eigenvectors=False)[0]
    w = hi - lo
    return float(lo - pad * w - 1e-3), float(hi + pad * w + 1e-3)


def chebyshev_step(matvec, psi, t, bounds, tol=1e-15):
    """exp(-i H t) psi by Chebyshev expansion; the spectrum of H must lie inside `bounds`."""
    lo, hi = bounds
    c, r = 0.5 * (hi + lo), 0.5 * (hi - lo)
    x = r * abs(t)
    kmax = int(x + 10 * x ** (1 / 3) + 40)
    J = jv(np.arange(kmax + 1), x)
    s = -1j if t >= 0 else 1j
    Ht = lambda v: (matvec(v) - c * v) / r
    t0 = np.array(psi, dtype=np.complex128)
    out = J[0] * t0
    t1 = Ht(t0)
    out += 2 * s * J[1] * t1
    for k in range(2, kmax + 1):
        t2 = 2 * Ht(t1) - t0
        out += 2 * s ** k * J[k] * t2
        t0, t1 = t1, t2
        if k > x and abs(J[k]) < tol and abs(J[k - 1]) < tol:
            break
    return np.exp(-1j * c * t) * out


def evolve_iter(H, psi0, times, bounds=None):
    """Generator of exp(-i H t) psi0 for t in `times` (one state alive at a time)."""
    times = np.asarray(times, dtype=float)
    psi = np.asarray(psi0, dtype=np.complex128)
    if isinstance(H, LinearOperator):
        if bounds is None:
            bounds = spectral_bounds(H)
        step = lambda v, dt: chebyshev_step(H.matvec, v, dt, bounds)
    else:
        A = -1j * (sp.csr_matrix(H) if sp.issparse(H) else np.asarray(H))
        step = lambda v, dt: expm_multiply(A * dt, v)
    tprev = 0.0
    for t in times:
        if t != tprev:
            psi = step(psi, t - tprev)
        yield psi
        tprev = t


def evolve(H, psi0, times, bounds=None):
    """States exp(-i H t) psi0 for each t in `times` (ascending, may start at 0).
    H may be a sparse matrix or a scipy LinearOperator. Returns array (len(times), dim).
    Sparse/dense H: scipy expm_multiply. LinearOperator: Chebyshev propagation
    (spectral bounds from Lanczos unless given)."""
    return np.array(list(evolve_iter(H, psi0, times, bounds)), dtype=np.complex128)
