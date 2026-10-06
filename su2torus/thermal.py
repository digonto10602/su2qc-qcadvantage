"""Thermal (late-time) anchors (docs/METHODS.md, section 5)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, expm_multiply

_EIG_CACHE: dict = {}


def _eig(H, vectors=True):
    key = (id(H), H.shape, vectors)
    if key not in _EIG_CACHE:
        _EIG_CACHE.clear()
        Hd = H.toarray() if sp.issparse(H) else np.asarray(H)
        _EIG_CACHE[key] = np.linalg.eigh(Hd) if vectors else (np.linalg.eigvalsh(Hd), None)
        _EIG_CACHE[key] = (_EIG_CACHE[key], H)   # keep H alive so id() stays unique
    return _EIG_CACHE[key][0]


def _weights(w, beta):
    x = -beta * w
    x = x - x.max()
    p = np.exp(x)
    return p / p.sum()


def _energy(w, beta):
    return float(_weights(w, beta) @ w)


def thermal_exact(H, beta: float, ops: dict) -> dict:
    """Canonical averages Tr(e^{-beta H} O)/Z by full diagonalization (n <= 14).
    ops maps name -> sparse/dense matrix. Returns {'energy': ..., name: ...}."""
    w, V = _eig(H, True)
    p = _weights(w, beta)
    out = {"energy": float(p @ w)}
    for name, O in ops.items():
        OV = O @ V
        diag = np.einsum("ij,ij->j", V.conj(), np.asarray(OV)).real
        out[name] = float(p @ diag)
    return out


def beta_for_energy_exact(H, energy: float) -> float:
    """beta (may be negative) with <H>_beta = energy, by full diagonalization."""
    if (id(H), H.shape, True) in _EIG_CACHE:
        w = _eig(H, True)[0]
    else:
        w = _eig(H, False)[0]
    if not (w[0] < energy < w[-1]):
        raise ValueError("energy outside the spectrum")
    f = lambda b: _energy(w, b) - energy
    lo, hi = -1.0, 1.0
    while f(lo) < 0:
        lo *= 2
    while f(hi) > 0:
        hi *= 2
    return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-12))


def _random_vectors(dim, n_samples, seed):
    rng = np.random.default_rng(seed)
    for _ in range(n_samples):
        r = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        yield r / np.linalg.norm(r)


def thermal_typicality(apply_h, dim: int, beta: float, ops: dict, n_samples: int = 10,
                       seed: int = 0) -> dict:
    """Quantum typicality: random vectors r, |r_b> = exp(-beta H/2)|r>,
    <O> ~ sum_r <r_b|O|r_b> / sum_r <r_b|r_b>.  Matrix-free: apply_h(psi) -> H psi;
    ops maps name -> callable psi -> O psi. Returns {'energy': ..., name: ...}."""
    A = LinearOperator((dim, dim), matvec=lambda v: -0.5 * beta * apply_h(v),
                       rmatvec=lambda v: -0.5 * beta * apply_h(v), dtype=complex)
    num = {k: 0.0 for k in ["energy", *ops]}
    den = 0.0
    for r in _random_vectors(dim, n_samples, seed):
        rb = expm_multiply(A, r, traceA=0.0) if beta != 0 else r
        den += np.vdot(rb, rb).real
        num["energy"] += np.vdot(rb, apply_h(rb)).real
        for name, op in ops.items():
            num[name] += np.vdot(rb, op(rb)).real
    return {k: float(v / den) for k, v in num.items()}


def beta_for_energy_typicality(apply_h, dim, energy, ops, beta0=0.0, step=0.1,
                               n_samples=4, seed=0, max_evals=8, rtol=0.02):
    """Bracketed false-position (Illinois) search for beta with <H>_beta = energy, using
    thermal_typicality with fixed random vectors (same seed for every evaluation).
    Stops when |E - energy| <= rtol*max(1,|energy|)/4, or after max_evals evaluations.
    Returns (beta, result_dict, history), the best evaluation found."""
    hist = []

    def ev(b):
        r = thermal_typicality(apply_h, dim, b, ops, n_samples, seed)
        hist.append((b, r))
        return r["energy"] - energy

    tol = 0.25 * rtol * max(1.0, abs(energy))
    a, fa = beta0, ev(beta0)
    if abs(fa) <= tol:
        return a, hist[-1][1], hist
    # E(beta) decreases with beta: walk downhill/uphill to bracket the root
    d = step if fa > 0 else -step
    b, fb = a + d, ev(a + d)
    while fa * fb > 0 and len(hist) < max_evals:
        a, fa = b, fb
        d *= 2
        b, fb = a + d, ev(a + d)
    while len(hist) < max_evals and abs(fb) > tol and fa * fb < 0:
        c = b - fb * (b - a) / (fb - fa)
        fc = ev(c)
        if fc * fb < 0:
            a, fa = b, fb
        else:
            fa *= 0.5            # Illinois modification
        b, fb = c, fc
    best = min(hist, key=lambda h: abs(h[1]["energy"] - energy))
    return best[0], best[1], hist
