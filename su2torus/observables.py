"""Observables on state vectors (docs/METHODS.md, section 4). psi is a normalized numpy vector."""
from __future__ import annotations

import numpy as np

from . import hamiltonian as hm


def _probs(psi):
    return np.abs(psi) ** 2


def _bit_sign(n, p, b=None):
    """(1 - 2 n_p) over the basis, as int8."""
    if b is None:
        b = np.arange(2 ** n, dtype=np.uint32)
    return (1 - 2 * ((b >> p) & 1)).astype(np.int8)


def z_expectations(lat, psi):
    """Array of <Z_p>, length lat.n."""
    n = lat.n
    pr = _probs(psi)
    out = np.empty(n)
    for p in range(n):
        m = pr.reshape(2 ** (n - 1 - p), 2, 2 ** p).sum(axis=(0, 2))
        out[p] = m[0] - m[1]
    return out


def zz_expectations(lat, psi, pairs):
    """Array of <Z_p Z_q> for the listed pairs."""
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
    idx = np.arange(psi.size, dtype=np.uint32 if lat.n < 32 else np.uint64)
    return float(np.vdot(psi, psi[idx ^ mask]).real)


def infinite_temperature(lat, g: float, abelian: bool = False) -> dict:
    """Fully decohered (maximally mixed) values:
    {'link_energy', 'zz_connected', 'electric_energy', 'magnetic_energy', 'hexagon_x_string'}."""
    return {"link_energy": 0.375, "zz_connected": 0.0,
            "electric_energy": 0.5 * g ** 2 * 0.375 * len(lat.bonds),
            "magnetic_energy": 0.0, "hexagon_x_string": 0.0, "z": 0.0}
