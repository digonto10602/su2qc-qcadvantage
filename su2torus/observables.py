"""Observables on state vectors (docs/METHODS.md, section 4). psi is a normalized numpy vector."""
from __future__ import annotations


def z_expectations(lat, psi):
    """Array of <Z_p>, length lat.n."""
    raise NotImplementedError


def link_energies(lat, psi):
    """Array over lat.bonds of <E^2_pq> = (3/8) <1 - Z_p Z_q>."""
    raise NotImplementedError


def zz_connected(lat, psi, p: int, q: int) -> float:
    """<Z_p Z_q> - <Z_p><Z_q>."""
    raise NotImplementedError


def electric_energy(lat, g: float, psi) -> float:
    """<H_E> = (g^2/2) * sum of link energies."""
    raise NotImplementedError


def magnetic_energy(lat, g: float, psi, abelian: bool = False) -> float:
    """<H_B>."""
    raise NotImplementedError


def hexagon_x_string(lat, psi, hexagon) -> float:
    """<prod_{p in hexagon} X_p>."""
    raise NotImplementedError


def infinite_temperature(lat, g: float, abelian: bool = False) -> dict:
    """Fully decohered (maximally mixed) values:
    {'link_energy', 'zz_connected', 'electric_energy', 'magnetic_energy', 'hexagon_x_string'}."""
    raise NotImplementedError
