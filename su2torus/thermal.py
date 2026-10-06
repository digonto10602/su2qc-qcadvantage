"""Thermal (late-time) anchors (docs/METHODS.md, section 5)."""
from __future__ import annotations


def thermal_exact(H, beta: float, ops: dict) -> dict:
    """Canonical averages Tr(e^{-beta H} O)/Z by full diagonalization (n <= 14).
    ops maps name -> sparse/dense matrix. Returns {'energy': ..., name: ...}."""
    raise NotImplementedError


def beta_for_energy_exact(H, energy: float) -> float:
    """beta (may be negative) with <H>_beta = energy, by full diagonalization."""
    raise NotImplementedError


def thermal_typicality(apply_h, dim: int, beta: float, ops: dict, n_samples: int = 10,
                       seed: int = 0) -> dict:
    """Quantum typicality: random vectors r, |r_b> = exp(-beta H/2)|r>,
    <O> ~ sum_r <r_b|O|r_b> / sum_r <r_b|r_b>.  Matrix-free: apply_h(psi) -> H psi;
    ops maps name -> callable psi -> O psi. Returns {'energy': ..., name: ...}."""
    raise NotImplementedError
