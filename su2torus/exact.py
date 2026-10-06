"""Exact (state-vector) tools."""
from __future__ import annotations


def product_state(lat, excited):
    """Complex128 basis vector with the listed plaquettes excited."""
    raise NotImplementedError


def ground_state(H):
    """(E0, psi0) of a sparse Hermitian matrix (Lanczos)."""
    raise NotImplementedError


def evolve(H, psi0, times):
    """States exp(-i H t) psi0 for each t in `times` (ascending, may start at 0).
    H may be a sparse matrix or a scipy LinearOperator. Returns array (len(times), dim)."""
    raise NotImplementedError
