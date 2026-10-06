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


def flip_amplitude(g: float, m: int, abelian: bool = False) -> float:
    """Off-diagonal element <s'|H|s> for a flip of a plaquette with m excited neighbours."""
    raise NotImplementedError


def build_sparse(lat, g: float, abelian: bool = False):
    """Full H as scipy.sparse.csr_matrix of shape (2**n, 2**n), float64."""
    raise NotImplementedError


def split_sparse(lat, g: float, abelian: bool = False):
    """(H_E, H_B^A, H_B^B) as csr matrices; their sum equals build_sparse."""
    raise NotImplementedError


def build_pauli(lat, g: float, abelian: bool = False):
    """Full H as qiskit.quantum_info.SparsePauliOp (simplified)."""
    raise NotImplementedError


def apply_h(lat, g: float, psi, abelian: bool = False):
    """Matrix-free H|psi> (numpy, no 2^n x 2^n matrix); needed for n = 24."""
    raise NotImplementedError
