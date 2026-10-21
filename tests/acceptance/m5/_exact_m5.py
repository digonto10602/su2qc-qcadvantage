"""Independent exact Trotter-circuit statevectors for the M5 acceptance tests (no su2torus code).

Built only from tests/acceptance/_ref.py (brute-force plaquette graph) and the formulas of
docs/METHODS.md sections 2-3. The state is held as a (2,)*n tensor whose axis n-1-p is qubit p, so
the C-order flattening has basis index b = sum_p n_p 2**p (METHODS section 2).

Trotter steps (METHODS section 3), U_X(t) = exp(-i H_X t):
    order 1:  U_B(dt) U_A(dt) U_E(dt)                                  (U_E acts first)
    order 2:  U_E(dt/2) U_A(dt/2) U_B(dt) U_A(dt/2) U_E(dt/2)
U_E is applied as the exact diagonal phase exp(-i t H_E) (differs from the rzz circuit only by a
global phase). Each flip exp(-i t a_p X_p) is applied exactly in every sector of the neighbour
occupations, which the flip does not change (neighbours lie on the other sublattice).
"""
import numpy as np

from tests.acceptance._ref import ref_bonds, ref_neighbors


def stripe(L1, L2):
    """METHODS section 4 'stripe': up(i,j) = 2(i L2 + j) for even i."""
    return sorted(2 * (i * L2 + j) for i in range(L1) for j in range(L2) if i % 2 == 0)


def exact_trotter(L1, L2, g, dt, steps, excited, order=1, abelian=False, observe=None):
    """Statevector (2**n,) after `steps` Trotter steps from the product state `excited` on the
    periodic L1 x L2 torus. If `observe` is given, return [observe(state) for t = 0..steps] instead
    (memory-light: only one 2**n state is kept)."""
    nb = ref_neighbors(L1, L2, True)
    n = len(nb)
    bonds = ref_bonds(nb)

    def ax(p):
        return n - 1 - p

    def zfield(p):
        shape = [1] * n
        shape[ax(p)] = 2
        return np.array([1.0, -1.0]).reshape(shape)

    he = np.zeros((2,) * n)
    for p, q in bonds:
        he += (3.0 * g ** 2 / 16.0) * (1.0 - zfield(p) * zfield(q))
    exc = set(excited)
    psi = np.zeros((2,) * n, dtype=complex)
    psi[tuple(1 if (n - 1 - a) in exc else 0 for a in range(n))] = 1.0

    def u_e(t):
        psi[...] *= np.exp(-1j * t * he)

    def flip(p, t):
        k = len(nb[p])
        for cfg in range(2 ** k):
            m = bin(cfg).count("1")
            a = -(1.0 / g ** 2) * (1.0 if abelian else 0.5 ** m)
            c, s = np.cos(t * a), np.sin(t * a)
            idx = [slice(None)] * n
            for r, q in enumerate(nb[p]):
                idx[ax(q)] = (cfg >> r) & 1
            i0 = list(idx)
            i0[ax(p)] = 0
            i1 = list(idx)
            i1[ax(p)] = 1
            v0 = psi[tuple(i0)].copy()
            v1 = psi[tuple(i1)].copy()
            psi[tuple(i0)] = c * v0 - 1j * s * v1
            psi[tuple(i1)] = c * v1 - 1j * s * v0

    def u_sub(sub, t):
        for p in range(n):
            if p % 2 == sub:
                flip(p, t)

    out = [observe(psi.reshape(-1))] if observe else None
    for _ in range(steps):
        if order == 1:
            u_e(dt)
            u_sub(0, dt)
            u_sub(1, dt)
        elif order == 2:
            u_e(dt / 2)
            u_sub(0, dt / 2)
            u_sub(1, dt)
            u_sub(0, dt / 2)
            u_e(dt / 2)
        else:
            raise ValueError(order)
        if observe:
            out.append(observe(psi.reshape(-1)))
    return out if observe else psi.reshape(-1)


def z_and_he(psi, L1, L2, g):
    """(<Z_p> for all p, <H_E>) of a normalised statevector, H_E = (3g^2/16) sum_bonds (1 - Z_p Z_q)."""
    nb = ref_neighbors(L1, L2, True)
    n = len(nb)
    prob = (np.abs(psi) ** 2).reshape((2,) * n)
    z = []
    for p in range(n):
        m = prob.sum(axis=tuple(a for a in range(n) if a != n - 1 - p))
        z.append(m[0] - m[1])
    zz = 0.0
    for p, q in ref_bonds(nb):
        m = prob.sum(axis=tuple(a for a in range(n) if a not in (n - 1 - p, n - 1 - q)))
        zz_pq = m[0, 0] + m[1, 1] - m[0, 1] - m[1, 0]
        zz += (3.0 * g ** 2 / 16.0) * (1.0 - zz_pq)
    return np.array(z), float(zz)
