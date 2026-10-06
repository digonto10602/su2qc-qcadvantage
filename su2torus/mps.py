"""MPS simulation of the Trotter circuits with quimb CircuitMPS (PLAN M5; API fixed by
tests/acceptance/m5). The circuit is exactly su2torus.trotter.trotter_circuit; gates are applied
one by one as explicit matrices (taken from Qiskit, so the conventions cannot drift), on an MPS
whose site k holds plaquette snake_order(lat)[k]."""
from __future__ import annotations

import time

import numpy as np

from .trotter import trotter_circuit, trotter_step

GPU_BACKENDS = ("cupy", "torch")
_Z = np.diag([1.0, -1.0])


def snake_order(lat) -> list[int]:
    """Boustrophedon over cells, lines along the shorter side; the two plaquettes of a cell on
    consecutive sites, (up, down) on even lines and (down, up) on odd lines."""
    out = []
    if lat.L1 <= lat.L2:
        for j in range(lat.L2):
            rng = range(lat.L1) if j % 2 == 0 else reversed(range(lat.L1))
            for i in rng:
                pair = [lat.index(0, i, j), lat.index(1, i, j)]
                out += pair if j % 2 == 0 else pair[::-1]
    else:
        for i in range(lat.L1):
            rng = range(lat.L2) if i % 2 == 0 else reversed(range(lat.L2))
            for j in rng:
                pair = [lat.index(0, i, j), lat.index(1, i, j)]
                out += pair if i % 2 == 0 else pair[::-1]
    return out


def _to_backend(backend):
    if backend == "numpy":
        return None
    if backend == "cupy":
        import cupy  # noqa: F401  (ImportError if missing: never fall back to CPU)
        return lambda x: cupy.asarray(x)
    if backend == "torch":
        import torch
        if not torch.cuda.is_available():
            raise RuntimeError("backend='torch' requested but no CUDA device is available")
        return lambda x: torch.tensor(x, device="cuda")
    raise ValueError(f"unknown backend {backend!r}; use 'numpy', 'cupy' or 'torch'")


_MAT_CACHE: dict = {}


def _gate_matrix(op):
    """Big-endian matrix (first qubit = most significant) of a Qiskit gate, for quimb."""
    key = (op.name, tuple(float(p) for p in op.params))
    if key not in _MAT_CACHE:
        from qiskit.quantum_info import Operator
        _MAT_CACHE[key] = np.asarray(Operator(op).reverse_qargs().data, dtype=np.complex128)
    return _MAT_CACHE[key]


def _to_numpy(x):
    if hasattr(x, "get"):
        return x.get()
    if hasattr(x, "detach"):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def _site_arrays(psi):
    """MPS tensors as (left, phys, right) numpy arrays (dummy bonds of size 1 at the ends)."""
    n = psi.L
    out = []
    for k in range(n):
        t = psi[k]
        left = [i for i in t.inds if k > 0 and i in psi[k - 1].inds]
        right = [i for i in t.inds if k < n - 1 and i in psi[k + 1].inds]
        arr = _to_numpy(t.transpose(*left, psi.site_ind(k), *right).data)
        if not left:
            arr = arr[None]
        if not right:
            arr = arr[..., None]
        out.append(np.asarray(arr, dtype=np.complex128))
    return out


def _transfer(X, A, O=None):
    if O is None:
        return np.einsum("ab,apc,bpd->cd", X, A, A.conj(), optimize=True)
    return np.einsum("ab,apc,qp,bqd->cd", X, A, O, A.conj(), optimize=True)


def _observables(lat, psi, site_of, g):
    """(<H_E>, <Z_p> in plaquette order) of an MPS, by left/right environment sweeps."""
    A = _site_arrays(psi)
    n = len(A)
    L = [np.ones((1, 1), complex)]
    for k in range(n):
        L.append(_transfer(L[-1], A[k]))
    R = [None] * (n + 1)
    R[n] = np.ones((1, 1), complex)
    for k in range(n - 1, -1, -1):
        R[k] = np.einsum("apc,bpd,cd->ab", A[k], A[k].conj(), R[k + 1], optimize=True)
    norm = np.real(L[n][0, 0])
    zsite = np.array([np.real(np.sum(_transfer(L[k], A[k], _Z) * R[k + 1]))
                      for k in range(n)])
    zsite = zsite / norm
    he = 0.0
    for p, q in lat.bonds:
        i, j = sorted((site_of[p], site_of[q]))
        X = _transfer(L[i], A[i], _Z)
        for k in range(i + 1, j):
            X = _transfer(X, A[k])
        zz = np.real(np.sum(_transfer(X, A[j], _Z) * R[j + 1])) / norm
        he += 1.0 - zz
    z = np.array([zsite[site_of[p]] for p in range(lat.n)])
    return 3.0 * g ** 2 / 16.0 * he, z


def simulate_mps(lat, g, dt, steps, order=1, abelian=False, initial=None, chi=256, cutoff=1e-12,
                 backend="numpy") -> dict:
    """Run trotter_circuit(lat, g, dt, steps, order, abelian, initial) as a CircuitMPS."""
    import quimb.tensor as qtn
    conv = _to_backend(backend)
    t0 = time.time()
    order_sites = snake_order(lat)
    site_of = {p: k for k, p in enumerate(order_sites)}
    circ = qtn.CircuitMPS(lat.n, max_bond=chi, cutoff=cutoff, to_backend=conv)

    def apply(qc):
        for inst in qc.data:
            if inst.operation.name == "barrier":
                continue
            qs = [site_of[qc.find_bit(q).index] for q in inst.qubits]
            circ.apply_gate(_gate_matrix(inst.operation), *qs)

    apply(trotter_circuit(lat, g, dt, 0, order, abelian, initial))   # preparation only
    step = trotter_step(lat, g, dt, order, abelian)
    he0, z0 = _observables(lat, circ.psi, site_of, g)
    ee, zs, err = [he0], [z0.tolist()], [0.0]
    for _ in range(steps):
        apply(step)
        he, z = _observables(lat, circ.psi, site_of, g)
        ee.append(he); zs.append(z.tolist())
        err.append(float(1.0 - _to_numpy(circ.fidelity_estimate())))
    psi = circ.psi
    fid = float(_to_numpy(circ.fidelity_estimate()))
    return {"n": lat.n, "chi": int(chi), "site_order": order_sites,
            "times": [k * dt for k in range(steps + 1)],
            "electric_energy": ee, "z": zs,
            "fidelity_estimate": fid, "error_estimate": 1.0 - fid,
            "error_estimate_per_step": err,
            "max_bond": int(psi.max_bond()), "wall_time_s": time.time() - t0, "mps": psi}


def mps_statevector(res) -> np.ndarray:
    """Dense normalised state in METHODS order b = sum_p n_p 2**p."""
    n = res["n"]
    v = _to_numpy(res["mps"].to_dense()).reshape((2,) * n)
    # site axis k holds plaquette order[k]; METHODS C-order axis of plaquette p is n-1-p
    perm = [0] * n
    for k, p in enumerate(res["site_order"]):
        perm[n - 1 - p] = k
    v = np.transpose(v, perm).reshape(-1)
    return v / np.linalg.norm(v)
