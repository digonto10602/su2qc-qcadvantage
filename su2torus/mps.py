"""MPS simulation of the Trotter circuits with quimb CircuitMPS (PLAN M5; API fixed by
tests/acceptance/m5). The circuit is exactly su2torus.trotter.trotter_circuit; gates are applied
one by one as explicit matrices (taken from Qiskit, so the conventions cannot drift), on an MPS
whose site k holds plaquette snake_order(lat)[k]."""
from __future__ import annotations

import time

import autoray as ar
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


def _site_array(psi, k):
    """Site tensor k as a (left, phys, right) array in its own backend (numpy/cupy/torch);
    dummy bonds of size 1 at the chain ends. No host copy."""
    n = psi.L
    t = psi[k]
    left = [i for i in t.inds if k > 0 and i in psi[k - 1].inds]
    right = [i for i in t.inds if k < n - 1 and i in psi[k + 1].inds]
    arr = t.transpose(*left, psi.site_ind(k), *right).data
    if not left:
        arr = arr[None]
    if not right:
        arr = arr[..., None]
    return arr


def _transfer(X, A, Z=False):
    """Left environment step X'_{cd} = sum X_{ab} A_{apc} conj(A)_{bpd}, with a Z = diag(1,-1)
    insertion on the physical index if Z. Backend-agnostic (autoray)."""
    T = ar.do("tensordot", X, A, axes=(0, 0))                    # (b, p, c)
    if Z:
        T = ar.do("concatenate", [T[:, :1], -T[:, 1:]], axis=1)
    return ar.do("tensordot", T, ar.do("conj", A), axes=([0, 1], [0, 1]))


def _close(X, R):
    return float(np.real(_to_numpy(ar.do("sum", X * R))))


def _observables(lat, psi, site_of, g):
    """(<H_E>, <Z_p> in plaquette order) of an MPS, by environment sweeps on the MPS's own
    backend (GPU arrays stay on the GPU). Memory: the n right environments (chi^2 each) plus
    O(1) left environments; site tensors are taken from the MPS as needed.
    Bonds are grouped by left site so each left environment is propagated once to all partners."""
    n = psi.L
    A0 = _site_array(psi, n - 1)
    ones = ar.do("ones", (1, 1), like=A0)
    R = [None] * (n + 1)
    R[n] = ones
    for k in range(n - 1, -1, -1):
        A = _site_array(psi, k)
        R[k] = ar.do("tensordot", ar.do("tensordot", A, R[k + 1], axes=(2, 0)),
                     ar.do("conj", A), axes=([1, 2], [1, 2]))
    norm = float(np.real(_to_numpy(R[0][0, 0])))
    partners: dict = {}
    for p, q in lat.bonds:
        i, j = sorted((site_of[p], site_of[q]))
        partners.setdefault(i, []).append(j)
    zsite = np.zeros(n)
    zzsum = 0.0
    L = ones
    for i in range(n):
        Ai = _site_array(psi, i)
        X = _transfer(L, Ai, Z=True)
        zsite[i] = _close(X, R[i + 1]) / norm
        k = i + 1
        for j in sorted(partners.get(i, [])):
            while k < j:
                X = _transfer(X, _site_array(psi, k))
                k += 1
            zzsum += _close(_transfer(X, _site_array(psi, j), Z=True), R[j + 1])
        L = _transfer(L, Ai)
        R[i] = None                                    # free right environments as we pass
    he = len(lat.bonds) - zzsum / norm
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
