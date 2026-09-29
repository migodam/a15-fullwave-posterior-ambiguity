"""A14 full-wave tangent actions using the validated A13 discrete solver.

The canonical state is real signed contrast chi.  All complex fields are
computed on CPU in complex128.  The same LU is reused for every illumination
and every tangent right-hand side at a scene.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator

from physics_backend import Solver, green_matrix


def pack_complex(x: np.ndarray) -> np.ndarray:
    """[Re(flat x), Im(flat x)] with a real Euclidean inner product."""
    v = np.asarray(x).reshape(-1)
    return np.concatenate((v.real, v.imag))


def unpack_complex(v: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    a = np.asarray(v, dtype=float).reshape(-1)
    n = int(np.prod(shape))
    if a.size != 2 * n:
        raise ValueError("packed complex vector has the wrong dimension")
    return (a[:n] + 1j * a[n:]).reshape(shape)


@dataclass(frozen=True)
class State:
    chi: np.ndarray
    lu: tuple[np.ndarray, np.ndarray]
    scene_hash: str
    condition_estimate: float | None


class PhysicalModel:
    def __init__(self, n: int = 32, k: float = 2 * np.pi,
                 n_rx: int = 32, rx_radius: float = .8):
        self.backend = Solver(n=n, k=k, n_rx=n_rx, n_illum=2,
                              rx_radius=rx_radius)
        self.n = n
        self.pixels = n * n
        self.k = k
        self.points = self.backend.points
        self.GD = self.backend.GD
        self.default_receivers = self.backend.receiver_pos

    def scene(self, seed: int, ood: bool = False) -> np.ndarray:
        return self.backend.scene(seed, ood=ood)

    def factorize(self, chi: np.ndarray, estimate_condition: bool = False) -> State:
        chi = np.asarray(chi, dtype=float).reshape(-1)
        if chi.size != self.pixels or not np.isfinite(chi).all():
            raise ValueError("chi must be a finite real contrast image")
        mat = np.eye(self.pixels, dtype=np.complex128) - chi[:, None] * self.GD
        condition = float(np.linalg.cond(mat)) if estimate_condition else None
        lu = lu_factor(mat, check_finite=False)
        key = hashlib.sha256(chi.tobytes()).hexdigest()
        return State(chi.copy(), lu, key, condition)

    def incident(self, angles: np.ndarray) -> np.ndarray:
        angles = np.atleast_1d(np.asarray(angles, dtype=float))
        directions = np.stack((np.cos(angles), np.sin(angles)), axis=1)
        return np.exp(1j * self.k * (self.points @ directions.T))

    def receiver_operator(self, receiver_positions: np.ndarray | None = None) -> np.ndarray:
        if receiver_positions is None:
            return self.backend.GS
        locations = np.asarray(receiver_positions, dtype=float)
        return (self.k**2 * self.backend.h**2 *
                green_matrix(self.points, locations, self.k))

    def shifted_receivers(self) -> np.ndarray:
        m = self.backend.n_rx
        angles = 2 * np.pi * (np.arange(m) + 1.) / m
        return self.backend.rx_radius * np.stack((np.cos(angles), np.sin(angles)), 1)

    def solve(self, state: State, illumination_angles: np.ndarray,
              receiver_positions: np.ndarray | None = None) -> dict:
        inc = self.incident(illumination_angles)
        j = lu_solve(state.lu, state.chi[:, None] * inc, check_finite=False)
        e = inc + self.GD @ j
        gs = self.receiver_operator(receiver_positions)
        y = gs @ j
        if not all(np.isfinite(a).all() for a in (j, e, y)):
            raise FloatingPointError("nonfinite physical solve")
        return {"chi": state.chi, "J": j, "E": e, "y": y,
                "incident": inc, "receiver_operator": gs,
                "linear_residual": float(np.linalg.norm(j-state.chi[:, None]*e)
                                         / max(np.linalg.norm(j), 1e-12))}

    def jvp(self, state: State, delta_chi: np.ndarray,
            illumination_angles: np.ndarray,
            task: str = "y", receiver_positions: np.ndarray | None = None,
            receiver_indices: np.ndarray | None = None,
            cached_fields: dict | None = None) -> np.ndarray:
        """Real material JVP.  Direction shape (P,) or (P,K)."""
        directions = np.asarray(delta_chi, dtype=float)
        single = directions.ndim == 1
        if single:
            directions = directions[:, None]
        if directions.ndim != 2 or directions.shape[0] != self.pixels:
            raise ValueError("delta_chi must have shape (pixels,) or (pixels,K)")
        fields = cached_fields or self.solve(state, illumination_angles, receiver_positions)
        e = fields["E"]
        rhs = (e[:, :, None] * directions[:, None, :]).reshape(self.pixels, -1)
        dj = lu_solve(state.lu, rhs, check_finite=False).reshape(self.pixels, e.shape[1], -1)
        if task == "J":
            ans = dj
        elif task == "E":
            ans = np.einsum("ij,jsk->isk", self.GD, dj, optimize=True)
        elif task == "y":
            gs = fields["receiver_operator"]
            if receiver_indices is not None:
                gs = gs[np.asarray(receiver_indices, dtype=int)]
            ans = np.einsum("ij,jsk->isk", gs, dj, optimize=True)
        else:
            raise ValueError(task)
        return ans[..., 0] if single else ans

    def vjp(self, state: State, cotangent: np.ndarray,
            illumination_angles: np.ndarray,
            task: str = "y", receiver_positions: np.ndarray | None = None,
            receiver_indices: np.ndarray | None = None,
            cached_fields: dict | None = None) -> np.ndarray:
        """Real adjoint: <cotangent, JVP(v)>_R = <VJP(cotangent), v>_R."""
        fields = cached_fields or self.solve(state, illumination_angles, receiver_positions)
        c = np.asarray(cotangent, dtype=np.complex128)
        if c.ndim != 2 or c.shape[1] != fields["E"].shape[1]:
            raise ValueError("cotangent shape mismatch")
        if task == "J":
            rhs = c
        elif task == "E":
            rhs = self.GD.conj().T @ c
        elif task == "y":
            gs = fields["receiver_operator"]
            if receiver_indices is not None:
                gs = gs[np.asarray(receiver_indices, dtype=int)]
            rhs = gs.conj().T @ c
        else:
            raise ValueError(task)
        q = lu_solve(state.lu, rhs, trans=2, check_finite=False)
        return np.real(np.sum(np.conj(q) * fields["E"], axis=1))

    def dense_jacobian(self, state: State, directions: np.ndarray,
                       illumination_angles: np.ndarray, task: str = "y",
                       receiver_positions: np.ndarray | None = None,
                       receiver_indices: np.ndarray | None = None,
                       block: int = 128) -> np.ndarray:
        directions = np.asarray(directions, dtype=float)
        fields = self.solve(state, illumination_angles, receiver_positions)
        out = []
        for lo in range(0, directions.shape[1], block):
            out.append(self.jvp(state, directions[:, lo:lo+block], illumination_angles,
                                task, receiver_positions, receiver_indices, fields))
        return np.concatenate(out, axis=-1)

    def linear_operator(self, state: State, illumination_angles: np.ndarray,
                        task: str = "y", receiver_positions: np.ndarray | None = None,
                        receiver_indices: np.ndarray | None = None) -> LinearOperator:
        fields = self.solve(state, illumination_angles, receiver_positions)
        shape = self.jvp(state, np.zeros(self.pixels), illumination_angles,
                         task, receiver_positions, receiver_indices, fields).shape
        def mv(v):
            return pack_complex(self.jvp(state, v, illumination_angles, task,
                                         receiver_positions, receiver_indices, fields))
        def rmv(c):
            return self.vjp(state, unpack_complex(c, shape), illumination_angles,
                            task, receiver_positions, receiver_indices, fields)
        return LinearOperator((2 * int(np.prod(shape)), self.pixels), matvec=mv,
                              rmatvec=rmv, dtype=float)


def physical_readout(chi: np.ndarray, conditions: dict,
                     model: PhysicalModel) -> dict:
    """No measurement enters the readout; J/E/y share one material realization."""
    return model.solve(model.factorize(chi), conditions["illumination_angles"],
                       conditions.get("receiver_positions"))
