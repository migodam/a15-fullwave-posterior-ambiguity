"""Self-contained 2D scalar Helmholtz direct solver for the A13 physics backend.

This module is a *bounded mechanical* forward-model backend for the
FlowMultiDecode work (A13).  It contains no research conclusions, no learned
components, and it never generates training data on import.  All computation is
local CPU NumPy/SciPy.

Conventions
-----------
Domain / grid
    D = [-0.5, 0.5]^2, uniform ``n x n`` cell-centre grid, cell side ``h = 1/n``,
    cell area ``h^2``.  Cell centres

        z_{i,j} = (-0.5 + (i + 0.5) h, -0.5 + (j + 0.5) h),

    flattened row-major with pixel index ``p = i*n + j`` (``i`` = x index).  The
    flattened grid has ``pixels = n*n`` entries.

Green function (background wavenumber ``k`` real positive)
    g(r, r') = (i/4) H_0^{(1)}(k |r - r'|).

Self-cell quadrature (the corrected v2 formula, copied verbatim from the
verified reference solver ``experiments/idea_loops/loop_2026-09-03_16-25-02/
experiment_pose_confounding_spectral_geometry/src/helmholtz.py``, SHA256
``ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585``)

    I_self = (i*pi*a/(2*k)) * H_1^{(1)}(k*a) - 1/k^2,      a = h/sqrt(pi),

    the exact integral of ``g`` over the equal-area disk of radius ``a``
    (area ``pi a^2 = h^2``).  The ``-1/k^2`` lower-endpoint term is required:
    without it the disk integral tends to ``1/k^2`` instead of ``0`` as
    ``h -> 0`` (see ``self_cell_green``).  The diagonal entry of ``GD`` is
    ``k^2 * I_self`` (``I_self`` already carries the cell-area factor).

Operators (dense ``complex128``)
    GD[p, q] = k^2 h^2 g(z_p, z_q)  off-diagonal,  diagonal = k^2 * I_self
    GS[a, q] = k^2 h^2 g(r_a, z_q)  (cell -> receiver propagator)
    incident[q, l] = exp(i k d_l . z_q)  (unit-amplitude plane wave, direction d_l)
    M(chi)   = I - diag(chi) GD
    J[:, l]  = solve(M, chi * incident[:, l])   (ONE LU reused for all illuminations)
    E[:, l]  = incident[:, l] + GD @ J[:, l]
    y[:, l]  = GS @ J[:, l]                     (scattered field at the receivers)
    rho      = k^2 * chi                        (contrast-source reflectivity)

Public API
    make_grid(n)                     -> (points, h)
    self_cell_green(k, h)            -> complex
    green_matrix(src, dst, k, ...)   -> complex ndarray
    Solver(n=32, k=2*pi, n_rx=32, n_illum=2, ...)
        .points, .GD, .GS, .incident, .geom, .indices, .receiver_indices
        .solve(chi)             -> dict (J, E, y, rho, residual, ...)
        .scene(seed, ood=False) -> real chi vector of length pixels
        .consistency(chi, J, E, y, indices=None) -> raw/normalized residual dict
    generate_dataset(path, count, seed=..., n=32, ...)   # NOT executed here

Numerics note: ``Solver.solve`` performs a single dense LU factorisation of
``M`` and reuses it for all illuminations.  It never forms a full SVD.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
from numpy.lib.format import open_memmap
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1

__all__ = [
    "SOURCE_PATH",
    "SOURCE_SHA256",
    "SELF_CELL_FORMULA",
    "make_grid",
    "self_cell_green",
    "green_matrix",
    "Solver",
    "generate_dataset",
]

# Provenance of the reference solver the corrected self-cell term is copied from.
SOURCE_PATH = (
    "experiments/idea_loops/loop_2026-09-03_16-25-02/"
    "experiment_pose_confounding_spectral_geometry/src/helmholtz.py"
)
SOURCE_SHA256 = "ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585"
SELF_CELL_FORMULA = (
    "I_self = (i*pi*a/(2*k))*H_1^(1)(k*a) - 1/k^2, a = h/sqrt(pi) "
    "(v2 corrected 2026-09-03: complete equal-area disk integral)"
)


# ---------------------------------------------------------------------------
# Grid and Green functions
# ---------------------------------------------------------------------------

def make_grid(n: int) -> tuple[np.ndarray, float]:
    """Return (points, h) for the D = [-0.5, 0.5]^2 ``n x n`` cell-centre grid.

    ``points`` has shape ``(n*n, 2)``; pixel ``p = i*n + j`` holds the cell
    centre ``(x_i, y_j)`` with row-major (C) flattening.
    """
    if not isinstance(n, (int, np.integer)) or int(n) < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}")
    n = int(n)
    h = 1.0 / n
    centres = -0.5 + (np.arange(n, dtype=float) + 0.5) * h
    X, Y = np.meshgrid(centres, centres, indexing="ij")
    points = np.empty((n * n, 2), dtype=float)
    points[:, 0] = X.ravel(order="C")
    points[:, 1] = Y.ravel(order="C")
    return points, h


def self_cell_green(k: float, h: float) -> complex:
    """Corrected equal-area disk integral of the Green function (area units).

    ``I_self = (i*pi*a/(2*k)) * H_1^{(1)}(k*a) - 1/k^2`` with ``a = h/sqrt(pi)``.

    Derivation of the ``-1/k^2`` term: with ``x = k*r`` and
    ``d[x H_1^{(1)}(x)]/dx = x H_0^{(1)}(x)``, the disk integral equals
    ``(i/4) 2*pi * (1/k^2) [x H_1^{(1)}(x)]_0^{k*a}``.  Because
    ``lim_{x->0} x H_1^{(1)}(x) = -2i/pi`` (not 0), the lower endpoint
    contributes ``(i/4) 2*pi * (1/k^2) * (2i/pi) = -1/k^2``.  Dropping it would
    leave an integral tending to ``1/k^2`` rather than ``0`` as ``h -> 0``.
    """
    k = float(k)
    a = float(h) / np.sqrt(np.pi)
    return (1j * np.pi * a / (2.0 * k)) * hankel1(1, k * a) - 1.0 / k**2


def _pair_diff(src: np.ndarray, dst: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """diff[d, s, :] = dst[d] - src[s]; R[d, s] = |diff[d, s, :]|."""
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    diff = dst[:, None, :] - src[None, :, :]
    return diff, np.linalg.norm(diff, axis=2)


def green_matrix(
    src: np.ndarray,
    dst: np.ndarray,
    k: float,
    include_self: bool = False,
    h: float | None = None,
) -> np.ndarray:
    """Point Green matrix ``G[d, s] = g(dst[d], src[s])``.

    Rows index the first argument (field/observation points ``dst``), columns
    the second argument (source/integration points ``src``).  Coincident pairs
    (``R = 0``) are set to 0; with ``include_self=True`` (requires ``h``) they
    are replaced by the disk-integrated ``self_cell_green(k, h)``, which already
    contains the cell-area factor, so callers must not multiply the diagonal by
    ``h^2`` again.
    """
    src = np.atleast_2d(np.asarray(src, dtype=float))
    dst = np.atleast_2d(np.asarray(dst, dtype=float))
    if src.shape[1] != 2 or dst.shape[1] != 2:
        raise ValueError("src and dst must have shape (..., 2)")
    _, R = _pair_diff(src, dst)
    safe = R > 0.0
    G = np.zeros(R.shape, dtype=np.complex128)
    if safe.any():
        G[safe] = (1j / 4.0) * hankel1(0, float(k) * R[safe])
    if include_self:
        if h is None:
            raise ValueError("include_self=True requires the cell side h")
        G[~safe] = self_cell_green(k, h)
    return G


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

class Solver:
    """Direct 2D scalar Helmholtz contrast-source solver on a fixed n x n grid.

    Illuminations are unit-amplitude plane waves; receivers are ``n_rx`` points
    equally spaced on a circle of radius ``rx_radius`` (> the domain half-width,
    so no receiver coincides with a grid cell).

    Geometry attributes (all fixed at construction, independent of ``chi``):
        n, k, n_rx, n_illum, pixels, h
        points           (pixels, 2)        cell centres
        GD               (pixels, pixels)   domain propagator (complex)
        GS               (n_rx, pixels)     cell -> receiver propagator (complex)
        incident         (pixels, n_illum)  plane-wave incident field (complex)
        receiver_pos     (n_rx, 2)          receiver coordinates
        illum_dir        (n_illum, 2)       plane-wave propagation directions
        geom             dict of the above geometry metadata
        indices          dict {"pixels","rx","illum"} full index ranges
        receiver_indices (n_rx,)            receiver ordering index
    """

    def __init__(
        self,
        n: int = 32,
        k: float = 2.0 * np.pi,
        n_rx: int = 32,
        n_illum: int = 2,
        rx_radius: float = 0.8,
        domain: tuple[float, float] = (-0.5, 0.5),
    ) -> None:
        if n_rx < 1 or n_illum < 1:
            raise ValueError("n_rx and n_illum must be positive")
        self.n = int(n)
        self.k = float(k)
        self.n_rx = int(n_rx)
        self.n_illum = int(n_illum)
        self.pixels = self.n * self.n
        self.rx_radius = float(rx_radius)
        self.domain = (float(domain[0]), float(domain[1]))

        self.points, self.h = make_grid(self.n)

        # Plane-wave illumination directions and incident field (pixels, n_illum).
        self.illum_angles = 2.0 * np.pi * np.arange(self.n_illum) / max(self.n_illum, 1)
        self.illum_dir = np.stack(
            [np.cos(self.illum_angles), np.sin(self.illum_angles)], axis=1
        )
        phase = self.points @ self.illum_dir.T  # (pixels, n_illum)
        self.incident = np.exp(1j * self.k * phase)

        # Receivers on a circle (fixed, illumination- and chi-independent).
        self.rx_angles = 2.0 * np.pi * (np.arange(self.n_rx) + 0.5) / self.n_rx
        self.receiver_pos = self.rx_radius * np.stack(
            [np.cos(self.rx_angles), np.sin(self.rx_angles)], axis=1
        )

        # Domain propagator GD (pixels, pixels) with corrected self-cell diagonal.
        raw = green_matrix(self.points, self.points, self.k)  # point values, diag = 0
        GD = (self.k**2) * (self.h**2) * raw
        np.fill_diagonal(GD, (self.k**2) * self_cell_green(self.k, self.h))
        self.GD = GD

        # Cell -> receiver propagator GS (n_rx, pixels); no coincident pairs.
        self.GS = (self.k**2) * (self.h**2) * green_matrix(
            self.points, self.receiver_pos, self.k
        )

        # Full index ranges + receiver ordering, for masking / subset selection.
        self.indices = {
            "pixels": np.arange(self.pixels),
            "rx": np.arange(self.n_rx),
            "illum": np.arange(self.n_illum),
        }
        self.receiver_indices = np.arange(self.n_rx)

        self.geom = {
            "domain": self.domain,
            "n": self.n,
            "pixels": self.pixels,
            "h": self.h,
            "k": self.k,
            "n_rx": self.n_rx,
            "n_illum": self.n_illum,
            "rx_radius": self.rx_radius,
            "rx_angles": self.rx_angles,
            "receiver_pos": self.receiver_pos,
            "illum_angles": self.illum_angles,
            "illum_dir": self.illum_dir,
            "receiver_indices": self.receiver_indices,
        }

    # -- public API ---------------------------------------------------------

    def solve(self, chi: np.ndarray) -> dict:
        """Solve the contrast-source system for one contrast ``chi``.

        Returns a dict with (shapes in the second column):
            chi          (pixels,)          input contrast (flattened grid)
            J            (pixels, n_illum)  induced current, one LU for all illuminations
            E            (pixels, n_illum)  total field
            y            (n_rx, n_illum)    scattered field at receivers
            rho          (pixels,)          k^2 * chi
            residual     (float)            ||M J - chi E_inc||_2 / ||chi E_inc||_2
            residual_abs (float)            ||M J - chi E_inc||_2

        ``chi`` may be real or complex; it is reshaped to a length-``pixels``
        vector.  All-zero ``chi`` gives ``J = 0``, ``y = 0``, ``E = incident``
        and a finite ``residual`` of exactly ``0.0`` (the denominator is zero).
        """
        chi = self._as_contrast(chi)
        rhs = chi[:, None] * self.incident
        M = np.eye(self.pixels, dtype=np.complex128) - chi[:, None] * self.GD
        lu, piv = lu_factor(M, check_finite=False)
        J = lu_solve((lu, piv), rhs, check_finite=False)
        E = self.incident + self.GD @ J
        y = self.GS @ J
        rho = (self.k**2) * chi

        dres = M @ J - rhs
        denom = float(np.linalg.norm(rhs))
        residual_abs = float(np.linalg.norm(dres))
        residual = residual_abs / denom if denom > 0.0 else 0.0
        return {
            "chi": chi,
            "J": J,
            "E": E,
            "y": y,
            "rho": rho,
            "residual": residual,
            "residual_abs": residual_abs,
        }

    def scene(self, seed: int, ood: bool = False) -> np.ndarray:
        """Return a random flattened contrast grid (length ``pixels``, real).

        The scene is the sum of 3-5 elliptical Gaussian blobs.  In-distribution
        (``ood=False``) the blob amplitude sign is *correlated with the blob
        size*: large (major axis >= 0.105) blobs are positive, small blobs are
        negative.  ``ood=True`` reverses that correlation (large negative, small
        positive).  Only ``chi`` is returned: no extra labels are invented.
        """
        xs = self.points[:, 0]
        ys = self.points[:, 1]
        chi = np.zeros(self.pixels, dtype=np.float64)
        for cx, cy, major, ratio, theta, amp in self._blob_params(seed, ood=ood):
            sx = major * np.sqrt(ratio)
            sy = major / np.sqrt(ratio)
            ct, st = np.cos(theta), np.sin(theta)
            dx = xs - cx
            dy = ys - cy
            xr = ct * dx + st * dy
            yr = -st * dx + ct * dy
            chi += amp * np.exp(-(xr**2 / (2.0 * sx**2) + yr**2 / (2.0 * sy**2)))
        return chi

    def _blob_params(self, seed: int, ood: bool = False) -> list[tuple]:
        """Blob parameters ``(cx, cy, major, ratio, theta, amp)`` for a scene.

        ``major`` is the semi-major axis, ``ratio = sx/sy`` the aspect ratio and
        ``amp`` the signed peak amplitude.  The sign encodes the amplitude/size
        correlation (reversed when ``ood=True``); this helper exposes the raw
        generator parameters so the correlation can be checked mechanically.
        """
        rng = np.random.default_rng(seed)
        n_blobs = int(rng.integers(3, 6))  # 3, 4 or 5
        params = []
        for _ in range(n_blobs):
            cx = float(rng.uniform(-0.28, 0.28))
            cy = float(rng.uniform(-0.28, 0.28))
            major = float(rng.uniform(0.05, 0.16))
            ratio = float(rng.uniform(0.65, 1.5))
            theta = float(rng.uniform(0.0, np.pi))
            amp = float(rng.uniform(0.5, 1.6))
            positive = major >= 0.105  # size/amplitude correlation
            if ood:
                positive = not positive
            if not positive:
                amp = -amp
            params.append((cx, cy, major, ratio, theta, amp))
        return params

    def consistency(
        self,
        chi: np.ndarray,
        J: np.ndarray,
        E: np.ndarray,
        y: np.ndarray,
        indices: dict | None = None,
    ) -> dict:
        """Raw and normalized physical-consistency residuals for a decode.

        Nothing is replaced: the supplied ``J``, ``E``, ``y`` are used as-is.
        Each residual is reported as a raw (complex) array plus a scalar norm and
        a normalized value (norm of the residual divided by the norm of its
        reference term).  Components:

            j_minus_chiE        J - chi * E               (contrast-source relation)
            E_minus_inc_GDJ     E - incident - GD @ J     (field equation)
            forward_minus_y     GS @ J - y                (data misfit of the given J)
            J_minus_resolved    J - J(chi)                (given J vs. re-solved J(chi))
            E_minus_resolved    E - E(chi)                (given E vs. re-solved E(chi))

        The last two re-solve the forward problem for the *same* ``chi`` purely as
        a diagnostic; they do not overwrite or return a replacement solution.

        ``indices`` optionally restricts the comparison, e.g.
        ``{"rx": np.arange(8), "illum": [0]}``; it is a dict with optional keys
        ``"pixels"``, ``"rx"``, ``"illum"`` (defaults: full ranges).
        """
        chi_v = self._as_contrast(chi)
        J = np.asarray(J, dtype=np.complex128)
        E = np.asarray(E, dtype=np.complex128)
        y = np.asarray(y, dtype=np.complex128)
        if J.shape != (self.pixels, self.n_illum):
            raise ValueError(f"J must have shape {(self.pixels, self.n_illum)}")
        if E.shape != J.shape:
            raise ValueError("E must have the same shape as J")
        if y.shape != (self.n_rx, self.n_illum):
            raise ValueError(f"y must have shape {(self.n_rx, self.n_illum)}")

        idx = self._resolve_indices(indices)
        pix, rx, ill = idx["pixels"], idx["rx"], idx["illum"]

        Jc = J[np.ix_(pix, ill)]
        Ec = E[np.ix_(pix, ill)]
        chia = chi_v[pix][:, None]
        inc = self.incident[np.ix_(pix, ill)]
        GDc = self.GD[pix, :]

        # 1. Contrast-source (current) relation.
        j_minus_chiE = Jc - chia * Ec
        # 2. Field equation.
        E_minus_inc_GDJ = Ec - inc - GDc @ J[:, ill]
        # 3. Forward data misfit of the supplied J.
        forward_minus_y = self.GS[rx, :] @ J[:, ill] - y[np.ix_(rx, ill)]
        # 4. Supplied fields vs. a re-solve at the same chi (diagnostic only).
        sol = self.solve(chi_v)
        J_minus_resolved = Jc - sol["J"][np.ix_(pix, ill)]
        E_minus_resolved = Ec - sol["E"][np.ix_(pix, ill)]

        def _norm(a):
            return float(np.linalg.norm(np.ravel(a)))

        def _entry(res, ref):
            d = _norm(ref)
            return {
                "raw": res,
                "norm": _norm(res),
                "relative": _norm(res) / max(d, 1e-12),
            }

        return {
            "j_minus_chiE": _entry(j_minus_chiE, chia * Ec),
            "E_minus_inc_GDJ": _entry(E_minus_inc_GDJ, Ec),
            "forward_minus_y": _entry(forward_minus_y, y[np.ix_(rx, ill)]),
            "material_forward_minus_y": _entry(sol["y"][np.ix_(rx, ill)] - y[np.ix_(rx, ill)], y[np.ix_(rx, ill)]),
            "J_minus_resolved": _entry(J_minus_resolved, sol["J"][np.ix_(pix, ill)]),
            "E_minus_resolved": _entry(E_minus_resolved, sol["E"][np.ix_(pix, ill)]),
            "indices": {key: val.copy() for key, val in idx.items()},
        }

    def flat_index(
        self, rx: np.ndarray | None = None, illum: np.ndarray | None = None
    ) -> np.ndarray:
        """Row-major flat indices into ``y.reshape(-1)`` for a receiver/illum mask."""
        rows = self.receiver_indices if rx is None else np.asarray(rx, dtype=int).reshape(-1)
        cols = np.arange(self.n_illum) if illum is None else np.asarray(illum, dtype=int).reshape(-1)
        rr, ll = np.meshgrid(rows, cols, indexing="ij")
        return np.ravel_multi_index((rr.ravel(), ll.ravel()), (self.n_rx, self.n_illum))

    # -- internals ----------------------------------------------------------

    def _as_contrast(self, chi: np.ndarray) -> np.ndarray:
        chi = np.asarray(chi)
        if chi.size != self.pixels:
            raise ValueError(f"chi must have {self.pixels} entries, got {chi.size}")
        if np.iscomplexobj(chi):
            return chi.astype(np.complex128).reshape(-1)
        return chi.astype(np.float64).reshape(-1)

    def _resolve_indices(self, indices: dict | None) -> dict:
        indices = indices or {}

        def _sel(key: str, size: int) -> np.ndarray:
            val = indices.get(key, None)
            if val is None:
                return np.arange(size)
            return np.asarray(val, dtype=int).reshape(-1)

        return {
            "pixels": _sel("pixels", self.pixels),
            "rx": _sel("rx", self.n_rx),
            "illum": _sel("illum", self.n_illum),
        }


# ---------------------------------------------------------------------------
# Dataset writer (implemented but NOT executed in this subtask)
# ---------------------------------------------------------------------------

def generate_dataset(
    path: str,
    count: int,
    seed: int = 0,
    n: int = 32,
    k: float = 2.0 * np.pi,
    n_rx: int = 32,
    n_illum: int = 2,
    chunk: int = 64,
    ood_fraction: float = 0.0,
) -> dict:
    """Write a labelled dataset in chunks with ``open_memmap`` plus a JSON manifest.

    Not called by this subtask (the Phase A gate has not passed); provided so the
    dataset layout is fixed and reviewable.  Arrays are written to

        chi.npy  (count, pixels)          float32
        y.npy    (count, n_rx, n_illum)   complex64
        J.npy    (count, pixels, n_illum) complex64
        E.npy    (count, pixels, n_illum) complex64
        rho.npy  (count, pixels)          float32

    and ``manifest.json`` records the geometry, provenance (``SOURCE_SHA256``,
    ``SELF_CELL_FORMULA``) and per-array shape/dtype/filename.  Scene ``i`` uses
    ``seed + i`` so the stream is stable under a change of ``count``.
    """
    solver = Solver(n=n, k=k, n_rx=n_rx, n_illum=n_illum)
    os.makedirs(path, exist_ok=True)
    count = int(count)
    chunk = max(1, int(chunk))
    pixels = solver.pixels
    wmode = "w+"

    chi_mm = open_memmap(os.path.join(path, "chi.npy"), mode=wmode,
                         dtype=np.float32, shape=(count, pixels))
    y_mm = open_memmap(os.path.join(path, "y.npy"), mode=wmode,
                       dtype=np.complex64, shape=(count, n_rx, n_illum))
    J_mm = open_memmap(os.path.join(path, "J.npy"), mode=wmode,
                       dtype=np.complex64, shape=(count, pixels, n_illum))
    E_mm = open_memmap(os.path.join(path, "E.npy"), mode=wmode,
                       dtype=np.complex64, shape=(count, pixels, n_illum))
    rho_mm = open_memmap(os.path.join(path, "rho.npy"), mode=wmode,
                         dtype=np.float32, shape=(count, pixels))

    started = time.perf_counter()
    diagnostics = np.zeros((count, 3), dtype=np.float64)
    for start in range(0, count, chunk):
        stop = min(start + chunk, count)
        for i in range(start, stop):
            ood = ood_fraction > 0.0 and ((i % 100) / 100.0) < ood_fraction
            chi = solver.scene(seed + i, ood=ood)
            try:
                out = solver.solve(chi)
                if not all(np.isfinite(out[key]).all() for key in ("J", "E", "y")):
                    raise FloatingPointError("nonfinite forward solution")
                if out["residual"] > 1e-8:
                    raise FloatingPointError(f"forward residual {out['residual']} exceeds 1e-8")
            except Exception as exc:
                np.savez_compressed(os.path.join(path, f"failed_scene_{i}.npz"), chi=chi, seed=seed+i)
                with open(os.path.join(path, "failures.jsonl"), "a") as fh:
                    fh.write(json.dumps({"index":i,"seed":seed+i,"error":str(exc)})+"\n")
                raise
            born = solver.GS @ (chi[:, None] * solver.incident)
            diagnostics[i] = [out["residual"], np.linalg.norm(born-out["y"])/max(np.linalg.norm(out["y"]),1e-12), np.abs(out["E"]).max()]
            chi_mm[i] = chi.astype(np.float32)
            y_mm[i] = out["y"].astype(np.complex64)
            J_mm[i] = out["J"].astype(np.complex64)
            E_mm[i] = out["E"].astype(np.complex64)
            rho_mm[i] = out["rho"].astype(np.float32)
        for mm in (chi_mm, y_mm, J_mm, E_mm, rho_mm):
            mm.flush()
        print(f"physics generation {stop}/{count} elapsed={time.perf_counter()-started:.1f}s", flush=True)

    np.savez_compressed(os.path.join(path,"solver_diagnostics.npz"), values=diagnostics,
                        columns=np.array(["relative_linear_residual","born_relative_discrepancy","max_abs_E"]))

    manifest = {
        "count": count,
        "seed": int(seed),
        "n": solver.n,
        "pixels": pixels,
        "k": solver.k,
        "n_rx": solver.n_rx,
        "n_illum": solver.n_illum,
        "h": solver.h,
        "domain": list(solver.domain),
        "rx_radius": solver.rx_radius,
        "source_path": SOURCE_PATH,
        "source_sha256": SOURCE_SHA256,
        "self_cell_formula": SELF_CELL_FORMULA,
        "chunk": chunk,
        "ood_fraction": ood_fraction,
        "generation_seconds": time.perf_counter()-started,
        "max_linear_residual": float(diagnostics[:,0].max()),
        "born_discrepancy_median": float(np.median(diagnostics[:,1])),
        "failure_policy": "retain failed scene and error, stop without completion marker; no silent replacement",
        "arrays": {
            "chi": {"file": "chi.npy", "shape": [count, pixels], "dtype": "float32"},
            "y": {"file": "y.npy", "shape": [count, n_rx, n_illum], "dtype": "complex64"},
            "J": {"file": "J.npy", "shape": [count, pixels, n_illum], "dtype": "complex64"},
            "E": {"file": "E.npy", "shape": [count, pixels, n_illum], "dtype": "complex64"},
            "rho": {"file": "rho.npy", "shape": [count, pixels], "dtype": "float32"},
        },
        "note": "mechanical forward-model dataset; no research conclusions",
    }
    with open(os.path.join(path, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return manifest
