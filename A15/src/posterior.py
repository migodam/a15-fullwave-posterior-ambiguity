"""Bounded-uniform 3D reference posterior, task-risk and tempered-SMC utilities.

This module is a *bounded mechanical* numerical toolkit for the A15
full-wave bright/dark task-ambiguity study.  It contains no physics, no
research conclusions and no learned components.  The forward response ``F``,
the scalar task readout ``Q`` and the log-likelihood are all supplied by the
caller; nothing here decides whether an ambiguity is physically real.

Design contract (from the A15 worker task):

* Prior: independent bounded uniform, ``theta = (s, Delta, gamma)`` on a box
  ``[a_i, b_i]``.  The reference grid is a tensor Gauss-Legendre rule; its
  weights are positive and normalised so that ``sum(w) == 1`` (pure prior
  probability weights, not quadrature weights in physical units).
* Likelihood: the measurement is complex, but the noise is specified
  independently on the **real and imaginary parts**, with one sigma for every
  *real* component.  For each complex channel, real and imaginary residuals
  are independent Gaussian.
* Normalisation always goes through :func:`scipy.special.logsumexp`.  Failed
  quadrature nodes (non-finite log-weight, or zero prior support = "missing")
  are **rejected**: the default is to raise, and even when rejection is
  explicitly allowed the node is dropped only when its *prior* mass is below a
  tolerance.  Failures are never silently renormalised away -- that would
  fabricate a posterior over a failed rule.
* SMC: adaptive tempering from an independent standard-Gaussian base reached
  by the exact uniform CDF transform ``theta_i = a_i + (b_i - a_i) Phi(u_i)``,
  with pCN rejuvenation and systematic (ESS-based) resampling.  This is a
  generic tempered sampler; it makes **no physical claim**.

All randomness is drawn from an explicit ``numpy.random.Generator``; results
are reproducible from ``seed``.
"""

from __future__ import annotations

import json
import os
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.special import logsumexp, ndtr, ndtri

__all__ = [
    "DEFAULT_CONFIG_PATH",
    "load_reference_tolerances",
    "gauss_legendre_box_nodes",
    "uniform_cdf_transform",
    "uniform_inverse_cdf",
    "complex_gaussian_log_likelihood",
    "posterior_log_weights",
    "normalize_log_weights",
    "ess",
    "conditional_ess",
    "systematic_resample",
    "weighted_mean",
    "weighted_variance",
    "weighted_std",
    "weighted_covariance",
    "weighted_quantiles",
    "weighted_credible_interval",
    "readout_statistics",
    "theta_statistics",
    "predictive_moments",
    "branch_weights",
    "posterior_summary",
    "population_risk",
    "probability_sd_above",
    "paired_bootstrap_risk",
    "compare_posterior_to_reference",
    "check_convergence_metrics",
    "adaptive_tempered_smc",
    "TemperedSMC",
]

# configs/discovery.json sits one level up from src/.
DEFAULT_CONFIG_PATH = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "configs", "discovery.json")
)

# Tolerance keys read from discovery.json for reference-grid comparison.
_REFERENCE_TOL_KEYS = (
    "reference_mean_tol",
    "reference_variance_rel_tol",
    "reference_variance_abs_tol",
    "reference_quantile_tol",
    "reference_branch_tol",
    "mesh_q_tol",
    "mesh_whitened_rms_tol",
    "mesh_whitened_max_tol",
    "linear_residual_tol",
)


# ---------------------------------------------------------------------------
# Discovery-config tolerances
# ---------------------------------------------------------------------------

def load_reference_tolerances(config_path: str | None = None) -> dict:
    """Load the frozen reference-grid tolerances from ``configs/discovery.json``.

    Returns a dict with the tolerance keys above plus ``quadrature_orders``,
    ``primary_noise`` and ``bootstrap_replicates``.  Missing keys raise, so a
    changed config can never silently loosen the comparison.
    """
    path = DEFAULT_CONFIG_PATH if config_path is None else config_path
    with open(path, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    tols: dict = {}
    for key in _REFERENCE_TOL_KEYS:
        if key not in cfg:
            raise KeyError(f"discovery.json is missing tolerance key {key!r}")
        tols[key] = float(cfg[key])
    tols["quadrature_orders"] = [int(o) for o in cfg.get("quadrature_orders", [])]
    tols["primary_noise"] = float(cfg.get("primary_noise", np.nan))
    tols["bootstrap_replicates"] = int(cfg.get("bootstrap_replicates", 0))
    tols["config_path"] = path
    return tols


# ---------------------------------------------------------------------------
# Tensor Gauss-Legendre rule for an independent bounded uniform prior
# ---------------------------------------------------------------------------

def gauss_legendre_box_nodes(bounds: Sequence | np.ndarray, order: int | Sequence[int]) -> dict:
    """Tensor Gauss-Legendre rule on a box ``[a_i, b_i]`` for a uniform prior.

    ``bounds`` is ``(d, 2)`` (or a length-``d`` sequence of ``(a_i, b_i)``);
    ``order`` is a scalar or per-dimension sequence.  Returns::

        {"nodes": (K, d), "weights": (K,), "log_weights": (K,), "raw_weights": (K,),
         "volume": float, "per_dim_nodes": [...], "per_dim_weights": [...],
         "order": [n_i]}

    ``weights`` are **normalised prior probabilities** (``K = prod(n_i)``,
    ``weights > 0``, ``sum(weights) == 1``).  ``raw_weights`` are the underlying
    quadrature weights in physical units (``sum == volume``).
    """
    b = np.asarray(bounds, dtype=float)
    if b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("bounds must be shape (d, 2)")
    if np.any(b[:, 1] <= b[:, 0]):
        raise ValueError("every box bound must satisfy b_i > a_i")
    d = b.shape[0]
    orders = [int(order)] * d if np.isscalar(order) else [int(o) for o in order]
    if len(orders) != d or any(o < 1 for o in orders):
        raise ValueError("order must be a positive scalar or a length-d sequence")

    per_nodes: list[np.ndarray] = []
    per_w: list[np.ndarray] = []
    for i in range(d):
        x, w = np.polynomial.legendre.leggauss(orders[i])  # x in [-1,1], w>0, sum(w)=2
        a_i, b_i = b[i, 0], b[i, 1]
        per_nodes.append(0.5 * (b_i - a_i) * x + 0.5 * (a_i + b_i))
        per_w.append(0.5 * (b_i - a_i) * w)

    grid_nodes = np.meshgrid(*per_nodes, indexing="ij")
    nodes = np.stack([g.ravel() for g in grid_nodes], axis=1)
    raw = np.ones(nodes.shape[0], dtype=float)
    for g in np.meshgrid(*per_w, indexing="ij"):
        raw = raw * g.ravel()

    volume = float(np.prod(b[:, 1] - b[:, 0]))
    weights = raw / volume
    return {
        "nodes": nodes,
        "weights": weights,
        "log_weights": np.log(weights),
        "raw_weights": raw,
        "volume": volume,
        "per_dim_nodes": per_nodes,
        "per_dim_weights": per_w,
        "order": orders,
    }


# ---------------------------------------------------------------------------
# Exact uniform <-> standard-normal CDF transform
# ---------------------------------------------------------------------------

def uniform_cdf_transform(u: np.ndarray, bounds: Sequence | np.ndarray) -> np.ndarray:
    """Exact map ``theta_i = a_i + (b_i - a_i) * Phi(u_i)`` from ``N(0, I)``.

    ``u`` is ``(N, d)`` (or ``(d,)``); the image lies exactly in the box, so all
    samples respect the physical prior family.  ``Phi`` is the standard normal
    CDF (``scipy.special.ndtr``).
    """
    b = np.asarray(bounds, dtype=float)
    if b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("bounds must be shape (d, 2)")
    u = np.asarray(u, dtype=float)
    if u.shape[-1] != b.shape[0]:
        raise ValueError("last axis of u must match the number of box dimensions")
    return b[:, 0] + (b[:, 1] - b[:, 0]) * ndtr(u)


def uniform_inverse_cdf(theta: np.ndarray, bounds: Sequence | np.ndarray) -> np.ndarray:
    """Inverse of :func:`uniform_cdf_transform`: ``u_i = Phi^{-1}(F_i(theta_i))``.

    Probabilities are clipped to ``[eps, 1 - eps]`` so exact box endpoints map to
    finite values instead of ``+/-inf``.
    """
    b = np.asarray(bounds, dtype=float)
    if b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("bounds must be shape (d, 2)")
    theta = np.asarray(theta, dtype=float)
    if theta.shape[-1] != b.shape[0]:
        raise ValueError("last axis of theta must match the number of box dimensions")
    p = (theta - b[:, 0]) / (b[:, 1] - b[:, 0])
    eps = float(np.finfo(float).eps)
    p = np.clip(p, eps, 1.0 - eps)
    return ndtri(p)


# ---------------------------------------------------------------------------
# Likelihood and posterior weights
# ---------------------------------------------------------------------------

def _sigma_per_real_component(sigma, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Split ``sigma`` into ``(sigma_real, sigma_imag)``, each length ``m``."""
    s = np.asarray(sigma, dtype=float)
    if s.ndim == 0:
        sr = si = np.full(m, float(s))
    elif s.shape == (m,):
        sr = si = s
    elif s.shape == (2 * m,):
        sr, si = s[:m], s[m:]
    else:
        raise ValueError(f"sigma must be scalar, ({m},) or ({2 * m},); got {s.shape}")
    if np.any(~np.isfinite(sr)) or np.any(~np.isfinite(si)) or np.any(sr <= 0) or np.any(si <= 0):
        raise ValueError("sigma entries must be finite and strictly positive")
    return sr, si


def complex_gaussian_log_likelihood(
    pred: np.ndarray,
    obs: np.ndarray,
    sigma,
    *,
    include_const: bool = True,
) -> np.ndarray:
    """Independent real/imag Gaussian log-likelihood for complex data.

    ``pred`` is ``(K, m)`` (or ``(m,)``) complex predictions, ``obs`` is ``(m,)``
    complex, and ``sigma`` sets one noise std for **every real component**:

    * scalar     -> same sigma for all real and imaginary parts;
    * ``(m,)``   -> ``sigma[j]`` used for both Re and Im of channel ``j``;
    * ``(2*m,)`` -> first ``m`` for Re, next ``m`` for Im.

    Returns ``(K,)`` log-likelihood with real and imaginary residuals treated as
    independent (``-0.5 sum (r/sigma)^2`` plus, when ``include_const``, the
    ``-sum log(2*pi*sigma^2)/2`` normalisation).
    """
    pred = np.atleast_2d(np.asarray(pred, dtype=complex))
    obs = np.asarray(obs, dtype=complex).reshape(-1)
    if pred.shape[1] != obs.shape[0]:
        raise ValueError("pred and obs disagree on the number of channels")
    m = obs.shape[0]
    sr, si = _sigma_per_real_component(sigma, m)
    res = pred - obs[None, :]
    rr, ri = res.real, res.imag
    ll = -0.5 * np.sum((rr / sr) ** 2, axis=1) - 0.5 * np.sum((ri / si) ** 2, axis=1)
    if include_const:
        ll = ll - 0.5 * np.sum(np.log(2.0 * np.pi * sr**2)) - 0.5 * np.sum(np.log(2.0 * np.pi * si**2))
    return ll


def posterior_log_weights(
    log_likelihood: np.ndarray,
    prior_weights: np.ndarray | None = None,
) -> np.ndarray:
    """Unnormalised posterior log-weights ``log(prior) + logL`` per node.

    ``prior_weights`` defaults to uniform ``1/K``.  Non-finite likelihood entries
    propagate (``nan``) and are caught later by :func:`normalize_log_weights`.
    """
    ll = np.asarray(log_likelihood, dtype=float).ravel()
    k = ll.size
    if prior_weights is None:
        log_prior = np.full(k, -np.log(k))
    else:
        pw = np.asarray(prior_weights, dtype=float).ravel()
        if pw.size != k:
            raise ValueError("prior_weights and log_likelihood must have the same length")
        if np.any(pw < 0) or not np.isfinite(pw).all() or pw.sum() <= 0:
            raise ValueError("prior_weights must be non-negative, finite and have positive sum")
        log_prior = np.where(pw > 0, np.log(np.where(pw > 0, pw, 1.0)), -np.inf)
    return log_prior + ll


def normalize_log_weights(
    log_weights: np.ndarray,
    prior_weights: np.ndarray | None = None,
    *,
    strict: bool = True,
    invalid_mass_tol: float = 1e-12,
) -> dict:
    """Normalise log-weights with ``logsumexp``, rejecting failed nodes.

    ``log_weights`` must already be the **unnormalised log posterior**
    (``log prior + log-likelihood``), as returned by :func:`posterior_log_weights`;
    the prior is *not* added again here.  ``prior_weights`` is used only to flag
    "missing" nodes (zero prior support) and to report the rejected prior mass.

    A node is invalid if its log-weight is non-finite (``nan``/``+/-inf``) or if
    its prior weight is zero.  Failed nodes are never silently renormalised away:

    * ``strict=True`` (default): any invalid node raises ``ValueError``;
    * ``strict=False``: invalid nodes are dropped **only** when their total prior
      mass ``<= invalid_mass_tol``; otherwise ``ValueError`` is raised.

    Returns ``{"weights", "log_weights", "log_evidence", "ess", "info"}`` where
    ``info`` records ``rejected_indices``, ``rejected_prior_mass`` and
    ``n_rejected``; rejected entries get weight exactly 0.
    """
    lw = np.asarray(log_weights, dtype=float).ravel()
    k = lw.size
    if prior_weights is None:
        pw = np.full(k, 1.0 / k)
    else:
        pw = np.asarray(prior_weights, dtype=float).ravel()
        if pw.size != k:
            raise ValueError("prior_weights and log_weights must have the same length")
        if np.any(pw < 0) or not np.isfinite(pw).all():
            raise ValueError("prior_weights must be non-negative and finite")

    invalid = (~np.isfinite(lw)) | (pw <= 0.0)
    info = {"rejected_indices": np.array([], dtype=int), "rejected_prior_mass": 0.0, "n_rejected": 0}
    if invalid.any():
        idx = np.nonzero(invalid)[0]
        mass = float(pw[invalid].sum())
        if strict or mass > invalid_mass_tol:
            raise ValueError(
                f"{idx.size} non-finite/missing posterior node(s) at {idx[:8].tolist()}"
                f"{'...' if idx.size > 8 else ''} with prior mass {mass:.3e}; "
                "refusing to renormalise a failed quadrature"
            )
        info = {"rejected_indices": idx, "rejected_prior_mass": mass, "n_rejected": int(idx.size)}

    keep = ~invalid
    log_num = np.where(keep, lw, -np.inf)
    log_z = logsumexp(log_num)
    if not np.isfinite(log_z):
        raise ValueError("posterior normalisation is not finite")
    log_w = log_num - log_z
    weights = np.exp(log_w)
    return {
        "weights": weights,
        "log_weights": log_w,
        "log_evidence": float(log_z),
        "ess": ess(weights),
        "info": info,
    }


# ---------------------------------------------------------------------------
# ESS and resampling
# ---------------------------------------------------------------------------

def ess(weights: np.ndarray) -> float:
    """Effective sample size ``1 / sum(w_norm^2)`` (``N`` for uniform weights)."""
    w = np.asarray(weights, dtype=float).ravel()
    if w.size == 0:
        raise ValueError("weights must be non-empty")
    s = w.sum()
    if not np.isfinite(s) or s <= 0:
        raise ValueError("weights must have a positive finite sum")
    w = w / s
    return float(1.0 / np.sum(w**2))


def conditional_ess(log_weights: np.ndarray, log_increments: np.ndarray) -> float:
    """CESS = N (sum w exp(delta))^2 / sum w exp(2 delta).

    Unlike the ESS of reweighted particles, this equals N for a zero increment
    even when current weights are nonuniform. This permits CESS=.8N with
    independent ESS=.5N resampling without a temperature stall.
    """
    lw=np.asarray(log_weights,float).ravel();li=np.asarray(log_increments,float).ravel()
    if lw.shape!=li.shape:raise ValueError("weight/increment shape mismatch")
    lw=lw-logsumexp(lw)
    return float(np.exp(np.log(lw.size)+2*logsumexp(lw+li)-logsumexp(lw+2*li)))


def systematic_resample(weights: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Systematic resampling indices for normalised ``weights`` (low variance)."""
    w = np.asarray(weights, dtype=float).ravel()
    n = w.size
    s = w.sum()
    if not np.isfinite(s) or s <= 0:
        raise ValueError("weights must have a positive finite sum")
    w = w / s
    positions = (rng.random() + np.arange(n)) / n
    cdf = np.cumsum(w)
    cdf[-1] = 1.0
    return np.searchsorted(cdf, positions, side="left")


# ---------------------------------------------------------------------------
# Weighted summaries
# ---------------------------------------------------------------------------

def _norm_weights(weights: np.ndarray, n: int) -> np.ndarray:
    w = np.asarray(weights, dtype=float).ravel()
    if w.size != n:
        raise ValueError("weights length must match the number of samples")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("weights must be finite and non-negative")
    s = w.sum()
    if s <= 0:
        raise ValueError("weights must have a positive sum")
    return w / s


def weighted_mean(x: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted mean of ``x`` (leading axis = samples)."""
    x = np.asarray(x, dtype=float)
    w = _norm_weights(weights, x.shape[0])
    return np.tensordot(w, x, axes=(0, 0))


def weighted_variance(x: np.ndarray, weights: np.ndarray) -> float:
    """Posterior variance ``sum w (x - mean)^2`` for scalar ``x``."""
    x = np.asarray(x, dtype=float).ravel()
    w = _norm_weights(weights, x.size)
    m = float(w @ x)
    return float(w @ (x - m) ** 2)


def weighted_std(x: np.ndarray, weights: np.ndarray) -> float:
    """Posterior standard deviation ``sqrt(Var)``."""
    return float(np.sqrt(max(weighted_variance(x, weights), 0.0)))


def weighted_covariance(x: np.ndarray, weights: np.ndarray, mean: np.ndarray | None = None) -> np.ndarray:
    """Weighted covariance matrix of ``(N, d)`` samples."""
    x = np.atleast_2d(np.asarray(x, dtype=float))
    w = _norm_weights(weights, x.shape[0])
    m = weighted_mean(x, w) if mean is None else np.asarray(mean, dtype=float)
    d = x - m
    return (d.T * w) @ d


def weighted_quantiles(x: np.ndarray, weights: np.ndarray, probs: Sequence[float] | np.ndarray) -> np.ndarray:
    """Weighted quantiles via the cumulative-weight midpoint rule.

    ``x`` is a scalar sample array, ``probs`` probabilities in ``[0, 1]``.
    """
    x = np.asarray(x, dtype=float).ravel()
    w = _norm_weights(weights, x.size)
    order = np.argsort(x, kind="mergesort")
    xs, ws = x[order], w[order]
    cw = np.cumsum(ws)
    mid = (cw - 0.5 * ws) / cw[-1]
    p = np.asarray(probs, dtype=float).ravel()
    return np.interp(p, mid, xs)


def weighted_credible_interval(x: np.ndarray, weights: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    """Equal-tailed weighted credible interval ``[q_{(1-l)/2}, q_{(1+l)/2}]``."""
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie in (0, 1)")
    lo = 0.5 * (1.0 - level)
    q = weighted_quantiles(x, weights, [lo, 1.0 - lo])
    return float(q[0]), float(q[1])


def readout_statistics(
    q: np.ndarray,
    weights: np.ndarray,
    *,
    probs: Sequence[float] = (0.025, 0.16, 0.5, 0.84, 0.975),
    level: float = 0.9,
) -> dict:
    """Posterior summary of the scalar task readout ``Q``.

    Returns mean, variance, std, coefficient of variation, requested quantiles,
    an equal-tailed credible interval and the interval width.
    """
    q = np.asarray(q, dtype=float).ravel()
    w = _norm_weights(weights, q.size)
    m = float(w @ q)
    var = float(w @ (q - m) ** 2)
    std = float(np.sqrt(max(var, 0.0)))
    lo, hi = weighted_credible_interval(q, w, level=level)
    return {
        "mean": m,
        "variance": var,
        "std": std,
        "cv": float(std / abs(m)) if m != 0.0 else float("inf"),
        "quantiles": dict(zip([float(p) for p in probs], [float(v) for v in weighted_quantiles(q, w, probs)])),
        "credible_interval": (lo, hi),
        "credible_level": float(level),
        "interval_width": float(hi - lo),
        "ess": ess(w),
        "n": int(q.size),
    }


def _corr_from_cov(cov: np.ndarray) -> np.ndarray:
    sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    denom = np.outer(sd, sd)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.where(denom > 0, cov / denom, 0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def theta_statistics(nodes: np.ndarray, weights: np.ndarray) -> dict:
    """Weighted posterior mean/std, covariance and correlation of ``theta``."""
    nodes = np.atleast_2d(np.asarray(nodes, dtype=float))
    w = _norm_weights(weights, nodes.shape[0])
    mean = weighted_mean(nodes, w)
    cov = weighted_covariance(nodes, w, mean=mean)
    return {
        "mean": mean,
        "std": np.sqrt(np.clip(np.diag(cov), 0.0, None)),
        "covariance": cov,
        "correlation": _corr_from_cov(cov),
        "n": int(nodes.shape[0]),
        "ess": ess(w),
    }


def predictive_moments(
    pred: np.ndarray,
    weights: np.ndarray,
    noise_cov: np.ndarray | None = None,
) -> dict:
    """Posterior predictive mean and covariance of a complex external response.

    ``pred`` is ``(K, m)`` complex; the response is embedded in the ``2m`` real
    representation ``R = [Re(pred), Im(pred)]``.  With ``noise_cov`` (``2m x 2m``)
    the measurement noise is added to the parameter-induced covariance.
    """
    pred = np.atleast_2d(np.asarray(pred, dtype=complex))
    w = _norm_weights(weights, pred.shape[0])
    m = pred.shape[1]
    r = np.concatenate([pred.real, pred.imag], axis=1)  # (K, 2m)
    mean = weighted_mean(r, w)
    cov = weighted_covariance(r, w, mean=mean)
    if noise_cov is not None:
        noise_cov = np.asarray(noise_cov, dtype=float)
        if noise_cov.shape != (2 * m, 2 * m):
            raise ValueError(f"noise_cov must be {(2 * m, 2 * m)}")
        cov = cov + noise_cov
    return {
        "mean_real": mean,
        "covariance_real": cov,
        "mean_complex": mean[:m] + 1j * mean[m:],
        "covariance_complex": cov[:m, :m] + 1j * cov[:m, m:],
        "n_channels": int(m),
        "ess": ess(w),
    }


def branch_weights(branch_ids: np.ndarray, weights: np.ndarray, n_branches: int | None = None) -> dict:
    """Aggregate posterior weight per **externally supplied** branch label.

    ``branch_ids`` are integer component labels given by the caller (e.g. from a
    sign / clustering rule); this function performs no clustering itself.  Returns
    ``{"probabilities": {label: p}, "counts": {label: n}}`` with probabilities
    summing to 1 over labels present in ``n_branches`` (default: the observed
    labels).
    """
    ids = np.asarray(branch_ids).ravel()
    w = _norm_weights(weights, ids.size)
    labels = np.arange(n_branches) if n_branches is not None else np.unique(ids)
    probs, counts = {}, {}
    for lab in labels:
        mask = ids == lab
        probs[int(lab)] = float(w[mask].sum())
        counts[int(lab)] = int(mask.sum())
    return {"probabilities": probs, "counts": counts}


def posterior_summary(
    nodes: np.ndarray,
    log_likelihood: np.ndarray,
    *,
    prior_weights: np.ndarray | None = None,
    q: np.ndarray | None = None,
    pred: np.ndarray | None = None,
    branch_ids: np.ndarray | None = None,
    noise_cov: np.ndarray | None = None,
    strict: bool = True,
    invalid_mass_tol: float = 1e-12,
    quantile_probs: Sequence[float] = (0.025, 0.16, 0.5, 0.84, 0.975),
    level: float = 0.9,
) -> dict:
    """One-call reference posterior summary over a precomputed tensor grid.

    ``nodes`` is ``(K, d)``; ``log_likelihood`` is ``(K,)``; optional ``q`` is
    the scalar task readout, ``pred`` the complex external response, and
    ``branch_ids`` externally supplied branch labels.  Returns normalised
    weights, evidence, theta/Q/predictive/branch summaries and rejection info.
    """
    nodes = np.atleast_2d(np.asarray(nodes, dtype=float))
    log_w = posterior_log_weights(log_likelihood, prior_weights)
    norm = normalize_log_weights(log_w, prior_weights, strict=strict, invalid_mass_tol=invalid_mass_tol)
    w = norm["weights"]
    out = {
        "weights": w,
        "log_weights": norm["log_weights"],
        "log_evidence": norm["log_evidence"],
        "ess": norm["ess"],
        "rejection": norm["info"],
        "theta": theta_statistics(nodes, w),
    }
    if q is not None:
        out["q"] = readout_statistics(q, w, probs=quantile_probs, level=level)
    if pred is not None:
        out["predictive"] = predictive_moments(pred, w, noise_cov=noise_cov)
    if branch_ids is not None:
        out["branch"] = branch_weights(branch_ids, w)
    return out


# ---------------------------------------------------------------------------
# Population risk, prevalence and paired bootstrap
# ---------------------------------------------------------------------------

def population_risk(variances: np.ndarray, query_weights: np.ndarray | None = None) -> float:
    """Population task risk ``R_Q = E_Y Var(Q | Y)`` from per-query variances."""
    v = np.asarray(variances, dtype=float).ravel()
    if v.size == 0:
        raise ValueError("variances must be non-empty")
    if query_weights is None:
        return float(np.mean(v))
    w = _norm_weights(query_weights, v.size)
    return float(w @ v)


def probability_sd_above(variances: np.ndarray, tau: float, query_weights: np.ndarray | None = None) -> float:
    """``P_Y[SD(Q | Y) > tau]`` from per-query variances (Gate A1 statistic)."""
    v = np.asarray(variances, dtype=float).ravel()
    sd = np.sqrt(np.clip(v, 0.0, None))
    hit = sd > float(tau)
    if query_weights is None:
        return float(np.mean(hit))
    w = _norm_weights(query_weights, v.size)
    return float(w @ hit.astype(float))


def paired_bootstrap_risk(
    variance_by_method: Mapping[str, np.ndarray],
    *,
    replicates: int = 5000,
    seed: int = 0,
    ci: float = 0.95,
    tau: float | None = None,
    baseline: str | None = None,
) -> dict:
    """Query-level **paired** bootstrap of population risk across methods.

    ``variance_by_method`` maps a method name to a length-``n_queries`` array of
    per-query ``Var(Q | Y)`` (paired: index ``j`` is the same query for every
    method).  Query indices are resampled with replacement once per replicate and
    reused for every method, so risk differences/reductions stay paired.

    Returns per-method risk with CI, optional prevalence ``P[SD > tau]`` with CI,
    and -- relative to ``baseline`` (default: first key) -- the paired risk
    reduction ``1 - R_m/R_base`` and difference ``R_base - R_m`` with CIs.
    """
    names = list(variance_by_method.keys())
    if not names:
        raise ValueError("variance_by_method must be non-empty")
    arrays = {k: np.asarray(variance_by_method[k], dtype=float).ravel() for k in names}
    n = arrays[names[0]].size
    for k, a in arrays.items():
        if a.size != n:
            raise ValueError(f"method {k!r} has {a.size} queries, expected {n}")
    base = baseline if baseline is not None else names[0]
    if base not in arrays:
        raise ValueError(f"baseline {base!r} not among methods {names}")

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(int(replicates), n))
    alpha = 0.5 * (1.0 - ci)

    def _ci(samples: np.ndarray) -> tuple[float, float]:
        return (float(np.quantile(samples, alpha)), float(np.quantile(samples, 1.0 - alpha)))

    risks = {k: float(np.mean(a)) for k, a in arrays.items()}
    rep_risk = {k: arrays[k][idx].mean(axis=1) for k in names}
    out = {
        "n_queries": int(n),
        "replicates": int(replicates),
        "ci": float(ci),
        "baseline": base,
        "paired": True,
        "seed": int(seed),
        "risk": risks,
        "risk_ci": {k: _ci(rep_risk[k]) for k in names},
    }

    if tau is not None:
        hit = {k: (np.sqrt(np.clip(arrays[k], 0.0, None)) > float(tau)).astype(float) for k in names}
        out["tau"] = float(tau)
        out["prevalence"] = {k: float(np.mean(hit[k])) for k in names}
        out["prevalence_ci"] = {k: _ci(hit[k][idx].mean(axis=1)) for k in names}

    reduction, reduction_ci, difference, difference_ci = {}, {}, {}, {}
    for k in names:
        if k == base:
            continue
        r = 1.0 - rep_risk[k] / rep_risk[base]
        dd = rep_risk[base] - rep_risk[k]
        reduction[k] = float(np.mean(r))
        reduction_ci[k] = _ci(r)
        difference[k] = float(np.mean(dd))
        difference_ci[k] = _ci(dd)
    out["risk_reduction"] = reduction
    out["risk_reduction_ci"] = reduction_ci
    out["risk_difference"] = difference
    out["risk_difference_ci"] = difference_ci
    return out


# ---------------------------------------------------------------------------
# Reference-grid comparison
# ---------------------------------------------------------------------------

def _as_1d(a) -> np.ndarray:
    return np.asarray(a, dtype=float).ravel()


def compare_posterior_to_reference(
    reference: Mapping[str, object],
    candidate: Mapping[str, object],
    tolerances: Mapping[str, float] | None = None,
    *,
    config_path: str | None = None,
) -> dict:
    """Compare a candidate posterior summary to a reference using frozen tolerances.

    Both mappings may contain ``mean``, ``variance``, ``quantiles`` and
    ``branch_weights`` (scalars or arrays).  Passes when the max absolute mean
    difference ``<= reference_mean_tol``, every variance satisfies the relative
    *or* absolute tolerance, the max quantile difference ``<=
    reference_quantile_tol`` and the max branch-weight difference ``<=
    reference_branch_tol``.
    """
    tol = dict(load_reference_tolerances(config_path)) if tolerances is None else dict(tolerances)
    checks: dict = {}
    ok = True

    if "mean" in reference and "mean" in candidate:
        val = float(np.max(np.abs(_as_1d(reference["mean"]) - _as_1d(candidate["mean"]))))
        passed = val <= tol["reference_mean_tol"]
        checks["mean"] = {"max_abs_diff": val, "tol": tol["reference_mean_tol"], "passed": bool(passed)}
        ok = ok and bool(passed)

    if "variance" in reference and "variance" in candidate:
        ref = _as_1d(reference["variance"])
        cand = _as_1d(candidate["variance"])
        rel = np.abs(cand - ref) / np.maximum(np.abs(ref), 1e-300)
        absd = np.abs(cand - ref)
        per = (rel <= tol["reference_variance_rel_tol"]) | (absd <= tol["reference_variance_abs_tol"])
        passed = bool(np.all(per))
        checks["variance"] = {
            "max_rel_diff": float(np.max(rel)) if rel.size else 0.0,
            "max_abs_diff": float(np.max(absd)) if absd.size else 0.0,
            "rel_tol": tol["reference_variance_rel_tol"],
            "abs_tol": tol["reference_variance_abs_tol"],
            "passed": passed,
        }
        ok = ok and passed

    if "quantiles" in reference and "quantiles" in candidate:
        val = float(np.max(np.abs(_as_1d(reference["quantiles"]) - _as_1d(candidate["quantiles"]))))
        passed = val <= tol["reference_quantile_tol"]
        checks["quantiles"] = {"max_abs_diff": val, "tol": tol["reference_quantile_tol"], "passed": bool(passed)}
        ok = ok and bool(passed)

    if "branch_weights" in reference and "branch_weights" in candidate:
        val = float(np.max(np.abs(_as_1d(reference["branch_weights"]) - _as_1d(candidate["branch_weights"]))))
        passed = val <= tol["reference_branch_tol"]
        checks["branch_weights"] = {"max_abs_diff": val, "tol": tol["reference_branch_tol"], "passed": bool(passed)}
        ok = ok and bool(passed)

    return {"checks": checks, "passed": bool(ok), "tolerances": tol}


def check_convergence_metrics(
    metrics: Mapping[str, float],
    tolerances: Mapping[str, float] | None = None,
    *,
    config_path: str | None = None,
) -> dict:
    """Upper-bound checks for mesh/linear convergence metrics from discovery.json.

    Recognised metric keys map to tolerances: ``mesh_q`` -> ``mesh_q_tol``,
    ``whitened_rms`` -> ``mesh_whitened_rms_tol``, ``whitened_max`` ->
    ``mesh_whitened_max_tol``, ``linear_residual`` -> ``linear_residual_tol``.
    """
    tol = dict(load_reference_tolerances(config_path)) if tolerances is None else dict(tolerances)
    mapping = {
        "mesh_q": "mesh_q_tol",
        "whitened_rms": "mesh_whitened_rms_tol",
        "whitened_max": "mesh_whitened_max_tol",
        "linear_residual": "linear_residual_tol",
    }
    checks, ok = {}, True
    for name, value in metrics.items():
        if name not in mapping:
            raise KeyError(f"unknown convergence metric {name!r}; expected one of {sorted(mapping)}")
        passed = float(value) <= tol[mapping[name]]
        checks[name] = {"value": float(value), "tol": tol[mapping[name]], "passed": bool(passed)}
        ok = ok and bool(passed)
    return {"checks": checks, "passed": bool(ok), "tolerances": tol}


# ---------------------------------------------------------------------------
# Adaptive tempering SMC (generic sampler; no physical claim)
# ---------------------------------------------------------------------------

def _validate_bounds(bounds) -> np.ndarray:
    b = np.asarray(bounds, dtype=float)
    if b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("bounds must be shape (d, 2)")
    if np.any(b[:, 1] <= b[:, 0]):
        raise ValueError("every box bound must satisfy b_i > a_i")
    return b


def adaptive_tempered_smc(
    log_likelihood: Callable[[np.ndarray], np.ndarray],
    bounds: Sequence | np.ndarray,
    *,
    n_particles: int = 1024,
    seed: int = 0,
    ess_fraction: float = 0.8,
    resample_threshold: float = 0.5,
    pcn_steps: int = 5,
    pcn_beta: float = 0.2,
    max_steps: int = 512,
    min_lambda_increment: float = 1e-9,
    readout: Callable[[np.ndarray], np.ndarray] | None = None,
    branch_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> dict:
    """Adaptive-tempering SMC on an independent standard-Gaussian base ``u``.

    The bounded uniform prior is realised exactly through the CDF transform
    ``theta(u) = a + (b - a) Phi(u)``; because the prior density is constant and
    the Jacobian is the Gaussian density, the SMC target is exactly
    ``propto N(u; 0, I) * L(theta(u))``.  Temperatures ``lambda`` on ``log L``
    are chosen adaptively so the (conditional) ESS hits ``ess_fraction *
    n_particles``; particles are systematically resampled when the ESS drops
    below ``resample_threshold * n_particles`` and rejuvenated with ``pcn_steps``
    preconditioned Crank-Nicolson sweeps.

    ``log_likelihood(theta)`` must accept an ``(N, d)`` array and return ``(N,)``;
    any non-finite value is rejected (``ValueError``).  Returns particles,
    weights, per-step diagnostics, the SMC log-evidence estimate and optional
    readout/branch summaries.  **This routine is a numerical tool only and makes
    no physical claim.**
    """
    b = _validate_bounds(bounds)
    d = b.shape[0]
    if n_particles < 2:
        raise ValueError("n_particles must be >= 2")
    if not 0.0 < ess_fraction <= 1.0:
        raise ValueError("ess_fraction must lie in (0, 1]")
    if not 0.0 < resample_threshold <= 1.0:
        raise ValueError("resample_threshold must lie in (0, 1]")
    rng = np.random.default_rng(seed)
    n = int(n_particles)
    target_ess = ess_fraction * n

    u = rng.standard_normal((n, d))
    theta = uniform_cdf_transform(u, b)
    log_l = np.asarray(log_likelihood(theta), dtype=float).ravel()
    if log_l.size != n:
        raise ValueError("log_likelihood must return one value per particle")
    if not np.isfinite(log_l).all():
        bad = np.nonzero(~np.isfinite(log_l))[0]
        raise ValueError(
            f"SMC received {bad.size} non-finite log-likelihood value(s) at {bad[:8].tolist()}"
            f"{'...' if bad.size > 8 else ''}; refusing to temper a failed likelihood"
        )
    n_evals = n

    log_w = np.full(n, -np.log(n))
    lam = 0.0
    log_evidence = 0.0
    diagnostics: list[dict] = []

    def _ess_for_increment(log_increments: np.ndarray) -> float:
        return conditional_ess(log_w, log_increments)

    ancestors=np.arange(n)
    step = 0
    while lam < 1.0 - 1e-15 and step < max_steps:
        remaining = 1.0 - lam
        if _ess_for_increment(remaining * log_l) >= target_ess:
            dlam = remaining
        else:
            lo, hi = 0.0, remaining
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if _ess_for_increment(mid * log_l) >= target_ess:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-12:
                    break
            dlam = lo
        if dlam < min_lambda_increment:
            break

        log_w_prev = log_w.copy()
        log_inc = dlam * log_l
        log_unnorm = log_w + log_inc
        log_inc_evidence = float(logsumexp(log_unnorm))  # prior weights sum to 1
        log_evidence += log_inc_evidence
        log_w = log_unnorm - log_inc_evidence
        lam += dlam
        w = np.exp(log_w)
        current_ess = ess(w)
        # Conditional ESS of the reweighted particles, computed from the actual
        # pre-increment weights; different from post-increment ESS.
        cond_ess = conditional_ess(log_w_prev, log_inc)

        resampled = False
        if current_ess < resample_threshold * n:
            idx = systematic_resample(w, rng)
            ancestors=ancestors[idx]
            u = u[idx]
            theta = theta[idx]
            log_l = log_l[idx]
            log_w = np.full(n, -np.log(n))
            resampled = True
            ess_after = float(n)
        else:
            ess_after = current_ess

        accepts = 0
        if pcn_steps > 0:
            scale = np.sqrt(1.0 - pcn_beta**2)
            for _ in range(int(pcn_steps)):
                prop = scale * u + pcn_beta * rng.standard_normal((n, d))
                theta_prop = uniform_cdf_transform(prop, b)
                log_l_prop = np.asarray(log_likelihood(theta_prop), dtype=float).ravel()
                n_evals += n
                if not np.isfinite(log_l_prop).all():
                    bad = np.nonzero(~np.isfinite(log_l_prop))[0]
                    raise ValueError(
                        f"SMC rejuvenation produced non-finite log-likelihood at {bad[:8].tolist()}"
                        f"{'...' if bad.size > 8 else ''}"
                    )
                accept = np.log(rng.random(n)) < lam * (log_l_prop - log_l)
                u = np.where(accept[:, None], prop, u)
                theta = np.where(accept[:, None], theta_prop, theta)
                log_l = np.where(accept, log_l_prop, log_l)
                accepts += int(accept.sum())

        diagnostics.append({
            "step": step,
            "lambda": float(lam),
            "delta_lambda": float(dlam),
            "ess": float(current_ess),
            "conditional_ess": float(cond_ess),
            "log_evidence_increment": float(log_inc_evidence),
            "resampled": bool(resampled),
            "unique_ancestors": int(np.unique(ancestors).size),
            "ess_after_resample": float(ess_after),
            "pcn_steps": int(pcn_steps),
            "pcn_acceptance": float(accepts / (int(pcn_steps) * n)) if pcn_steps > 0 else float("nan"),
            "log_likelihood_evals": int(n_evals),
        })
        step += 1

    weights = np.exp(log_w)
    weights = weights / weights.sum()
    out = {
        "u": u,
        "theta": theta,
        "log_weights": log_w,
        "weights": weights,
        "log_likelihood": log_l,
        "lambda": float(lam),
        "log_evidence": float(log_evidence),
        "ess": float(ess(weights)),
        "converged": bool(lam >= 1.0 - 1e-12),
        "step_count": int(step),
        "n_particles": int(n),
        "n_likelihood_evals": int(n_evals),
        "diagnostics": diagnostics,
    }
    if readout is not None:
        q = np.asarray(readout(theta), dtype=float).ravel()
        out["readout_statistics"] = readout_statistics(q, weights)
    if branch_fn is not None:
        out["branch_weights"] = branch_weights(branch_fn(theta), weights)
    return out


class TemperedSMC:
    """Thin reusable wrapper around :func:`adaptive_tempered_smc`.

    Stores the caller-supplied log-likelihood and box bounds; ``sample`` runs one
    tempered-SMC pass.  No physical claim is made anywhere in this class.
    """

    def __init__(
        self,
        log_likelihood: Callable[[np.ndarray], np.ndarray],
        bounds: Sequence | np.ndarray,
    ) -> None:
        self.log_likelihood = log_likelihood
        self.bounds = _validate_bounds(bounds)

    def sample(self, *, seed: int = 0, n_particles: int = 1024, readout=None, branch_fn=None, **kwargs) -> dict:
        return adaptive_tempered_smc(
            self.log_likelihood,
            self.bounds,
            n_particles=n_particles,
            seed=seed,
            readout=readout,
            branch_fn=branch_fn,
            **kwargs,
        )
