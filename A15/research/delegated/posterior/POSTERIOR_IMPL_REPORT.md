# A15 posterior utilities — bounded implementation worker report

Scope: implement **only** `src/posterior.py` and `tests/test_posterior.py`. No
physics run, no gate/protocol/threshold change, no Flow, no torch, no training.
This is a mechanical numerics delivery; it makes no physical claim.

## What was delivered

- `src/posterior.py` — numpy/scipy-only toolkit:
  - tensor Gauss–Legendre rule for an independent bounded uniform 3D prior
    (`gauss_legendre_box_nodes`): positive weights normalised to `sum(w)==1`, plus
    raw physical quadrature weights.
  - exact uniform↔normal CDF transform (`uniform_cdf_transform`,
    `uniform_inverse_cdf`) for the SMC base and for prior-constrained decoding.
  - real/imag independent Gaussian log-likelihood with sigma per **real**
    component (`complex_gaussian_log_likelihood`; scalar / `(m,)` / `(2m,)`).
  - `posterior_log_weights` + `normalize_log_weights` (logsumexp; strict rejection
    of non-finite/missing nodes, never silent renormalisation).
  - weighted Q mean/var/std/CV/quantiles/credible interval, theta mean/std/
    covariance/correlation, posterior predictive mean/covariance (real 2m
    embedding, optional noise), externally-supplied branch weights.
  - population task risk, `P[SD(Q|Y)>tau]`, and query-level **paired** bootstrap
    of risk / risk-reduction / risk-difference across methods.
  - reference-grid comparison tolerances loaded from `configs/discovery.json`
    (mean/variance/quantile/branch) plus mesh/linear convergence upper-bound checks.
  - adaptive-tempering SMC with independent-Gaussian base, exact CDF transform,
    pCN rejuvenation, correctly-computed conditional ESS, ESS-based systematic
    resampling and complete per-step diagnostic records; log-likelihood and Q
    readout are caller-supplied.
- `tests/test_posterior.py` — 34 `unittest` tests (analytic Gaussian,
  zero-information uniform, bimodal 1D embedded in 3D, rejection, risk/bootstrap,
  reference tolerances, SMC).

## Environment / how to run

Path resolution: from `Flow Matching/A15` the valid interpreter is
`../../Gaussian/.venv_nn/bin/python` (the task's `../../../` variant does **not**
resolve). Python 3.13.13, NumPy 2.5.3, SciPy 1.18.1. `pytest` is **not**
installed in that venv, so tests are `unittest`-based (pytest-collectible too).

```text
cd "Flow Matching/A15"
../../Gaussian/.venv_nn/bin/python -m unittest -v tests.test_posterior
```

Result: `Ran 34 tests ... OK`. Full log: `research/delegated/posterior/test_run.log`.

## Bugs found and fixed during testing (implementation, not physics)

1. `normalize_log_weights` was adding the prior a second time even though
   `posterior_log_weights` already returns `log prior + logL` (weights ∝ prior²;
   corner nodes over-penalised, posterior variance biased to 0.562 vs 0.75).
   Fixed: `normalize_log_weights` now normalises the supplied log-posterior and uses
   `prior_weights` only to flag missing nodes and account rejected prior mass.
2. Adaptive tempering stalled at the first temperature when
   `resample_threshold == ess_fraction` (target ESS reached exactly, no resample,
   next increment becomes ~0). Fixed: default `resample_threshold = 2/3 >
   ess_fraction = 1/2`, with an explicit guard that raises if
   `resample_threshold <= ess_fraction`.

## Independent cross-checks (beyond unit tests)

See `research/delegated/posterior/DIAGNOSTICS.json`.

- SMC log-evidence for an analytic 3D Gaussian (`-5.7575`) agrees with an
  independent 161³ Gauss–Legendre grid marginal likelihood (`-5.7679`) to
  `|Δ| = 0.0104`; uniform prior integrates to `logZ == 0`.
- Rejection: non-finite node raises under `strict=True`; under `strict=False` it
  raises if the node's prior mass exceeds tolerance, and drops it (weight 0,
  weights still sum to 1) only when negligible.
- Determinism: identical `seed` reproduces particles and weights exactly.
- Dimensionality: SMC works for `d=1,2,3` (d=1 mean 0.211 vs target 0.2; d=2 mean
  `[0.097, -0.184]` vs target `[0.1, -0.2]`).
- Predictive covariance is symmetric and positive definite after adding noise.
- `torch` is never imported.

## Uncertainties / limitations (do not over-read)

- "Conditional ESS" here is the ESS of the reweighted particles (equal to the
  post-reweight ESS by definition); because the default resamples at every
  tempering step, this equals the common formulation. It is a numerical gauge, not
  a physical quantity.
- The weighted credible interval is nominal/equal-tailed (midpoint cumulative
  weight interpolation); no finite-sample or HPD correction.
- Reference-grid *smooth* functionals (moments) are near-exact (~1e-12); the GL
  rule is only ~1e-2 accurate for discontinuous functionals (indicator/CDF
  evaluations), as expected for Gauss–Legendre.
- SMC log-evidence is a stochastic estimate and is offered only as a numerical
  diagnostic; it carries no physical meaning.
- This module computes posteriors/risks from caller-supplied likelihoods and
  readouts; it cannot and does not establish whether any ambiguity is real.

## Scope compliance

Only `src/posterior.py` and `tests/test_posterior.py` were created/edited.
`configs/discovery.json` is read-only input (tolerances); Gate A/B, thresholds,
protocols, A14 files and Flow state were not touched.
