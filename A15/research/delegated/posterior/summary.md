Done. Delivered only the two assigned files plus delegated artifacts.

**Findings (bounded mechanical work; no physics, no gate change)**

1. Implemented [src/posterior.py](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/src/posterior.py>): tensor Gauss–Legendre nodes with positive weights normalised to `sum(w)==1`; exact uniform↔normal CDF transform; real/imag independent Gaussian likelihood with sigma per real component (scalar/`(m,)`/`(2m,)`); logsumexp normalisation that rejects non-finite/missing nodes instead of renormalising; weighted Q mean/var/std/quantiles/CI, theta covariance/correlation, predictive mean/covariance, externally-supplied branch weights; population risk, `P[SD>τ]`, and query-level paired bootstrap; discovery.json reference/convergence tolerances; adaptive-tempering SMC (Gaussian base via CDF transform, pCN rejuvenation, conditional ESS, ESS resampling, full per-step diagnostics). NumPy/SciPy only — no torch, no training.

2. Implemented [tests/test_posterior.py](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/tests/test_posterior.py>) with 34 `unittest` tests covering analytic Gaussian, zero-information uniform, and bimodal 1D embedded in 3D, plus rejection, risk/bootstrap, tolerances and SMC.

3. Interpreter resolved as instructed: from A15 the valid path is `../../Gaussian/.venv_nn/bin/python` (the `../../../` variant does not resolve). That venv has Python 3.13.13 / NumPy 2.5.3 / SciPy 1.18.1 and no pytest, so tests are `unittest`-based (also pytest-collectible).

4. Test result: `Ran 34 tests … OK` (verbose log at [test_run.log](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/research/delegated/posterior/test_run.log>)). Both `python -m unittest` and direct `python tests/test_posterior.py` run clean; `py_compile` clean.

5. Two real implementation bugs were found by the tests and fixed: prior double-counted in `normalize_log_weights` (weights ∝ prior², biased variance 0.562 vs 0.75), and an adaptive-tempering stall when `resample_threshold == ess_fraction` (now default `2/3 > 1/2` with an explicit guard that raises otherwise).

6. Independent cross-checks (not just unit tests), in [DIAGNOSTICS.json](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/research/delegated/posterior/DIAGNOSTICS.json>): SMC log-evidence `-5.7575` vs an independent 161³ Gauss–Legendre grid marginal `-5.7679` (|Δ|=0.0104); uniform prior integrates to `logZ==0`; same-seed runs reproduce particles/weights exactly; SMC works for d=1,2,3; predictive covariance symmetric PSD after noise.

7. Full report at [POSTERIOR_IMPL_REPORT.md](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/research/delegated/posterior/POSTERIOR_IMPL_REPORT.md>). Scope respected: only the two assigned files were created/edited; `configs/discovery.json` is read-only input; gates, thresholds, protocol, A14 files and Flow state untouched. Other `src/*` files that appeared concurrently belong to other agents and were left alone.

**Uncertainties**
- "Conditional ESS" here is the ESS of the reweighted particles (equals post-reweight ESS by definition); because the default resamples each tempering step it matches the common formulation. It is a numerical gauge, not physical.
- The weighted credible interval is nominal/equal-tailed (midpoint cumulative-weight interpolation); no finite-sample or HPD correction.
- Gauss–Legendre is near-exact (~1e-12) for smooth functionals but only ~1e-2 accurate for discontinuous ones (indicator/CDF evaluations), so reference-grid CDF comparisons need that tolerance.
- SMC log-evidence is a stochastic diagnostic with no physical meaning; this module cannot establish whether any ambiguity is physically real — that remains the scientific owner's call.

