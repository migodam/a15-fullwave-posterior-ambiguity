"""Mechanical tests for src/posterior.py.

Covers: tensor Gauss-Legendre prior nodes/weights; the exact uniform CDF
transform; real/imag Gaussian likelihood; logsumexp normalisation with rejection
of failed nodes; weighted Q/theta/predictive summaries; externally-supplied branch
weights; population risk and query-level paired bootstrap; discovery-config
reference tolerances; and adaptive-tempering SMC (analytic Gaussian, zero
information, bimodal 1D embedded in 3D).  No physics, no torch, no training.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import posterior as P  # noqa: E402


BOX = [[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]]


class GaussLegendrePriorTest(unittest.TestCase):
    def test_nodes_positive_normalised_weights(self):
        bounds = [[0.0, 1.0], [0.0, 2.0], [-1.0, 1.0]]
        rule = P.gauss_legendre_box_nodes(bounds, [9, 7, 5])
        n = rule["nodes"]
        w = rule["weights"]
        self.assertEqual(n.shape, (9 * 7 * 5, 3))
        self.assertTrue(np.all(w > 0))
        self.assertAlmostEqual(float(w.sum()), 1.0, places=12)
        self.assertAlmostEqual(float(rule["raw_weights"].sum()), rule["volume"], places=12)
        for i, (a, b) in enumerate(bounds):
            self.assertTrue(np.all(n[:, i] >= a) and np.all(n[:, i] <= b))

    def test_exact_integration_of_moments(self):
        bounds = [[0.0, 1.0], [0.0, 2.0], [-1.0, 1.0]]
        rule = P.gauss_legendre_box_nodes(bounds, [9, 7, 5])
        n, w = rule["nodes"], rule["weights"]
        self.assertAlmostEqual(float(w.sum()), 1.0, places=12)
        self.assertAlmostEqual(float(w @ n[:, 0]), 0.5, places=12)
        self.assertAlmostEqual(float(w @ n[:, 0] ** 2), 1.0 / 3.0, places=12)
        self.assertAlmostEqual(float(w @ n[:, 1]), 1.0, places=12)
        self.assertAlmostEqual(float(w @ n[:, 1] ** 2), 4.0 / 3.0, places=12)
        self.assertAlmostEqual(float(w @ (n[:, 0] * n[:, 1])), 0.5, places=12)


class CdfTransformTest(unittest.TestCase):
    def test_roundtrip_inside_box(self):
        rng = np.random.default_rng(3)
        u = rng.standard_normal((200, 3))
        theta = P.uniform_cdf_transform(u, BOX)
        self.assertTrue(np.all(theta > np.array([b[0] for b in BOX])))
        self.assertTrue(np.all(theta < np.array([b[1] for b in BOX])))
        u2 = P.uniform_inverse_cdf(theta, BOX)
        np.testing.assert_allclose(u2, u, atol=1e-9)

    def test_endpoints_map_finite(self):
        u = np.array([[40.0, -40.0, 0.0]])
        theta = P.uniform_cdf_transform(u, BOX)
        self.assertTrue(np.isfinite(theta).all())
        self.assertTrue(np.isfinite(P.uniform_inverse_cdf(theta, BOX)).all())


class LikelihoodTest(unittest.TestCase):
    def test_real_imag_gaussian_matches_manual(self):
        obs = np.array([0.7 - 0.2j, -0.4 + 0.9j])
        pred = np.array([[0.6 - 0.1j, -0.4 + 0.8j], [0.8 - 0.3j, -0.3 + 1.0j]])
        sigma = 0.15
        ll = P.complex_gaussian_log_likelihood(pred, obs, sigma)
        manual = np.zeros(2)
        for k in range(2):
            r = pred[k] - obs
            manual[k] = -0.5 * np.sum(r.real**2 + r.imag**2) / sigma**2
            manual[k] -= 0.5 * 2 * 2 * np.log(2 * np.pi * sigma**2)
        np.testing.assert_allclose(ll, manual, atol=1e-12)

    def test_sigma_shapes_and_zero_at_observation(self):
        obs = np.array([1.0 + 2.0j, -1.0 - 1.0j])
        sigma = np.array([0.1, 0.2, 0.3, 0.4])  # per real component (Re, Re, Im, Im)
        pred = obs[None, :]
        ll = P.complex_gaussian_log_likelihood(pred, obs, sigma)
        expected = -0.5 * (
            np.log(2 * np.pi * 0.1**2)
            + np.log(2 * np.pi * 0.2**2)
            + np.log(2 * np.pi * 0.3**2)
            + np.log(2 * np.pi * 0.4**2)
        )
        self.assertAlmostEqual(float(ll[0]), float(expected), places=12)
        # no-const likelihood is 0 at the observation
        ll0 = P.complex_gaussian_log_likelihood(pred, obs, sigma, include_const=False)
        self.assertAlmostEqual(float(ll0[0]), 0.0, places=12)


class NormalizeRejectTest(unittest.TestCase):
    def test_rejects_nonfinite_strict(self):
        lw = np.zeros(50)
        lw[3] = np.nan
        with self.assertRaises(ValueError):
            P.normalize_log_weights(lw)

    def test_rejects_mass_above_tol_when_not_strict(self):
        lw = np.zeros(50)
        lw[3] = np.nan
        with self.assertRaises(ValueError):
            P.normalize_log_weights(lw, strict=False)  # mass 0.02 > 1e-12

    def test_drops_only_when_negligible(self):
        lw = np.zeros(1000)
        lw[7] = np.nan
        out = P.normalize_log_weights(lw, strict=False, invalid_mass_tol=1e-2)
        self.assertEqual(out["info"]["n_rejected"], 1)
        self.assertEqual(int(out["info"]["rejected_indices"][0]), 7)
        self.assertAlmostEqual(out["info"]["rejected_prior_mass"], 0.001, places=12)
        self.assertEqual(out["weights"][7], 0.0)
        self.assertAlmostEqual(float(out["weights"].sum()), 1.0, places=12)

    def test_zero_prior_weight_is_missing_node(self):
        pw = np.full(10, 0.1)
        pw[2] = 0.0
        with self.assertRaises(ValueError):
            P.normalize_log_weights(np.zeros(10), pw)

    def test_ess_uniform_and_concentrated(self):
        self.assertAlmostEqual(P.ess(np.full(20, 1 / 20)), 20.0, places=9)
        concentrated = np.zeros(20)
        concentrated[0] = 1.0
        self.assertAlmostEqual(P.ess(concentrated), 1.0, places=9)

    def test_conditional_ess_definition(self):
        rng = np.random.default_rng(5)
        log_w = np.log(rng.dirichlet(np.ones(30)))
        log_inc = rng.normal(scale=0.5, size=30)
        unnorm = log_w + log_inc
        w = np.exp(unnorm - P.logsumexp(unnorm))
        self.assertAlmostEqual(P.conditional_ess(log_w, log_inc), len(log_w)*(np.exp(log_w)@np.exp(log_inc))**2/(np.exp(log_w)@np.exp(2*log_inc)), places=9)
        self.assertAlmostEqual(P.conditional_ess(log_w, np.zeros_like(log_inc)), len(log_w), places=9)

    def test_systematic_resample_frequency(self):
        w = np.array([0.7, 0.2, 0.1])
        idx = P.systematic_resample(w, np.random.default_rng(0))
        self.assertEqual(idx.size, 3)
        self.assertTrue(np.all((idx >= 0) & (idx < 3)))


class GridPosteriorTest(unittest.TestCase):
    def _grid(self, order=(81, 9, 9)):
        rule = P.gauss_legendre_box_nodes(BOX, list(order))
        return rule["nodes"], rule["weights"]

    def test_analytic_gaussian_section(self):
        nodes, pw = self._grid()
        mu, s = 0.6, 0.25
        lik = -0.5 * ((nodes[:, 0] - mu) / s) ** 2
        norm = P.normalize_log_weights(P.posterior_log_weights(lik, pw), pw)
        stats = P.readout_statistics(nodes[:, 0], norm["weights"])
        self.assertLess(abs(stats["mean"] - mu), 2e-3)
        self.assertLess(abs(stats["variance"] - s**2) / s**2, 5e-3)

    def test_zero_information_is_prior(self):
        nodes, pw = self._grid()
        norm = P.normalize_log_weights(P.posterior_log_weights(np.zeros(nodes.shape[0]), pw), pw)
        np.testing.assert_allclose(norm["weights"], pw, atol=1e-15)
        stats = P.readout_statistics(nodes[:, 0], norm["weights"])
        self.assertAlmostEqual(stats["mean"], 0.0, places=10)
        self.assertAlmostEqual(stats["variance"], 0.75, places=9)  # (3)^2 / 12
        # The prior weights are Gauss-Legendre weights (not uniform over nodes), so
        # the ESS is the (finite) ESS of the prior rule, never larger than K.
        self.assertAlmostEqual(norm["ess"], P.ess(pw), places=9)
        self.assertLessEqual(norm["ess"], nodes.shape[0] + 1e-9)

    def test_weighted_quantiles_discrete_exact(self):
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        w = np.ones(5)
        q = P.weighted_quantiles(x, w, [0.1, 0.5, 0.9])
        np.testing.assert_allclose(q, [0.0, 2.0, 4.0], atol=1e-12)
        lo, hi = P.weighted_credible_interval(x, w, level=0.8)
        self.assertAlmostEqual(lo, 0.0, places=12)
        self.assertAlmostEqual(hi, 4.0, places=12)
        # concentrated mass: 90% at x=0 pulls the (midpoint-rule) median to 0.1
        self.assertAlmostEqual(float(P.weighted_quantiles([0.0, 1.0], [0.9, 0.1], [0.5])[0]), 0.1, places=12)

    def test_prior_cdf_sanity(self):
        nodes, pw = self._grid()
        # Indicator quadrature is only ~1% accurate on Gauss-Legendre nodes.
        self.assertLess(abs(float(pw[nodes[:, 0] <= -1.2].sum()) - 0.1), 0.02)
        self.assertLess(abs(float(pw[nodes[:, 0] <= 0.6].sum()) - 0.7), 0.02)

    def test_theta_covariance_uniform(self):
        nodes, pw = self._grid(order=(41, 9, 9))
        w = P.normalize_log_weights(P.posterior_log_weights(np.zeros(nodes.shape[0]), pw), pw)["weights"]
        ts = P.theta_statistics(nodes, w)
        np.testing.assert_allclose(ts["mean"], np.zeros(3), atol=1e-10)
        np.testing.assert_allclose(ts["covariance"], np.eye(3) * 0.75, atol=1e-9)
        np.testing.assert_allclose(ts["correlation"], np.eye(3), atol=1e-9)

    def test_predictive_moments(self):
        nodes, pw = self._grid(order=(41, 9, 9))
        w = P.normalize_log_weights(P.posterior_log_weights(np.zeros(nodes.shape[0]), pw), pw)["weights"]
        pred = (nodes[:, 0] + 1j * nodes[:, 1])[:, None]  # (K, 1) complex
        pm = P.predictive_moments(pred, w, noise_cov=np.diag([0.01, 0.02]))
        np.testing.assert_allclose(pm["mean_real"], [0.0, 0.0], atol=1e-10)
        np.testing.assert_allclose(pm["covariance_real"], np.diag([0.76, 0.77]), atol=1e-9)
        self.assertEqual(pm["mean_complex"].shape, (1,))

    def test_bimodal_branch_weights_external_labels(self):
        rule = P.gauss_legendre_box_nodes([[-2.5, 2.5], [0.0, 1.0], [0.0, 1.0]], [121, 5, 5])
        nodes, pw = rule["nodes"], rule["weights"]
        w0 = 0.2
        lik = np.log(np.exp(-0.5 * ((nodes[:, 0] - 1.0) / w0) ** 2)
                     + np.exp(-0.5 * ((nodes[:, 0] + 1.0) / w0) ** 2))
        norm = P.normalize_log_weights(P.posterior_log_weights(lik, pw), pw)
        labels = (nodes[:, 0] > 0).astype(int)
        bw = P.branch_weights(labels, norm["weights"], n_branches=2)["probabilities"]
        self.assertLess(abs(bw[0] - 0.5), 0.02)
        self.assertLess(abs(bw[1] - 0.5), 0.02)
        # unequal amplitudes -> weight ratio tracks the amplitudes
        lik2 = np.log(np.exp(-0.5 * ((nodes[:, 0] - 1.0) / w0) ** 2)
                      + 3.0 * np.exp(-0.5 * ((nodes[:, 0] + 1.0) / w0) ** 2))
        norm2 = P.normalize_log_weights(P.posterior_log_weights(lik2, pw), pw)
        bw2 = P.branch_weights(labels, norm2["weights"], n_branches=2)["probabilities"]
        # branch 1 is x>0, the amplitude-1 mode -> 1/4; branch 0 is x<0 -> 3/4
        self.assertLess(abs(bw2[1] - 0.25), 0.03)
        self.assertLess(abs(bw2[0] - 0.75), 0.03)

    def test_posterior_summary_keys(self):
        nodes, pw = self._grid(order=(41, 9, 9))
        lik = np.zeros(nodes.shape[0])
        pred = (nodes[:, 0] + 1j * nodes[:, 1])[:, None]
        labels = (nodes[:, 0] > 0).astype(int)
        summ = P.posterior_summary(nodes, lik, q=nodes[:, 0], pred=pred, branch_ids=labels)
        for key in ("weights", "log_evidence", "ess", "rejection", "theta", "q", "predictive", "branch"):
            self.assertIn(key, summ)
        self.assertEqual(summ["rejection"]["n_rejected"], 0)


class RiskBootstrapTest(unittest.TestCase):
    def test_population_risk_and_prevalence(self):
        v = np.array([0.01, 0.04, 0.09, 0.16])
        self.assertAlmostEqual(P.population_risk(v), 0.075, places=12)
        self.assertAlmostEqual(P.probability_sd_above(v, 0.25), 0.5, places=12)
        self.assertAlmostEqual(P.probability_sd_above(v, 0.1), 0.75, places=12)

    def test_paired_bootstrap_reduction_and_ci(self):
        rng = np.random.default_rng(7)
        base = rng.uniform(0.05, 0.2, size=400)
        new = 0.6 * base
        out = P.paired_bootstrap_risk(
            {"base": base, "targeted": new}, replicates=2000, seed=11, tau=0.3, baseline="base"
        )
        self.assertTrue(out["paired"])
        self.assertLess(abs(out["risk_reduction"]["targeted"] - 0.4), 0.05)
        lo, hi = out["risk_reduction_ci"]["targeted"]
        self.assertTrue(lo <= 0.4 <= hi)
        self.assertIn("prevalence", out)
        self.assertEqual(out["risk"]["targeted"], float(new.mean()))

    def test_paired_bootstrap_identical_methods_zero(self):
        v = np.linspace(0.02, 0.3, 200)
        out = P.paired_bootstrap_risk({"a": v, "b": v.copy()}, replicates=500, seed=1)
        self.assertAlmostEqual(out["risk_reduction"]["b"], 0.0, places=12)
        self.assertEqual(out["risk_reduction_ci"]["b"], (0.0, 0.0))


class ReferenceToleranceTest(unittest.TestCase):
    def test_loads_discovery_tolerances(self):
        tol = P.load_reference_tolerances()
        self.assertAlmostEqual(tol["reference_mean_tol"], 0.02, places=12)
        self.assertAlmostEqual(tol["reference_quantile_tol"], 0.03, places=12)
        self.assertAlmostEqual(tol["reference_branch_tol"], 0.02, places=12)
        self.assertAlmostEqual(tol["reference_variance_rel_tol"], 0.05, places=12)
        self.assertAlmostEqual(tol["reference_variance_abs_tol"], 1e-5, places=12)
        self.assertEqual(tol["quadrature_orders"], [33, 41])
        self.assertEqual(tol["bootstrap_replicates"], 5000)
        self.assertTrue(Path(tol["config_path"]).name == "discovery.json")

    def test_compare_identical_and_perturbed(self):
        ref = {"mean": 0.5, "variance": 0.04, "quantiles": [0.1, 0.5, 0.9], "branch_weights": [0.5, 0.5]}
        same = P.compare_posterior_to_reference(ref, dict(ref))
        self.assertTrue(same["passed"])

        bad = dict(ref)
        bad["mean"] = 0.55  # exceeds 0.02
        res = P.compare_posterior_to_reference(ref, bad)
        self.assertFalse(res["passed"])
        self.assertFalse(res["checks"]["mean"]["passed"])

    def test_variance_abs_tolerance_path(self):
        ref = {"variance": 1e-6}
        cand = {"variance": 2e-6}  # rel 1.0 > 0.05, abs 1e-6 <= 1e-5
        res = P.compare_posterior_to_reference(ref, cand)
        self.assertTrue(res["checks"]["variance"]["passed"])

    def test_convergence_metrics(self):
        good = P.check_convergence_metrics({"mesh_q": 0.01, "whitened_rms": 0.05, "whitened_max": 0.2, "linear_residual": 1e-11})
        self.assertTrue(good["passed"])
        bad = P.check_convergence_metrics({"linear_residual": 1e-3})
        self.assertFalse(bad["passed"])
        with self.assertRaises(KeyError):
            P.check_convergence_metrics({"unknown_metric": 0.0})


class SmcTest(unittest.TestCase):
    def test_analytic_gaussian(self):
        mu = np.array([0.3, -0.4, 0.5])
        s = 0.35

        def loglik(theta):
            return -0.5 * np.sum(((theta - mu) / s) ** 2, axis=1)

        res = P.adaptive_tempered_smc(loglik, [[-3, 3], [-3, 3], [-3, 3]], n_particles=4000, seed=1)
        self.assertTrue(res["converged"])
        self.assertAlmostEqual(res["lambda"], 1.0, places=9)
        mean = res["weights"] @ res["theta"]
        np.testing.assert_allclose(mean, mu, atol=0.04)
        var = res["weights"] @ (res["theta"] - mean) ** 2
        np.testing.assert_allclose(var, np.full(3, s**2), atol=0.03)
        self.assertTrue(np.isfinite(res["log_evidence"]))
        self.assertGreater(len(res["diagnostics"]), 0)

    def test_zero_information_single_step(self):
        const = 5.0

        def loglik(theta):
            return np.full(theta.shape[0], const)

        res = P.adaptive_tempered_smc(loglik, BOX, n_particles=2000, seed=2)
        self.assertEqual(res["step_count"], 1)
        self.assertAlmostEqual(res["log_evidence"], const, places=6)
        self.assertAlmostEqual(res["ess"], 2000.0, places=3)
        mean = res["weights"] @ res["theta"]
        # Monte-Carlo mean of the uniform prior: SE ~ sqrt(0.75/2000) ~ 0.019.
        np.testing.assert_allclose(mean, np.zeros(3), atol=0.08)
        self.assertAlmostEqual(res["diagnostics"][0]["pcn_acceptance"], 1.0, places=9)

    def test_diagnostics_complete(self):
        rng = np.random.default_rng(0)
        mu = rng.normal(size=3)

        def loglik(theta):
            return -0.5 * np.sum((theta - mu) ** 2, axis=1) / 0.3**2

        res = P.adaptive_tempered_smc(loglik, [[-2, 2]] * 3, n_particles=1500, seed=3, pcn_steps=3)
        keys = {
            "step", "lambda", "delta_lambda", "ess", "conditional_ess",
            "log_evidence_increment", "resampled", "ess_after_resample",
            "pcn_steps", "pcn_acceptance", "log_likelihood_evals",
        }
        for rec in res["diagnostics"]:
            self.assertTrue(keys.issubset(rec.keys()))
            self.assertLessEqual(rec["ess"], 1500 + 1e-6)
            self.assertGreaterEqual(rec["ess"], 0.0)
            self.assertTrue(0.0 <= rec["pcn_acceptance"] <= 1.0)
        lam = [rec["lambda"] for rec in res["diagnostics"]]
        self.assertTrue(all(b >= a for a, b in zip(lam, lam[1:])))

    def test_bimodal_1d_embedded_3d_branch_weights(self):
        w0 = 0.2

        def loglik(theta):
            x = theta[:, 0]
            return np.log(np.exp(-0.5 * ((x - 1.0) / w0) ** 2) + np.exp(-0.5 * ((x + 1.0) / w0) ** 2))

        res = P.adaptive_tempered_smc(
            loglik,
            [[-2.5, 2.5], [0.0, 1.0], [0.0, 1.0]],
            n_particles=4000,
            seed=4,
            readout=lambda theta: theta[:, 0],
            branch_fn=lambda theta: (theta[:, 0] > 0).astype(int),
        )
        probs = res["branch_weights"]["probabilities"]
        self.assertLess(abs(probs[0] - 0.5), 0.06)
        self.assertLess(abs(probs[1] - 0.5), 0.06)
        self.assertLess(abs(res["readout_statistics"]["mean"]), 0.1)

    def test_rejects_nonfinite_likelihood(self):
        def loglik(theta):
            out = -np.sum(theta**2, axis=1)
            out[theta[:, 0] < 0.0] = np.nan
            return out

        with self.assertRaises(ValueError):
            P.adaptive_tempered_smc(loglik, [[-2, 2]] * 3, n_particles=500, seed=5)

    def test_wrapper_matches_function(self):
        def loglik(theta):
            return -0.5 * np.sum(theta**2, axis=1) / 0.4**2

        smc = P.TemperedSMC(loglik, [[-3, 3]] * 3)
        res = smc.sample(seed=6, n_particles=1000)
        self.assertTrue(res["converged"])
        self.assertEqual(res["n_particles"], 1000)
        self.assertAlmostEqual(res["lambda"], 1.0, places=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
