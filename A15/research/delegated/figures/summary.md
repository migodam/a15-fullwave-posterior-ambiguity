Done. Script runs clean under the required interpreter, writes only to the allowed paths, and touches nothing in `src/`, `configs/`, `EXPERIMENT_PROTOCOL.md`, or `results/`.

## Artifact list

Script: [plot_results.py](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/experiments/plot_results.py>)

Figures ([figures/](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/figures>), duplicated under `research/delegated/figures/`):

| id | file | status |
| --- | --- | --- |
| Fig1 | fig01_geometry_roi_acquisition.png | GENERATED (protocol.json: geometry/ROI/receivers/roles) |
| Fig2 | fig02_frequency_spectrum.png | PROVISIONAL (discovery sweep; bright/dark labels withheld pending modal audit) |
| Fig3 | fig03_mode_diagnostic.png | NOT RUN |
| Fig4 | fig04_parameter_map.png | PROVISIONAL (discovery 9³ screening; not prevalence) |
| Fig5 | fig05_parameter_posterior.png | NOT RUN |
| Fig6 | fig06_Q_posterior.png | NOT RUN |
| Fig7 | fig07_material_field_pair.png | NOT RUN |
| Fig8 | fig08_query_sd_distribution.png | NOT RUN |
| Fig9 | fig09_risk_comparison.png | NOT RUN |
| Fig10 | fig10_coupling_ablation.png | NOT RUN |
| Fig11 | fig11_modal_decomposition.png | NOT RUN |

Manifest: [FIGURES_MANIFEST.md](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Flow Matching/A15/research/delegated/figures/FIGURES_MANIFEST.md>) and `FIGURES_MANIFEST.json`.

## Missing data (why the rest is NOT RUN)

- Fig3, Fig11: `results/modal_audit.npz` (+ `modal_audit.json`)
- Fig5, Fig6, Fig7: `results/table_41_coupled/{prior.npz,Q.npy,y.npy}`, `results/queries_validation_11_coupled.npz`, `results/validation_41_coupled/summary.json`
- Fig8, Fig9: `results/validation_41_coupled/summary.json`
- Fig10: `results/validation_41_coupled/summary.json`, `results/validation_41_uncoupled/summary.json`

Each missing case is written as a clearly labelled `NOT RUN` placeholder PNG (no invented data). When the confirmatory files appear, rerunning the script auto-fills them and marks them `PROVISIONAL` until `reference_validation_gate.json` passes (currently absent → gate reads false).

## Errors

None. No runtime errors, no missing-glyph/font warnings, and every one of the eleven figure code paths was exercised: Figs 5/6/7/8/9/10 and 3/11 were dry-run against synthetic fixtures in a temp dir (Fig7 including real `CylinderModel.solve(fields=True)` / `.field` solves) and returned clean, so the parent rerun should not hit an unfilled path.

Behavioural notes preserved from the task: representative posterior uses the maximum-variance validation query (never test), weights are `log prior − 0.5·‖y_obs−y_pred‖²/σ²` with the protocol σ, Fig7 selects the 10%/90%-Q nodes by highest posterior weight and is labelled "low-Q node / high-Q node" (explicitly not asserted as branches).