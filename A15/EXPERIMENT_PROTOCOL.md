# A15 frozen confirmatory protocol
SHA256: eddc3a97a3c39f0bdead44956a3bffadaaa103d2a8bb2568e0d1b9b9658cb38b

Source of requirements: A15_prompt.md; all thresholds inherited from the approved plan and configs/discovery.json.

Selected two-circle family: radius=0.1, epsilon_base=12.0, target center=(0.18,0), bright center=(target_x-2r-s,0). Independent uniform bounds [[0.02, 0.16], [-3.0, 3.0], [0.001, 0.08]]. Full width retained: mass=1 within selected family; 1/6 under an explicitly uniform six-template discovery mixture. No universal Maxwell prior probability is asserted.

Time exp(-i omega t), epsilon=epsilon_real+i gamma k_ref/k; unit E_z plane wave propagating +x. ROI is fixed target disk. Q=natural log of mean ROI |E|^2 relative to incident.

k0=8.625; frequencies=[8.19375, 8.625, 9.05625]; RMS per REAL component=0.2025837700889245; primary sigma=0.00405167540177849. Shared absolute noise for all channels, 1/2/5% boundaries.

Base channels 0/1; same-frequency fixed 2/3; optimized [4, 8]; random [9, 8]; neighbor frequency index 2; near 16/17. Every addition has exactly 2 complex outputs. Receiver positions and design-proxy predictions in JSON. Local full-wave proxy prediction is not a successful two-mode explanation.

seeds 11/29/47; validation 128 per seed; locked test 512 per seed. Truth is continuous prior, not quadrature nodes. 33^3/41^3 positive-weight quadrature, independent SMC diagnostic. No table/node failure deletion.

Gate A: SD(Q|y)>ln1.25 AND 90% energy interval ratio>=1.5. Each seed prevalence >=.1 and pooled Wilson lower95>=.1, conditional on valid physical/reference accuracy. Gate B also requires two-mode errors below noise/task margin, predicted targeted advantage >=10 percentage points vs fixed AND optimized same-frequency baselines, paired query bootstrap with Holm adjustment and consistent seed effects, plus coupling attribution.

Frozen numerical backend adjustment: VIE retained and audited; exact-circle expansion order20 with independent order24 checks selected because VIE n<=128 did not meet the noise error budget. This does not alter physical family or gates.

Flow blocked until Gate B PASS and cheap simple posterior models fail. All readouts use same physical theta. Failure/incomplete results must remain visible. Resume always rereads full A15_prompt.md before this protocol and RESUME_STATE.md.
