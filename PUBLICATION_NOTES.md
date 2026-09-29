# Publication notes

This curated copy preserves numerical results and scientific gates. Research originals remain the local source of truth. No new experiment was run for publication.

Included: theory, frozen protocol, source, tests, reports, eleven figures, per-query metrics, SMC/adaptive diagnostics, and all three raw full-wave tables.

Excluded: bytecode, redundant source ZIPs (which embed local transport configuration), full worker transcripts, duplicate worker PNGs and build helper. No credentials are published. Worker review summaries and literature retrieval errors remain.

Local LAN login/address and SSH identity filename are replaced with placeholders in the following files; scientific data/configuration unchanged:
- `reports/PLAN_AND_AUDIT.md`
- `src/remote.py`

Original RUN_MANIFEST and table metadata describe the original execution, including original source hashes; PUBLIC_MANIFEST describes this curated snapshot. Do not mistake publication substitutions for original run provenance. Original source archives are excluded; local transport must be configured before remote use. Local numerical runs do not require remote.py.

All original result tables retain their original file SHA256. Tests in the copy are rerun for packaging integrity, not as new scientific acceptance. The archive is about 200 MB; binary files cannot necessarily be interpreted by a web-only reader. Reports and JSON expose the necessary evidence without binary execution.
