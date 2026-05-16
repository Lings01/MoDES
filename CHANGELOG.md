# Changelog

All notable changes to MoDES will be documented in this file.

## v0.5.0-beta

### Added
- CLI: `modes run` and `modes validate-input` with JSON output
- Single-cell pseudobulk workflow (`from_pseudobulk` with sparse matrix support)
- Spatial region-pseudobulk (`from_spatial_pseudobulk`)
- MuData input (`from_mudata`) with pseudobulk aggregation
- cis-ATAC score gene-level aggregation (`conditional_mode="cis_score"`)
- `MoDESResult.save()` / `MoDESResult.load()` round-trip
- `run_by_context()` for multi-cell-type batch analysis
- Enhanced `filter()`: `states`, `genes`, `peaks`, `context`, `min_quality_score`
- Deterministic event_id via SHA1 hash
- `contrast=(ref,tgt)` parameter for explicit reference level
- `allow_poisson_fallback` / `allow_simplified_fallback` gating
- Enhanced artifact risk: `low_atac_depth`, `low_rna_depth`, `borderline_quality`
- Run params: version info, link counts in output
- Benchmark suites: simulated, semi-real spike-in, negative control, baseline comparison
- Benchmark `--quick` mode and standardized outputs
- PBMC 10k multiome real-data test
- Biological demo: chromatin priming detection
- CI: pytest (3.10/3.11), lint (ruff), build (wheel + twine)
- `docs/api_reference.md`, `docs/install_review.md`
- `docs/release_checklist.md`, `docs/v1.1_checklist.md`
- `pyproject.toml` with ruff config
- `requirements-dev.txt` with pytest, pytest-cov, ruff
- Output schema frozen (26 columns, user-friendly order)

### Changed
- Repository moved from alpha prototype to beta release candidate
- `confidence` → `state_confidence` (all outputs, API, tests)
- `artifact_like` removed as primary biological state; replaced by `artifact_risk` + `artifact_reason`
- Direction based on coefficient sign, not hardcoded FDR threshold
- Design matrix rank-deficient errors now include specific causes
- Link source normalized: `promoter`, `distal_250kb`

### Fixed
- `from_pseudobulk()` _group loss bug
- DataFrame boolean evaluation in `external_links`
- Silent empty events now raise `ValueError`
- Poisson fallback no longer mislabeled as NB in diagnostics
- HTML report escaping for user-controlled strings
- CPM normalization zero-library-size check
- `get_library_sizes()` zero-check

### Known Limitations
- RNA+ATAC only (no native protein layer)
- Binary condition only
- Native spatial graph planned (region-pseudobulk available)
- `state_confidence` is not a calibrated posterior
- Fixed-effect donor/batch (no random effects)

---

## v0.1.0-alpha

Initial release with core RNA+ATAC regulatory event state inference.
