# Changelog

All notable changes to MoDES will be documented in this file.

## v2.0.0

### Added
- Multi-modal data structures: `MoDEData.modalities` dict + `modality_specs`
- Grammar-driven `StateClassifier` with 17 possible states across 5 grammars
- RNA+ATAC core (5 states): concordant, chromatin_primed, rna_only, discordant_opposite, null
- Epigenomic extension (4 states): epigenomic_concordant, active_enhancer_primed, repressive_concordant, repressive_primed, mark_only
- Protein extension (4 states): full_activation, protein_buffered, protein_memory, protein_opposite
- Spatial extension (4 states): spatial_region_specific, spatial_niche_driven, cell_intrinsic, spatial_edge_artifact
- CUT&Tag target registry: H3K27ac, H3K4me1, H3K4me3, H3K27me3, H3K9me3, H3K36me3, CTCF, RAD21, TF
- Generic modality effect estimation (`estimate_modality_effects`)
- `EventCandidateBuilder`: O(G log P) interval index for event construction
- Long-format `event_modality_evidence.tsv` output
- CUT&Tag benchmark: RNA + H3K27ac CUT&Tag state recovery
- Protein benchmark: RNA + ATAC + Protein state recovery
- Tutorial notebook (Jupyter) with executed cells
- Real-data tests: 10x PBMC, GSE166188, Visium spatial
- `ModalitySpec` base class for abstract modal specification

### Changed
- Core pipeline: `self.effects` dict replaces hardcoded atac_effects/rna_effects
- `EvidenceBuilder`: accepts extra_modality_effects for CUT&Tag/Protein
- `StateClassifier`: grammar-driven, dynamic state selection based on available modalities
- Dedup key now includes tf_name (peak+gene+tf) to preserve TF-specific events
- Quality score: depth_score uses percentile-based reference (fixed saturation bug)
- Quality combination: proper weighted average (0.6*atac + 0.4*rna, fixed 0.76 max bug)

### Fixed
- TF dedup: `peak_id|gene` → `peak_id|gene|tf_name` preserves TF-specific events
- depth_score saturation: `log(mean+1)/log(ref+1)` with 95th percentile reference
- Quality combination: `0.6*0.6*atac + 0.4*rna` → `0.6*atac + 0.4*rna`

### Known Limitations
- RNA+ATAC core is the only validated layer; CUT&Tag/Protein states are experimental
- Binary condition only (two-group comparison)
- `state_confidence` is not a calibrated posterior
- Fixed-effect donor/batch (no random effects)
- Conditional decomposition (Step 3) only decomposes RNA after ATAC

---

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

## v1.0.0-rc.1

### Added
- Release candidate for stable RNA+ATAC MoDES-RA
- API freeze candidate (semantic versioning from v1.0.0)
- Output schema freeze (26 columns, user-friendly order)
- Fresh clone / install review evidence
- Benchmark and biological demo documentation
- `pyproject.toml` as primary build entry

### Known limitations
- RNA+ATAC only (no native protein layer)
- Binary condition only (two-group comparison)
- Native spatial graph planned (region-pseudobulk available)
- `state_confidence` is not a calibrated posterior
- Fixed-effect donor/batch (no random effects)
