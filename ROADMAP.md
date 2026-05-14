# Roadmap

## v0.1.0-alpha (current)
RNA+ATAC bulk / pseudobulk MVP with 76 passing tests, state classification,
artifact risk, and full output schema.

## v0.2.0-alpha
- Real-data smoke tests (PBMC multiome)
- Synthetic benchmark with confusion matrix and per-state metrics
- Semi-real spike-in benchmark

## v0.3.0-alpha
- Single-cell pseudobulk workflow hardening
- `from_pseudobulk()` robustification (edge cases, MuData support)
- `run_by_context()` for multi-cell-type batch analysis
- `examples/singlecell_pseudobulk/`

## v0.4.0-alpha
- Spatial region-pseudobulk support (`from_spatial_pseudobulk()`)
- `context = spatial_region` in output
- Spatial documentation and example

## v0.5.0-beta
- CLI interface (`modes run`, `modes validate-input`)
- `pyproject.toml` migration
- Documentation freeze (`docs/` directory)
- API stability policy published
- Coverage ≥ 80%

## v0.8.0-beta
- Full benchmark suite (synthetic, semi-real, negative control, baseline comparison)
- Lint / type-check / format CI gates
- Package build verification
- Real biological demo (time-course, differentiation, or treatment dataset)

## v1.0.0
Stable RNA+ATAC release. Requirements:
- [ ] Bulk RNA+ATAC stable
- [ ] Single-cell pseudobulk stable
- [ ] Spatial region-pseudobulk stable (or explicitly out of scope)
- [ ] CLI operational
- [ ] Benchmark suite complete
- [ ] Real data demo done
- [ ] Output schema frozen
- [ ] API frozen
- [ ] Coverage ≥ 85%
- [ ] CI green on Python 3.10/3.11/3.12
- [ ] Full documentation (`docs/`)

## v1.1
- Protein layer (MoDES-RAP): `full_activation`, `protein_buffered`, `protein_memory`
- `MoDEData` protein support

## v1.2
- Native spatial graph model
- `SpatialMoDEData` with coordinates and neighborhood graph
- Spatial states: `spatial_region_specific`, `spatial_niche_driven`, `spatial_artifact_edge`

## v1.3
- Multi-condition contrasts (beyond binary)
- Continuous pseudotime / time-course lag inference
- ATAC→RNA delay estimation

## v2.0
- Full multi-layer event model (RNA + ATAC + Protein + Methylation + Spatial)
- Calibrated Bayesian posterior probabilities
- Best-in-class peak-gene linking integration
- Cell-level mixed models
