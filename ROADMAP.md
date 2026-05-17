# Roadmap

## Current: v1.0.0-rc.1

Completed:
- RNA+ATAC regulatory event-state inference
- Bulk RNA+ATAC support
- Single-cell pseudobulk support
- Spatial region-pseudobulk (experimental)
- MuData input (experimental)
- CLI (`modes run`, `modes validate-input`)
- Benchmark suites (4 types: synthetic, semi-real, negative control, baseline)
- CI (pytest 3.10/3.11 + lint + build)
- Output schema frozen (26 columns)
- Public API reference documented
- Fresh install review completed

## v1.0.0 Release Candidate

- [ ] Version/docs consistency audit
- [ ] Full benchmark summary with metrics
- [ ] Real-data demo documentation
- [ ] External install review sign-off
- [ ] Release checklist all-green

## v1.0.0

Stable RNA+ATAC release. Requirements:
- All v1.0.0-rc.1 gates passed
- CI green on Python 3.10/3.11
- All benchmarks produce standardized outputs
- External review confirmed

## Post-1.0

### v1.1
- Protein layer (MoDES-RAP): `full_activation`, `protein_buffered`, `protein_memory`
- `MoDEData` protein support

### v1.2
- Native spatial graph model (`SpatialMoDEData`)
- Spatial states: `spatial_region_specific`, `spatial_niche_driven`, `spatial_artifact_edge`

### v1.3
- Multi-condition contrasts (beyond binary)
- Continuous pseudotime / time-course lag inference

### v2.0
- Full multi-layer event model (RNA + ATAC + Protein + Methylation + Spatial)
- Calibrated Bayesian posteriors
