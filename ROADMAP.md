# Roadmap

## Current: v2.0.0

Completed:
- RNA+ATAC regulatory event-state inference (validated core)
- Bulk RNA+ATAC support
- Single-cell pseudobulk support (experimental)
- Spatial region-pseudobulk (experimental)
- Multi-modal data structures (CUT&Tag, Protein, Spatial)
- Grammar-driven state classification for up to 17 states
- CUT&Tag pipeline (effect estimation + state classification)
- Protein pipeline (effect estimation + state classification)
- MuData (.h5mu) input (experimental)
- Interval-index event construction (O(G log P))
- CLI (`modes run`, `modes validate-input`)
- Benchmark suites (6 types)
- CI (pytest 3.10/3.11 + lint + build)
- Output schema frozen (26 columns) + long-format modality evidence

## v2.0 Hardening

- [ ] Multi-condition pipeline validation
- [ ] Conditional decomposition for CUT&Tag/Protein (RNA ~ Condition + ExtraModality)
- [ ] Pseudotime/continuous-dose support
- [ ] Calibrated posterior probabilities (state_confidence)
- [ ] Real multi-modal validation datasets (multi-omics PBMC, spatial + cut&tag)

## v2.1+

- Multi-condition contrasts (>2 groups) with contrast matrix
- Native spatial graph states validated on real spatial transcriptomics
- Full multi-layer event model (RNA + ATAC + Protein + Epigenomic + Spatial)
- Calibrated Bayesian posteriors with external benchmarks
