# Changelog

All notable changes to MoDES will be documented in this file.

## v0.1.0-alpha

### Added
- RNA+ATAC regulatory event state inference framework (MoDES-RA)
- Five biological states: `concordant`, `chromatin_primed`, `rna_only`,
  `discordant_opposite`, `null`
- Three-tier artifact risk assessment: `low`, `medium`, `high`
- `artifact_reason` field with human-readable risk explanations
- Event-level p-value and BH-corrected FDR (`event_pval`, `event_fdr`)
- NB GLM with empirical Bayes variance moderation for effect estimation
- Conditional decomposition (RNA after linked-peak adjustment)
- Four-tier GLM fallback strategy with full diagnostic tagging
- `model_diagnostics.tsv` output (marginal ATAC/RNA diagnostics)
- Design matrix rank deficiency detection
- Binary condition enforcement with clear error messages
- Zero library size detection
- Gene/peak coordinate coverage warnings
- `MoDEData` unified data container (TSV, DataFrame, AnnData)
- `MoDEData.from_matrices()` — bulk TSV/CSV/DataFrame input
- `MoDEData.from_anndata()` — single AnnData input
- `MoDEData.from_pseudobulk()` — pseudobulk aggregation from single-cell
- `MoDESResult.to_tsv()` — multi-file TSV export
- `MoDESResult.to_graphml()` — Cytoscape/Gephi-compatible network export
- `MoDESResult.to_report()` — self-contained HTML report
- `MoDESResult.filter()` — result filtering by state, confidence, FDR, artifact risk
- `MoDESResult.summary()` — text summary
- State confidence (`state_confidence`) replacing uncalibrated "posterior" naming
- Strict integration test with deterministic ground truth recovery
- Direct StateClassifier unit test (GLM-independent)
- CI import smoke test
- GitHub Actions CI workflow (Python 3.10, 3.11)
- `examples/minimal_bulk/` with real gene symbols and external links
- Split `requirements.txt` (runtime) / `requirements-dev.txt` (testing)

### Limitations
- Binary condition only (two-group comparison)
- RNA + ATAC only (no native protein or methylation layer)
- Bulk or externally-generated pseudobulk recommended
- Event candidates require `external_links`, `genome_annotation`, or `tss_map`
- No native spatial graph model (spatial data usable only via pseudobulk)
- No native cell-level mixed model
- No CLI interface
- Fixed-effect donor/batch (no random effects)
- Single linked-peak conditional model (no cis-ATAC score aggregation)
- `state_confidence` is not a calibrated posterior probability
