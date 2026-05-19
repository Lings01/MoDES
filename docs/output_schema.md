# Output Schema

MoDES v2.0.0 produces a fixed-schema event table, a long-format modality evidence
table, conditional effects, model diagnostics, and run parameters.

## 1. event_table.tsv (Fixed Schema)

Main output: one row per candidate regulatory event.

| Field | Type | Description |
|---|---|---|
| `event_id` | str | Unique event identifier |
| `tf_name` | str | TF annotation (from motifs) |
| `peak_id` | str | Genomic interval (e.g. chr1:100-200) |
| `gene` | str | Target gene |
| `context` | str | Context label |
| `link_source` | str | Link provenance (promoter, distal_250kb, external) |
| `link_score` | float | Peak-gene link confidence (0-1) |
| `state_family` | str | Coarse state group (concordant, chromatin_primed, etc.) |
| `state` | str | Specific state (e.g. concordant_activation) |
| `state_assignment_score` | float | Evidence score × link_score (ranking only) |
| `raw_state_assignment_score` | float | Evidence score before link_score adjustment |
| `state_support_score` | float | Directional evidence score (NOT formal FDR) |
| `state_support_adjusted_score` | float | BH-style ranking adjustment |
| `supporting_modalities` | str | Modalities supporting this state (;) |
| `absent_modalities` | str | Modalities required to be absent (;) |
| `conflicting_modalities` | str | Modalities conflicting with this state (;) |
| `missing_modalities` | str | Required modalities not measured (;) |
| `atac_coef` | float | ATAC condition effect (logFC) |
| `atac_se` | float | ATAC standard error |
| `atac_pval` | float | ATAC p-value |
| `atac_fdr` | float | ATAC BH-corrected FDR |
| `atac_direction` | int | ATAC direction (+1/-1/0) |
| `rna_coef` | float | RNA condition effect (logFC) |
| `rna_se` | float | RNA standard error |
| `rna_pval` | float | RNA p-value |
| `rna_fdr` | float | RNA BH-corrected FDR |
| `rna_direction` | int | RNA direction (+1/-1/0) |
| `rna_after_atac_coef` | float | RNA effect after ATAC adjustment |
| `rna_after_atac_se` | float | Conditional SE |
| `rna_after_atac_pval` | float | Conditional p-value |
| `rna_after_atac_fdr` | float | Conditional FDR |
| `artifact_risk` | str | low / medium / high |
| `artifact_reason` | str | Semicolon-separated reasons |
| `quality_score` | float | Event quality (0-1) |

### Deprecated Fields (backward compatibility only)

| Field | Replacement |
|---|---|
| `state_confidence_deprecated` | Use `state_assignment_score` |
| `event_pval_deprecated` | Use `state_support_score` |
| `event_fdr_deprecated` | Use `state_support_adjusted_score` |

These are not calibrated probabilities. `state_confidence` is not a posterior.

## 2. event_modality_evidence.tsv (Long Format)

One row per event × modality combination. Core multi-modal evidence table.

| Field | Description |
|---|---|
| `event_id` | Links to event_table |
| `modality` | Modality name (rna, atac, h3k27ac_cuttag, protein, etc.) |
| `assay` | Assay type (RNA, ATAC, CUTTAG, PROTEIN, SPATIAL) |
| `target` | Molecular target (H3K27ac, H3K27me3, etc.) |
| `feature_id` | Feature measured (gene, peak, protein ID) |
| `role` | Regulatory role (transcript_output, chromatin_accessibility, activating_mark, protein_output, etc.) |
| `coef` | Condition effect (logFC) |
| `se` | Standard error |
| `pval` | Nominal p-value |
| `fdr` | BH-corrected FDR |
| `direction` | Effect direction (+1/-1/0) |
| `quality_score` | Per-modality quality |
| `model_used` | GLM model used |
| `converged` | Whether GLM converged |

## 3. conditional_effects.tsv (Diagnostic)

Multi-model conditional decomposition results. **Diagnostic only; not used for state assignment.**

| Field | Description |
|---|---|
| `event_id` | Links to event_table |
| `model_name` | e.g. rna_after_atac, rna_after_h3k27ac, protein_after_rna |
| `response_modality` | Response modality |
| `conditioning_modalities` | Conditioning modalities (;) |
| `condition_coef` | Condition coefficient after adjustment |
| `condition_se` | Standard error |
| `condition_pval` | P-value |
| `condition_fdr` | FDR |
| `attenuation` | Coef change relative to marginal |
| `converged` | GLM convergence |
| `model_used` | GLM family |

## 4. model_diagnostics.tsv

Per-feature GLM diagnostics for all modalities.

| Field | Description |
|---|---|
| `feature_id` | Feature identifier |
| `modality` | RNA, ATAC, CUTTAG, PROTEIN, etc. |
| `model_used` | nb_default_alpha, nb_fixed_alpha, poisson_fallback, etc. |
| `family` | negative_binomial, poisson |
| `converged` | Whether GLM converged |
| `dropped_covariates` | Whether covariates were dropped |
| `warning` | Model warning message |

## 5. run_params.tsv

Key-value run parameters (condition, FDR threshold, versions, sample counts).

## Optional Outputs

- `event_state_confidence.tsv` — State classifier output (backward compat)
- `event_layer_effects.tsv` — Per-layer ATAC/RNA effects
- `event_evidence_vectors.tsv` — D_e evidence vectors
- `event_network.graphml` — TF-peak-gene network
- `report.html` — HTML summary report
