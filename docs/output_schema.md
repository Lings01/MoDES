# Output Schema

MoDES v2.0.0 produces TSV files, an optional GraphML network, and an optional HTML report.

The `event_table.tsv` schema is **frozen** for the v2.0 release.

## File Inventory

| File | Format | Rows | Description |
|---|---|---|---|
| `event_table.tsv` | TSV | n_events | Main output table |
| `event_state_confidence.tsv` | TSV | n_events | State confidence and artifact risk |
| `event_layer_effects.tsv` | TSV | n_events | Per-layer effect estimates |
| `event_evidence_vectors.tsv` | TSV | n_events | Evidence vectors (z-scores) |
| `model_diagnostics.tsv` | TSV | n_features | GLM diagnostics per feature |
| `run_params.tsv` | TSV | n_params | Run parameters |
| `event_network.graphml` | GraphML | n_edges | TF-peak-gene network (optional) |
| `report.html` | HTML | — | Visual summary (optional) |

## event_table.tsv

One row per candidate regulatory event (peak-gene pair in a given context).

| Field | Type | Range | Nullable | Description |
|---|---|---|---|---|
| `event_id` | str | — | No | Unique identifier |
| `tf_name` | str | — | Yes | TF motif annotation |
| `peak_id` | str | — | No | ATAC peak identifier |
| `gene` | str | — | No | Target gene name |
| `context` | str | — | Yes | Cell type / condition label |
| `state` | str | {concordant, chromatin_primed, rna_only, discordant_opposite, null} | No | Biological state |
| `state_confidence` | float | [0, 1] | No | Rule-refined confidence score |
| `artifact_risk` | str | {low, medium, high} | No | Technical risk flag |
| `artifact_reason` | str | — | Yes | Semicolon-separated reasons |
| `event_pval` | float | [0, 1] | No | Event-level p-value |
| `event_fdr` | float | [0, 1] | No | Event-level BH-corrected FDR |
| `quality_score` | float | [0, 1] | No | Event quality score |
| `atac_coef` | float | ℝ | No | Condition effect on ATAC (logFC scale) |
| `atac_se` | float | (0, ∞) | No | Standard error |
| `atac_pval` | float | [0, 1] | No | ATAC p-value |
| `atac_fdr` | float | [0, 1] | No | ATAC BH-corrected FDR |
| `atac_direction` | int | {-1, 0, 1} | No | ATAC effect direction |
| `rna_coef` | float | ℝ | No | Condition effect on RNA (logFC scale) |
| `rna_se` | float | (0, ∞) | No | Standard error |
| `rna_pval` | float | [0, 1] | No | RNA p-value |
| `rna_fdr` | float | [0, 1] | No | RNA BH-corrected FDR |
| `rna_direction` | int | {-1, 0, 1} | No | RNA effect direction |
| `rna_after_atac_coef` | float | ℝ | No | RNA effect after linked-peak adjustment |
| `rna_after_atac_se` | float | (0, ∞) | No | Standard error |
| `rna_after_atac_pval` | float | [0, 1] | No | Conditional p-value |
| `rna_after_atac_fdr` | float | [0, 1] | No | Conditional BH-corrected FDR |

### State Definitions

| State | ATAC | RNA | Direction | Meaning |
|---|---|---|---|---|
| `concordant` | sig | sig | same | Intact cis-regulatory activation |
| `chromatin_primed` | sig | not sig | — | Chromatin open, transcription not started |
| `rna_only` | not sig | sig | — | RNA change not from local chromatin |
| `discordant_opposite` | sig | sig | opposite | Opposite direction across layers |
| `null` | not sig | not sig | — | No significant change |

### Event-Level p-value

Computed from state-specific combination of ATAC/RNA p-values:

- `concordant` / `discordant_opposite`: `max(atac_pval, rna_pval)`
- `chromatin_primed`: `atac_pval`
- `rna_only`: `rna_pval`
- `null`: `1.0`

BH correction is applied across all events to produce `event_fdr`.

## model_diagnostics.tsv

One row per fitted feature (peak or gene).

| Field | Type | Description |
|---|---|---|
| `feature_id` | str | Peak or gene identifier |
| `modality` | str | ATAC or RNA |
| `model_used` | str | GLM variant used (nb_default_alpha, nb_fixed_alpha, poisson_fallback, nb_simple_fallback) |
| `family` | str | Distribution family (negative_binomial, poisson) |
| `alpha` | float | NB dispersion parameter (None for Poisson) |
| `alpha_estimated` | bool | Whether alpha was estimated per-feature |
| `converged` | bool | Whether the GLM converged |
| `dropped_covariates` | bool | Whether covariates were dropped (simplified fallback) |
| `warning` | str | Diagnostic warning message |

> `model_diagnostics.tsv` currently covers marginal ATAC and RNA GLMs only.
> Conditional RNA-after-ATAC diagnostics are summarized in `event_table.tsv`.

## event_network.graphml

A Cytoscape/Gephi-compatible network with nodes (gene, peak, tf) and edges
(peak→gene event, tf→peak motif). Edge attributes include `state`,
`state_confidence`, `artifact_risk`, `event_fdr`, `atac_coef`, and `rna_coef`.
