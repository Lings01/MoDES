# MoDES: Multi-Omics Discordance/Event State inference

Multi-omics discordance-guided decomposition of regulatory event states.

## Overview

MoDES is a statistical framework that treats regulatory events — not genes, peaks, or
clusters — as the primary unit of analysis. Instead of compressing multi-omics layers
into a single embedding, MoDES uses cross-modality concordance and discordance to
classify each event into an interpretable regulatory state.

MoDES reports two layers of information:

1. **biological state**: `concordant`, `chromatin_primed`, `rna_only`,
   `discordant_opposite`, `null`
2. **artifact risk**: `low`, `medium`, `high`

The current release is **MoDES-RA v0.1.0-alpha** (RNA + ATAC).

## Current status

MoDES-RA is an alpha-stage RNA+ATAC bulk/pseudobulk prototype.

Current scope:
- RNA + ATAC only
- binary condition only
- bulk or externally generated pseudobulk recommended
- event candidates require `external_links`, `genome_annotation`, or `tss_map`
- native protein model is planned
- native spatial graph model is planned

## Installation

```bash
pip install -e .
pip install -r requirements-dev.txt  # for running tests
```

## Quick start

### Option A: external peak-gene links (recommended)

```python
import pandas as pd
from modes import MoDES, MoDEData

data = MoDEData.from_matrices(
    rna_counts="rna_counts.tsv",
    atac_counts="atac_peaks.tsv",
    metadata="sample_metadata.tsv",
    condition_col="condition",
    index_col=0,
)

links = pd.read_csv("peak_gene_links.tsv", sep="\t")
modes = MoDES(
    data=data,
    condition_col="condition",
    external_links=links,
)
result = modes.run()
```

### Option B: GTF genome annotation

```python
modes = MoDES(
    data=data,
    condition_col="condition",
    genome_annotation="genes.gtf",
)
result = modes.run()
```

### View and export results

```python
print(result.summary())
result.to_tsv("output/")
result.to_graphml("output/event_network.graphml")
result.to_report("output/report.html")
```

### Minimal example

```bash
python examples/minimal_bulk/run_minimal.py
```

## Algorithm

1. **Event Candidate Construction** — link peaks to target genes via proximity,
   external links, or motif annotation
2. **Effect Size Estimation** — NB GLM with empirical Bayes variance moderation
   for ATAC and RNA condition effects
3. **Conditional Decomposition** — RNA condition effect after adjustment for the
   linked ATAC peak
4. **Evidence Vector Construction** — `D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality]`
5. **State Classification** — rule-based initialization followed by empirical
   Bayes confidence refinement

## Input formats

- **Bulk**: TSV/CSV count matrices + metadata (recommended)
- **Single-cell**: experimental in v0.1.0. Aggregate externally to pseudobulk
  first, then load with `MoDEData.from_matrices()`.
  `MoDEData.from_pseudobulk()` supports simple AnnData objects.
- **Spatial**: planned. v0.1.0 does not model spatial coordinates or
  neighborhood graphs.

## Output files

| File | Description |
|---|---|
| `event_table.tsv` | Main output: one row per peak-gene regulatory event |
| `event_state_confidence.tsv` | State, state_confidence, artifact_risk per event |
| `event_layer_effects.tsv` | ATAC/RNA layer effect estimates |
| `event_evidence_vectors.tsv` | Evidence vectors used for state classification |
| `model_diagnostics.tsv` | Marginal ATAC/RNA GLM diagnostics (conditional model diagnostics planned) |
| `run_params.tsv` | Run parameters and data summary |
| `event_network.graphml` | Optional: TF-peak-gene network |
| `report.html` | Optional: HTML report |

### `event_table.tsv` fields

| Field | Description |
|---|---|
| `event_id` | Unique event identifier |
| `tf_name` | Optional TF annotation |
| `peak_id` | ATAC peak |
| `gene` | Target gene |
| `state` | Biological event state |
| `state_confidence` | State confidence, range 0 to 1 |
| `artifact_risk` | Technical risk flag: `low` / `medium` / `high` |
| `artifact_reason` | Semicolon-separated risk reasons |
| `event_pval` | Event-level p-value |
| `event_fdr` | Event-level BH-adjusted FDR |
| `atac_coef` | Condition effect on ATAC peak |
| `atac_pval` / `atac_fdr` | ATAC layer significance |
| `rna_coef` | Condition effect on RNA gene |
| `rna_pval` / `rna_fdr` | RNA layer significance |
| `rna_after_atac_coef` | RNA condition effect after linked-peak adjustment |
| `rna_after_atac_pval` / `rna_after_atac_fdr` | Conditional decomposition significance |
| `quality_score` | Event quality score |

> `model_diagnostics.tsv` currently reports marginal ATAC and RNA GLM diagnostics.
> Conditional RNA-after-ATAC diagnostics are summarized in `event_table.tsv`;
> full conditional diagnostics will be added in a future version.

## Citation

MoDES: Multi-Omics Discordance-guided decomposition of regulatory event states.

## License

MIT
