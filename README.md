# MoDES: Multi-Omics Discordance/Event State inference

![tests](https://github.com/Lings01/MoDES/actions/workflows/tests.yml/badge.svg)

A statistical framework that treats regulatory events — not genes, peaks, or clusters —
as the primary unit of analysis. Instead of compressing multi-omics layers into a single
embedding, MoDES uses cross-modality concordance and discordance to classify each event
into an interpretable regulatory state with artifact risk assessment.

---

## 1. Principles

### 1.1 Core Idea

Traditional multi-omics analysis typically performs differential analysis on each modality
separately and then cross-validates — "RNA found DEGs, let's check if nearby peaks are
open." This implicitly treats one modality as primary and others as validation layers.

MoDES takes a different approach:

> **RNA, ATAC, and protein are not primary vs. secondary layers; they are different
> observation layers of the same regulatory event.**

- **Concordance** — complete regulatory chain activation: chromatin → transcription → protein
- **ATAC-only** (chromatin priming) — chromatin is open but transcription has not started
- **RNA-only** (trans-driven) — RNA change is not explained by local chromatin, possibly
  from trans regulation, RNA stability changes, or unmeasured regulatory layers
- **RNA-protein discordance** — post-transcriptional buffering, protein memory, or
  technical artifacts
- **Spatial local concordance** — niche-driven regulation

The fundamental improvement is in **data utilization**: rather than compressing multi-omics
into an embedding, MoDES uses discordance patterns between layers to infer *at which layer*
a regulatory event occurs.

### 1.2 Unit of Analysis: the Regulatory Event, not the Gene

A regulatory event is defined as:

```
e = (TF, enhancer/peak, target_gene, context)
```

For example:

> **STAT1 motif / enhancer chr1:100-200 accessibility ↑**
> → **IFIT3 RNA ↑**
> → occurring in **disease monocytes**

This is closer to real regulatory biology than gene-level differential expression or
peak-level differential accessibility alone.

### 1.3 Biological State Classification

MoDES-RA currently classifies RNA+ATAC regulatory events into five biological states:

| State | Pattern | Biological Interpretation |
|---|---|---|
| `concordant` | ATAC↑ RNA↑ | Local chromatin opening drives transcriptional activation (intact cis-regulatory chain) |
| `chromatin_primed` | ATAC↑ RNA→ | Chromatin is open but transcription has not started (epigenetic priming) |
| `rna_only` | ATAC→ RNA↑ | RNA change not explained by local chromatin (trans regulation, RNA stability, unmeasured layer) |
| `discordant_opposite` | ATAC↑ RNA↓ or ATAC↓ RNA↑ | Opposite directions across layers; may indicate complex regulation or technical issues |
| `null` | ATAC→ RNA→ | No significant change under the current condition |

Each event also carries an **artifact_risk** flag:

| Risk Level | Meaning |
|---|---|
| `low` | Good data quality; results are trustworthy |
| `medium` | Some quality concerns present |
| `high` | Single-modality signal with low quality score; interpret with caution |

This two-layer design (biological state + artifact risk) is more informative than labeling
events as "artifact": an event can be both `chromatin_primed` *and* `artifact_risk = high`,
preserving the biological signal while flagging quality concerns.

---

## 2. Current Status

**MoDES-RA v0.5.0-beta** — RNA+ATAC regulatory event-state inference.

> MoDES v0.5.0-beta currently implements the RNA+ATAC layer only.
> Protein-related and native spatial-graph states are conceptual extensions
> planned for post-1.0 releases.

| Capability | Status |
|---|---|
| Bulk RNA+ATAC | ✅ Supported |
| Single-cell pseudobulk | ⚠️ Experimental |
| Spatial region-pseudobulk | ⚠️ Experimental |
| MuData (.h5mu) input | ⚠️ Experimental |
| Sparse matrix aggregation | ⚠️ Experimental |
| Native spatial graph | 🔮 Planned (v1.2) |
| Protein layer | 🔮 Planned (v1.1) |
| Multi-class condition | 🔮 Planned |
| Time / pseudotime delay | 🔮 Planned |

---

## 3. Installation

```bash
# Base installation
pip install -e .

# Run tests (optional)
pip install -r requirements-dev.txt
python -m pytest -q
```

Dependencies: `numpy`, `scipy`, `pandas`, `statsmodels`, `anndata`, `matplotlib`, `seaborn`,
`networkx`.

---

## 4. Usage

### 4.1 Quick Start

```python
import pandas as pd
from modes import MoDES, MoDEData

# 1. Load data
data = MoDEData.from_matrices(
    rna_counts="rna_counts.tsv",
    atac_counts="atac_peaks.tsv",
    metadata="sample_metadata.tsv",
    condition_col="condition",
    index_col=0,
)

# 2. Load peak-gene links (recommended approach)
links = pd.read_csv("peak_gene_links.tsv", sep="\t")

# 3. Run MoDES
modes = MoDES(
    data=data,
    condition_col="condition",
    external_links=links,
)
result = modes.run()

# 4. View results
print(result.summary())

# 5. Export
result.to_tsv("output/")
result.to_graphml("output/network.graphml")
result.to_report("output/report.html")
```

### 4.2 Step-by-Step Execution

Each step can be run independently for debugging and interactive analysis:

```python
modes = MoDES(data=data, condition_col="condition", external_links=links)

# Step 1: Build candidate events
events = modes.build_events()

# Step 2: Estimate ATAC and RNA effects
atac_effects, rna_effects = modes.estimate_effects()

# Step 3: Conditional decomposition (RNA after ATAC)
conditional = modes.decompose()

# Step 4: Build evidence vectors
evidence = modes.build_evidence()

# Step 5: Classify states
states = modes.classify_states()

# Assemble final results
result = modes._assemble_results()
```

### 4.3 Filtering Results

```python
# View only concordant events
conc = result.filter(state="concordant")

# Exclude high artifact-risk events
clean = result.filter(exclude_high_artifact=True)

# Filter by event-level FDR
sig = result.filter(max_event_fdr=0.1)

# Combine filters
trusted = result.filter(
    state="concordant",
    min_confidence=0.8,
    max_event_fdr=0.05,
    exclude_high_artifact=True,
)
```

### 4.4 Run the Example

```bash
python examples/minimal_bulk/run_minimal.py
```

---

## 5. Parameter Reference

### 5.1 MoDES Constructor

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | `MoDEData` | required | Input data container |
| `condition_col` | `str` | required | Column in `data.obs` specifying the binary condition |
| `covariate_cols` | `list[str]` | `[]` | Additional covariate column names |
| `donor_col` | `str` | `None` | Donor/replicate identifier (treated as fixed effect) |
| `batch_col` | `str` | `None` | Batch identifier (treated as fixed effect) |
| `fdr_threshold` | `float` | `0.1` | FDR threshold for significance calls in state classification |
| `genome_annotation` | `str` | `None` | Path to GTF/GFF file for gene TSS coordinates |
| `external_links` | `DataFrame` | `None` | Pre-computed peak-to-gene links (recommended) |
| `motif_annotation` | `DataFrame` | `None` | Peak-to-TF motif mapping |
| `tss_map` | `dict` | `None` | Manual `gene → (name, chr, tss_pos)` mapping |

### 5.2 `external_links` Format

| Column | Required | Description |
|---|---|---|
| `peak_id` | ✅ | Peak identifier, format `chr:start-end` |
| `gene` | ✅ | Gene name |
| `tf_name` | optional | TF annotation |
| `source` | optional | Link provenance (e.g. `scenic`, `scarlink`, `archr`) |

### 5.3 `MoDEData.from_matrices()` Parameters

| Parameter | Type | Description |
|---|---|---|
| `rna_counts` | `str` or `DataFrame` | RNA count matrix (samples × genes) |
| `atac_counts` | `str` or `DataFrame` | ATAC peak count matrix (samples × peaks) |
| `metadata` | `str` or `DataFrame` | Sample metadata; must contain the condition column |
| `condition_col` | `str` | Condition column name |
| `donor_col` | `str` | Donor column name (optional) |
| `batch_col` | `str` | Batch column name (optional) |
| `index_col` | `int` | Row index column in TSV files (default 0) |

---

## 6. Input File Formats

### 6.1 Bulk Count Matrix (TSV)

**rna_counts.tsv** — rows = samples, columns = genes:

| sample | STAT1 | GZMB | IL7R |
|---|---|---|---|
| ctrl_1 | 120 | 85 | 200 |
| ctrl_2 | 98 | 92 | 180 |
| trt_1 | 350 | 90 | 220 |

**atac_counts.tsv** — rows = samples, columns = peaks:

| sample | chr1:100-200 | chr2:300-400 |
|---|---|---|
| ctrl_1 | 50 | 60 |
| trt_1 | 180 | 55 |

### 6.2 Sample Metadata

**metadata.tsv**:

| sample | condition | batch |
|---|---|---|
| ctrl_1 | control | A |
| trt_1 | treatment | B |

The `condition` column must be binary (e.g. control / treatment).

### 6.3 Peak-Gene Links

**peak_gene_links.tsv**:

| peak_id | gene | tf_name |
|---|---|---|
| chr1:100-200 | STAT1 | IRF1 |
| chr2:300-400 | GZMB | |

> **Important**: Plain gene symbols (e.g. STAT1, GZMB) lack genomic coordinates. If you
> do not provide `external_links`, `genome_annotation`, or `tss_map`, MoDES will raise
> an error because no candidate events can be generated.

---

## 7. Algorithm Details

### Step 1: Event Candidate Construction

For each gene, candidate regulatory elements are identified from:

- **Promoter peaks**: within ±2 kb of the TSS
- **Distal peaks**: within ±250 kb of the TSS
- **External links**: pre-computed links from SCENIC+, SCARlink, ArchR, etc.
- **Motif annotation**: TF binding motifs on peaks

### Step 2: Effect Size Estimation

For each peak and gene, a **Negative Binomial GLM** (log link) estimates the condition
effect:

```
log(E[Y_u]) = α + β_cond × C_u + X_u × γ + offset(log(libsize_u))
```

where:
- `C_u` = condition indicator (0/1)
- `X_u` = covariate matrix (batch, donor, age, etc.)
- `offset` = log library size (DESeq2-style median-of-ratios normalization)

**Limma-style empirical Bayes variance moderation** is applied post-hoc:
- Variance estimates are pooled across all features to estimate a prior distribution
- Posterior variances produce more stable moderated t-statistics
- Particularly effective for small-sample scenarios

Multiple testing correction: **Benjamini-Hochberg FDR** (per modality).

### Step 3: Conditional Decomposition

The core statistical question: **can the RNA condition effect be explained by local
chromatin accessibility?**

Two models are compared:

**Model 0 (RNA-only)**:
```
log(E[RNA_g]) ~ Condition + Covariates
```

**Model 1 (RNA | ATAC)**:
```
log(E[RNA_g]) ~ Condition + ATAC_peak + Covariates
```

The key quantity is the attenuation of β_cond from Model 0 to Model 1:

| β_cond Change | Interpretation |
|---|---|
| Significant in M0, not in M1 | RNA change explained by local chromatin (concordant evidence) |
| Significant in M0, still significant in M1, little attenuation | RNA has residual effect not explained by ATAC (rna_only evidence) |
| Not significant in M0 | No RNA-layer effect |

> **Note**: The conditional model uses a single linked peak as covariate. The correct
> interpretation is "condition effect after adjustment for the linked ATAC peak," not
> "ATAC explains RNA."

### Step 4: Evidence Vector Construction

For each event `e`, a standardized evidence vector is built:

```
D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality_score]
```

where:
- `z_m = coef_m / SE(coef_m)` — standardized effect (z-score)
- `quality_score ∈ [0, 1]` — data quality score (from sequencing depth, detection rate,
  batch association, etc.)

### Step 5: State Classification

Two-stage classification:

**Stage 1 — Rule-based**:
```
if ATAC_sig and RNA_sig and same_direction → concordant
elif ATAC_sig and RNA_sig and opposite → discordant_opposite
elif ATAC_sig and not RNA_sig → chromatin_primed
elif not ATAC_sig and RNA_sig → rna_only
else → null
```

**Stage 2 — Empirical Bayes refinement**:
- Rule-based labels serve as initial assignments
- Per-state Gaussian distributions are fitted in evidence-vector space
- State confidence is computed for each event (rule-refined confidence, not a
  calibrated posterior)

**Artifact risk** is computed independently:
```
quality_score < threshold
  → check if single-modality significant
    yes → artifact_risk = high
    no  → artifact_risk = medium
else → artifact_risk = low
```

---

## 8. Statistical Models

### 8.1 GLM Fitting Strategy

Each feature (peak or gene) is fit with a multi-level fallback:

| Priority | Model | Family | When Used |
|---|---|---|---|
| 1 | NB with default α | NegativeBinomial | Preferred; most flexible |
| 2 | NB with α = 1.0 | NegativeBinomial | Fallback 1: more stable |
| 3 | Poisson | Poisson | Fallback 2: no overdispersion |
| 4 | Simplified NB | NegativeBinomial (α=1.0) | Fallback 3: intercept + condition only |

All fallback outcomes are recorded in `model_diagnostics.tsv` with the specific
`model_used`, `family`, and whether covariates were dropped.

### 8.2 Design Matrix

The design matrix includes:
- Intercept
- Condition (binary, dummy-coded)
- Covariates (numeric or one-hot encoded categorical)
- Batch (one-hot encoded, first category as reference)
- Donor (one-hot encoded, first category as reference)

A rank check is performed; rank-deficient designs (e.g. condition perfectly confounded
with batch) raise a `ValueError`.

### 8.3 Event-Level Significance

Event-level p-values are computed based on state:
- `concordant` / `discordant_opposite`: `max(atac_pval, rna_pval)`
- `chromatin_primed`: `atac_pval`
- `rna_only`: `rna_pval`
- `null`: 1.0

BH correction is applied across all events to produce `event_fdr`.

---

## 9. Output Files

### 9.1 File Inventory

| File | Format | Description |
|---|---|---|
| `event_table.tsv` | TSV | **Main output**: one row per peak-gene event |
| `event_state_confidence.tsv` | TSV | State, state_confidence, artifact_risk per event |
| `event_layer_effects.tsv` | TSV | ATAC/RNA layer effect estimates |
| `event_evidence_vectors.tsv` | TSV | Evidence vectors used for state classification |
| `model_diagnostics.tsv` | TSV | Marginal ATAC/RNA GLM diagnostics |
| `run_params.tsv` | TSV | Run parameters and data summary |
| `event_network.graphml` | GraphML | TF-peak-gene network (optional) |
| `report.html` | HTML | Visualizations and statistical summary (optional) |

> `model_diagnostics.tsv` currently reports marginal ATAC and RNA GLM diagnostics only.
> Conditional RNA-after-ATAC diagnostics are summarized in `event_table.tsv`;
> full conditional diagnostics will be added in a future version.

### 9.2 event_table.tsv — Complete Field Reference

| Field | Type | Description |
|---|---|---|
| `event_id` | str | Unique event identifier |
| `tf_name` | str/NA | TF name associated with this peak (from motif annotation) |
| `peak_id` | str | ATAC peak identifier (e.g. `chr1:100-200`) |
| `gene` | str | Target gene name |
| `context` | str | Context label (e.g. cell type + condition) |
| `state` | str | Biological state: `concordant` / `chromatin_primed` / `rna_only` / `discordant_opposite` / `null` |
| `state_confidence` | float | State confidence, range 0–1 (rule-refined; not a calibrated posterior) |
| `artifact_risk` | str | Technical risk: `low` / `medium` / `high` (not a biological state) |
| `artifact_reason` | str | Risk reasons, semicolon-separated (e.g. `low_quality_score;single_modality_low_quality`) |
| `event_pval` | float | Event-level p-value (computed from state-specific combination of ATAC/RNA p-values) |
| `event_fdr` | float | Event-level BH-corrected FDR |
| `atac_coef` | float | Condition effect on ATAC peak (log fold change scale) |
| `atac_se` | float | Standard error of ATAC effect |
| `atac_pval` | float | ATAC-layer p-value |
| `atac_fdr` | float | ATAC-layer BH-corrected FDR |
| `atac_direction` | int | ATAC effect direction (+1 / -1 / 0) |
| `rna_coef` | float | Condition effect on RNA gene (log fold change scale) |
| `rna_se` | float | Standard error of RNA effect |
| `rna_pval` | float | RNA-layer p-value |
| `rna_fdr` | float | RNA-layer BH-corrected FDR |
| `rna_direction` | int | RNA effect direction (+1 / -1 / 0) |
| `rna_after_atac_coef` | float | Condition effect after adjusting for the linked ATAC peak |
| `rna_after_atac_se` | float | Standard error of conditional effect |
| `rna_after_atac_pval` | float | Conditional effect p-value |
| `rna_after_atac_fdr` | float | Conditional effect BH-corrected FDR |
| `quality_score` | float | Event quality score (0–1; based on sequencing depth, detection rate, etc.) |

---

## 10. Interpretation Guide

### 10.1 How to Read States

**concordant** — High-confidence cis-regulatory event. Local chromatin opening is
accompanied by transcriptional increase. Suitable for building core regulatory networks.

**chromatin_primed** — Chromatin has changed but transcription has not. May represent:
- Lineage priming during development
- Early epigenetic response to drug treatment
- Poised enhancers requiring a second stimulus

**rna_only** — RNA change not explained by local ATAC. May arise from:
- Trans regulation via distal enhancers
- RNA stability changes
- TF activity changes without chromatin remodeling
- Peak calling that missed key regulatory regions

**discordant_opposite** — ATAC and RNA move in opposite directions. Common explanations:
- Negative feedback regulation
- Repressor binding
- Requires further experimental validation

**null** — No significant change under the current condition for this event.

### 10.2 Credibility Assessment

1. **Check `artifact_risk` first**: prioritize `low` and `medium`; treat `high` with caution
2. **Then check `state_confidence`**: events with > 0.8 are more reliable
3. **Finally check `event_fdr`**: event-level FDR < 0.1 provides statistical guarantee

Recommended filtering pipeline:

```python
trusted = result.filter(
    min_confidence=0.7,
    max_event_fdr=0.1,
    exclude_high_artifact=True,
)
```

### 10.3 FAQ

**Q: Why do I have very few concordant events?**

A: Possible reasons:
- Small sample size (n < 5 per group) — insufficient statistical power
- Genuinely small effect sizes (biological)
- Strict FDR threshold — try increasing `fdr_threshold`
- Low-quality peak-gene links — try high-quality external links from SCENIC+ / SCARlink

**Q: Why are many events null?**

A: Null is expected — most peak-gene pairs are not differentially regulated in any given
condition. Focus on functional enrichment of non-null events.

**Q: How should I handle artifact_risk = high events?**

A: Do not discard them outright — they may still contain biological signal. Recommendations:
- Inspect the raw count data quality
- Check the `artifact_reason` field for specific causes
- Flag them as "requires validation" in downstream analyses

---

## 11. Testing

![tests](https://github.com/Lings01/MoDES/actions/workflows/tests.yml/badge.svg)

```bash
python -m pytest -q              # Run all tests
python -m pytest -k "states"     # Run state classification tests only
python -m pytest -k "integration" # Run integration tests
```

---

## 12. Citation

MoDES: Multi-Omics Discordance-guided decomposition of regulatory event states.

## 13. License

MIT
