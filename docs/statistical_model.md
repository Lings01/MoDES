# Statistical Model

MoDES-RA v0.1.0 implements a 5-step statistical pipeline for RNA+ATAC
regulatory event state inference.

## Step 1: Event Candidate Construction

For each gene, candidate regulatory elements are identified from:

- Promoter peaks: ±2 kb of TSS
- Distal peaks: ±250 kb of TSS
- External links: pre-computed from SCENIC+, SCARlink, ArchR
- Motif annotation: TF binding motifs on peaks

## Step 2: Marginal Effect Estimation

Per-feature Negative Binomial GLM (log link):

```
log(E[Y_u]) = α + β × C_u + X_u·γ + offset(log(L_u))
```

where:
- `Y_u` = count for sample u
- `C_u` = condition indicator (0/1)
- `X_u` = covariate vector (batch, donor, etc.)
- `L_u` = library size (median-of-ratios normalization)

### GLM Fallback Strategy

Features are fit with a multi-level fallback:

| Priority | Family | α | When Used |
|---|---|---|---|
| 1 | NegativeBinomial | default | Most flexible; estimated by statsmodels |
| 2 | NegativeBinomial | 1.0 | Fallback: fixed α for stability |
| 3 | Poisson | — | Fallback: no overdispersion |
| 4 | NegativeBinomial (simplified) | 1.0 | Final fallback: intercept + condition only |

All outcomes are recorded in `model_diagnostics.tsv`.

### Empirical Bayes Variance Moderation

Limma-style moderation pools variance estimates across features:

```
σ²_posterior = (df₀·s₀² + n·s²) / (df₀ + n)
```

This stabilizes standard errors, especially for small sample sizes.

### Multiple Testing Correction

Benjamini-Hochberg FDR is applied per modality (ATAC and RNA separately).

## Step 3: Conditional Decomposition

To determine whether the RNA condition effect is explained by local chromatin,
two models are compared:

**Model 0** (RNA only):
```
log(E[RNA_g]) = α₀ + δ_RNA × C + X·γ₀
```

**Model 1** (RNA adjusted for ATAC):
```
log(E[RNA_g]) = α₁ + ρ × C + η × ATAC_peak + X·γ₁
```

Interpretation:
- `δ_RNA` significant, `ρ` not significant → RNA effect explained by local chromatin
- `δ_RNA` significant, `ρ` still significant → RNA has residual effect not explained by ATAC
- `δ_RNA` not significant → no RNA-layer effect

> The conditional model uses a single linked ATAC peak as covariate.
> Interpretation: "condition effect after adjustment for the linked ATAC peak,"
> not "ATAC explains RNA."

## Step 4: Evidence Vector

```
D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality]
```

where `z = coef / SE` and `quality ∈ [0, 1]`.

## Step 5: State Classification

### Stage 1 — Rule-Based

```
concordant:        ATAC_sig ∧ RNA_sig ∧ same_direction
chromatin_primed:  ATAC_sig ∧ ¬RNA_sig
rna_only:          ¬ATAC_sig ∧ RNA_sig
discordant_opposite: ATAC_sig ∧ RNA_sig ∧ ¬same_direction
null:              otherwise
```

### Stage 2 — Empirical Bayes Refinement

Per-state Gaussian distributions are fitted to evidence vectors, and
`state_confidence` is computed for each event. This is a rule-refined
confidence score, not a calibrated posterior probability.

### Artifact Risk

Independent of biological state:
- `quality < threshold` and single-modality significant → `high`
- `quality < threshold` → `medium`
- otherwise → `low`

### Event-Level Significance

- `concordant` / `discordant_opposite`: `event_pval = max(atac_pval, rna_pval)`
- `chromatin_primed`: `event_pval = atac_pval`
- `rna_only`: `event_pval = rna_pval`
- `null`: `event_pval = 1.0`

BH correction is applied across all events.

## Design Matrix

The design matrix `X` contains:
- Intercept
- Condition variable (binary, dummy-coded: 0 = reference, 1 = target)
- Optional covariates (numeric or one-hot encoded categorical)
- Optional batch (one-hot, first category as reference)
- Optional donor (one-hot, first category as reference)

A rank check ensures the design is not rank-deficient (e.g., condition perfectly
confounded with batch raises a `ValueError`).

## Limitations

- Binary condition only (two-group comparison)
- Fixed-effect donor/batch (no random effects)
- Single linked-peak conditional model (not cis-ATAC score aggregation)
- `state_confidence` is not a calibrated posterior
- No multi-condition contrast matrix
- No continuous pseudotime / time-lag modeling
