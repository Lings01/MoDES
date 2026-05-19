# Statistical Model

MoDES v2.0.0 is a multi-omics evidence-scoring and event-state annotation pipeline
for candidate regulatory events. This document describes the statistical components
and their interpretation.

## 1. Scope and Claims

MoDES provides:

- Candidate regulatory event construction from peak-gene links
- Per-modality marginal effect estimation (NB GLM + empirical Bayes moderation)
- Multi-modal evidence scoring via state-rule grammar
- State support scores for ranking (NOT calibrated FDR)
- Artifact and quality diagnostics

MoDES does NOT:

- Infer causality
- Provide calibrated posterior probabilities
- Provide formal post-selection FDR unless calibration is performed
- Claim validated multi-omics regulatory inference for experimental modalities

## 2. Candidate Event Definition

A candidate regulatory event is a tuple `e = (peak, gene, context)` where:

- `peak`: a genomic region measured by ATAC-seq, CUT&Tag, etc.
- `gene`: a target gene measured by RNA-seq
- `context`: cell type, condition, or spatial context

Events are constructed from genomic proximity (promoter ±2kb, distal ±250kb)
and/or externally provided peak-gene links (SCENIC+, SCARlink, ArchR, etc.).

## 3. Long-Format Event-Modality Evidence

All per-modality effects are stored in `event_modality_evidence.tsv`:

```
event_id | modality | assay | target | feature_id | role |
coef | se | pval | fdr | direction | quality_score | model_used | converged
```

One event can have multiple rows — one per measured modality (RNA, ATAC,
CUT&Tag H3K27ac, CUT&Tag H3K27me3, protein, spatial, etc.).

## 4. Per-Modality Marginal GLMs

For each feature in each modality, a Negative Binomial GLM estimates condition
effects:

log(E[Y]) = α + β_cond × Condition + Xγ + offset(log(lib_size))

- Multi-level fallback: NB → NB(α=1.0) → Poisson
- Empirical Bayes variance moderation (Limma-style)
- Benjamini-Hochberg FDR correction per modality
- Directed evidence score: when effect direction matches expected, score = -log10(pval/2); when opposite, score = 0

## 5. StateRule Grammar

States are defined declaratively via `StateRule` objects with six evidence types:

| Type | Semantics |
|---|---|
| `required` | Must be significant in the specified direction |
| `required_absent` | Must be measured but NOT statistically significant |
| `forbidden` | Must NOT be significant (or not in a given direction) |
| `optional` | Bonus if significant; no penalty if absent |
| `missing_policy` | What happens when a modality is not available |

Each StateRule has a `state_family` (coarse grouping like "concordant") and
a specific `state` name (e.g., "concordant_activation", "concordant_repression").

## 6. Required / Absent / Forbidden / Optional / Missing Evidence

The StateClassifier scores all applicable rules simultaneously for each event:

1. **required**: All required modalities must be significant in correct direction
2. **required_absent**: All absent-required modalities must be measured but not significant
3. **forbidden**: No forbidden pattern may be present
4. **optional**: Bonus multiplier (1.0 + bonus) for each optional match
5. **missing_policy**: Penalty if required modalities are not available

If required or required_absent conditions fail, the rule score is 0 (invalid).

## 7. State Assignment Score

assignment_score = support_score × quality × conflict_penalty × missing_penalty
                   × specificity_bonus × optional_bonus × link_score

- **support_score**: sum of directed evidence scores for satisfied required evidence
- **quality**: per-event quality score from detection rate, depth, batch effects
- **link_score**: peak-gene link confidence (promoter=1.0, external=0.85, distal=distance-decay)

Higher scores indicate stronger multi-modal evidence for that state.
Scores are for ranking, not probabilistic interpretation.

## 8. State Support Score (NOT Formal FDR)

`state_support_score` is a ranking-oriented directional evidence score.
It is derived from the modalities that define each state, converted from
p-values using `-log10(pval/2)` when direction matches.

`state_support_adjusted_score` is a BH-style ranking adjustment applied
across events. **It is not a formal FDR** — the independence assumptions
of BH are not met for event-level data with correlated modalities.

For formal FDR control, users must perform calibration (permutation or
sample-shuffling) to estimate empirical false state rates.

## 9. Conditional Decompositions as Diagnostics

Multi-model conditional decompositions are provided in `conditional_effects.tsv`
as diagnostic summaries. They model the response modality after adjusting for
conditioning modalities (e.g., RNA after controlling for ATAC, or protein after
controlling for RNA).

**Conditional models are not used as primary state assignment evidence in v2.0.**
They indicate potential multi-layer dependencies but do not establish causality.

## 10. Artifact Risk

Artifact risk is assessed independently from biological state:

- **low**: Good data quality across modalities
- **medium**: Some quality concerns (borderline depth, detection)
- **high**: Single-modality signal with low quality; interpret with caution

Artifact reasons are concatenated from per-modality diagnostics:
low depth, low detection, batch association, library size outliers.

## 11. Experimental Modality Limitations

The RNA+ATAC core is the validated layer. CUT&Tag, protein, spatial, and
dynamic modality modules are experimental extensions:

- CUT&Tag states depend on interval-overlap matching of peak regions
- Protein states require explicit protein_gene_links (no fuzzy matching)
- Spatial states depend on spatial autocorrelation and neighbor effect estimates
- Dynamic/pseudotime states are conceptual skeletons

## 12. What MoDES Does Not Claim

- MoDES does not infer regulatory causality.
- MoDES does not provide calibrated posterior probabilities.
- MoDES does not provide formal post-selection FDR without user calibration.
- State assignment scores are ranking heuristics, not probabilities.
- Conditional decomposition attenuation does not prove mediation.
