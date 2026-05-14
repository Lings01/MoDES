# Frequently Asked Questions

## General

**Q: What is MoDES?**

A: MoDES (Multi-Omics Discordance/Event State inference) is a statistical framework
that treats regulatory events — not genes, peaks, or clusters — as the primary
unit of analysis. It classifies cross-modality concordance and discordance into
interpretable regulatory states.

**Q: Is MoDES a multiome integration method?**

A: No. MoDES is an event-state decomposition method. Unlike integration methods
that produce embeddings or clusters, MoDES classifies each peak-gene regulatory
event into a biological state based on cross-modality consistency patterns.

**Q: What does v0.1.0 support?**

A: RNA + ATAC, binary condition, bulk or externally-generated pseudobulk.
Protein, spatial graph, and multi-class conditions are planned for future releases.

## Usage

**Q: Can I run MoDES on cells directly?**

A: Not recommended. Use donor-aware pseudobulk aggregation
(`donor × condition × cell_type`). Treating cells as independent replicates
produces inflated significance.

**Q: How do I provide peak-gene links?**

A: Three options:
1. `external_links` DataFrame (recommended): pre-computed links from SCENIC+, SCARlink, etc.
2. `genome_annotation`: GTF/GFF file for automatic TSS-based linking
3. `tss_map`: manual gene → (name, chr, tss_pos) dictionary

**Q: My gene symbols don't have coordinates. What should I do?**

A: Use `external_links`. Plain gene symbols (STAT1, GZMB) lack genomic coordinates.
Without annotation, MoDES cannot generate candidate events.

**Q: Can I run multiple cell types at once?**

A: Yes. Include `cell_type` in the pseudobulk groupby columns. The `context` field
in the output captures the cell type information.

## Interpretation

**Q: Is state_confidence a posterior probability?**

A: No. `state_confidence` is an empirical confidence score derived from rule-based
initialization and optional empirical Bayes refinement. It should not be interpreted
as a calibrated posterior probability.

**Q: What does artifact_risk = high mean?**

A: The event shows a significant signal in only one modality and has a low quality
score. It may still contain biological signal but should be interpreted with caution.
Check the `artifact_reason` field for specific causes.

**Q: How should I filter results?**

A: Recommended filtering:
```python
result.filter(
    min_confidence=0.7,
    max_event_fdr=0.1,
    exclude_high_artifact=True,
)
```

**Q: Does MoDES infer peak-gene links?**

A: MoDES can generate proximity-based candidates (promoter ±2kb, distal ±250kb),
but external links from dedicated tools (SCENIC+, SCARlink, ArchR) are recommended
for production use.

## Performance

**Q: How many samples do I need?**

A: At least 3 per group for the NB GLM to be estimable. More than 5 per group is
recommended for reliable EB variance moderation.

**Q: How many events can MoDES handle?**

A: Effect estimation scales with the number of unique peaks and genes (not events).
Conditional decomposition scales with the number of events. For 100K events,
expect ~10-30 minutes on a modern CPU.

## Troubleshooting

**Q: I get "No candidate events were generated" error.**

A: Provide `external_links`, `genome_annotation`, or `tss_map`. Plain gene symbols
cannot be matched to peak coordinates without annotation.

**Q: I get "Design matrix is rank deficient" error.**

A: Check for confounding: condition perfectly correlated with donor or batch.
Each donor should appear in both conditions.

**Q: All my events are null.**

A: Possible causes: (1) effects are genuinely small, (2) sample size is too low,
(3) FDR threshold is too strict (try increasing `fdr_threshold`).

**Q: My concordant events have low state_confidence.**

A: This is normal for events where the conditional decomposition shows partial
attenuation. Check `rna_after_atac_coef` vs `rna_coef` to understand the pattern.
