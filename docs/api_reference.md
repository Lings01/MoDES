# API Reference

## Public API (v1.0.0-rc.1)

```python
from modes import MoDES, MoDEData, MoDESResult, EventCandidateBuilder
```

### MoDES

```python
modes = MoDES(
    data: MoDEData,
    condition_col: str,
    covariate_cols: list[str] | None = None,
    donor_col: str | None = None,
    batch_col: str | None = None,
    fdr_threshold: float = 0.1,
    genome_annotation: str | None = None,
    external_links: pd.DataFrame | None = None,
    motif_annotation: pd.DataFrame | None = None,
    tss_map: dict | None = None,
    contrast: tuple | None = None,
    allow_poisson_fallback: bool = True,
    allow_simplified_fallback: bool = False,
    conditional_mode: str = "auto",
)

result: MoDESResult = modes.run()

# Step-by-step access
events = modes.build_events()
atac_effects, rna_effects = modes.estimate_effects()
conditional = modes.decompose()
evidence = modes.build_evidence()
states = modes.classify_states()
result = modes._assemble_results()
```

### MoDEData

```python
# Bulk input
data = MoDEData.from_matrices(
    rna_counts, atac_counts, metadata,
    condition_col, donor_col=None, batch_col=None, index_col=0,
)

# AnnData input
data = MoDEData.from_anndata(adata, condition_col, donor_col=None, batch_col=None)

# Pseudobulk aggregation
data = MoDEData.from_pseudobulk(
    adata, groupby, condition_col, donor_col,
    min_cells_per_group=20, atac_layer=None,
)

# MuData input
data = MoDEData.from_mudata(
    mdata, rna_mod="rna", atac_mod="atac",
    groupby=None, condition_col="condition", donor_col=None,
)

# Spatial region-pseudobulk
data = MoDEData.from_spatial_pseudobulk(
    rna_counts, atac_counts, metadata,
    region_col="region", sample_col="sample", condition_col="condition",
)
```

### MoDESResult

```python
result.summary()                    # str summary
result.filter(state=, states=, genes=, peaks=, context=,
              min_confidence=, max_event_fdr=,
              exclude_high_artifact=, max_artifact_risk=)  # filter
result.to_tsv("output/")            # write TSV files
result.to_graphml("network.graphml") # Cytoscape/Gephi network
result.to_report("report.html")      # HTML report
result.save("output/")               # alias for to_tsv + JSON params
MoDESResult.load("output/")          # reload from directory
```

### EventCandidateBuilder

```python
builder = EventCandidateBuilder(promoter_window=2000, distal_window=250000)
events = builder.build(gene_names, peak_names, external_links=None,
                       motif_annotation=None, genome_annotation=None, tss_map=None)
```

### CLI

```bash
modes run --rna rna.tsv --atac atac.tsv --metadata meta.tsv \
  --condition condition --external-links links.tsv --out output/
modes validate-input --rna rna.tsv --atac atac.tsv \
  --metadata meta.tsv --condition condition --out validation.json
```

## API Stability
From v1.0.0-rc.1 onward, MoDES follows semantic versioning.
Breaking API changes require a major version bump.
This public API is **frozen** as of v1.0.0-rc.1.

