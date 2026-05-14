# Input Formats

MoDES v0.1.0 accepts bulk RNA+ATAC count matrices in TSV/CSV format, with optional
single-cell AnnData input via `from_pseudobulk()`.

## Required Inputs

| Input | Format | Description |
|---|---|---|
| RNA count matrix | TSV/CSV | `samples × genes` integer count matrix |
| ATAC count matrix | TSV/CSV | `samples × peaks` integer count matrix |
| Sample metadata | TSV/CSV | Sample annotations with a binary condition column |
| Event links | TSV/CSV | Peak-to-gene candidate links (or GTF / TSS map) |

## RNA Count Matrix

Rows = samples, columns = genes. Values = integer counts.

```
sample	STAT1	GZMB	IL7R
ctrl_1	120	85	200
ctrl_2	98	92	180
trt_1	350	90	220
trt_2	380	85	190
```

### Requirements

- All values must be non-negative integers (or floats castable to integers)
- Column names must be unique gene identifiers
- Row names must match the sample metadata index
- At least 3 samples (2 groups of at least 1 each, but >=3 per group is recommended)
- No NaN or infinite values

## ATAC Count Matrix

Rows = samples, columns = peaks. Values = integer counts.

```
sample	chr1:100-200	chr2:300-400
ctrl_1	50	60
ctrl_2	45	55
trt_1	180	55
trt_2	200	50
```

### Peak Name Format

Peaks must be named in `chr:start-end` format:

```
chr1:499000-501000
chr2:300000-302000
```

Underscore (`chr1_499000_501000`) and dash (`chr1-499000-501000`) variants are accepted
with a warning. Unparseable peak names are excluded from coordinate-based event generation.

## Sample Metadata

```
sample	condition	batch
ctrl_1	control	A
ctrl_2	control	A
trt_1	treatment	B
trt_2	treatment	B
```

### Requirements

- Row names must match the RNA and ATAC matrix row names
- `condition` column must contain exactly 2 unique values (binary comparison)
- Optional: `donor`, `batch`, and numeric/ categorical covariate columns

## Event Links

MoDES requires peak-to-gene links to define candidate regulatory events. Provide one of:

| Method | Description |
|---|---|
| `external_links` DataFrame | Pre-computed peak-gene links (SCENIC+, SCARlink, etc.) |
| `genome_annotation` GTF/GFF | Automatically build promoter (±2kb) and distal (±250kb) links |
| `tss_map` dict | Manual gene → (name, chr, tss_pos) mapping |

### external_links Format

```
peak_id	gene	tf_name	source
chr1:100-200	STAT1	IRF1	scenic
chr2:300-400	GZMB		scarlink
```

Required columns: `peak_id`, `gene`.
Optional columns: `tf_name`, `source`, `score`, `distance`, `link_type`.

### When No Links Are Found

If no candidate events can be generated, MoDES raises:

```
ValueError: No candidate events were generated. For normal gene symbols,
provide one of: genome_annotation, tss_map, or external_links.
```

Plain gene symbols (e.g. STAT1, GZMB) lack genomic coordinates. Without annotation,
coordinate-based event generation cannot match peaks to genes.

## Single-Cell Pseudobulk Input (Experimental)

```python
MoDEData.from_pseudobulk(
    adata,
    groupby=["donor", "condition", "cell_type"],
    condition_col="condition",
    donor_col="donor",
    min_cells_per_group=20,
)
```

Requirements:
- RNA in `adata.X`
- ATAC in `adata.obsm["atac"]` or `adata.layers[atac_layer]`
- ATAC feature names in `adata.uns["atac_var_names"]`
- Metadata columns (`donor`, `condition`, `cell_type`) in `adata.obs`

## Validation

Use the CLI to validate inputs:

```bash
modes validate-input \
  --rna rna_counts.tsv \
  --atac atac_peaks.tsv \
  --metadata metadata.tsv \
  --condition condition \
  --external-links peak_gene_links.tsv \
  --out validation_report.txt
```

Or programmatically:

```python
data = MoDEData.from_matrices(...)
issues = data.validate()
for issue in issues:
    print(issue)
```
