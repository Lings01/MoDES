# Single-Cell Pseudobulk Analysis

MoDES v0.1.0 supports single-cell multiome data through donor-aware pseudobulk
aggregation. **Cells are never treated as independent biological replicates.**

## Rationale

Single-cell differential analysis that treats each cell as an independent replicate
produces inflated significance and donor-level pseudoreplication artifacts. The
recommended approach is pseudobulk aggregation:

1. Group cells by `donor × condition × cell_type`
2. Sum counts within each group
3. Treat each pseudobulk sample as an independent observation

## Usage

```python
from modes.data import MoDEData

data = MoDEData.from_pseudobulk(
    adata,
    groupby=["donor", "condition", "cell_type"],
    condition_col="condition",
    donor_col="donor",
    min_cells_per_group=20,
)
```

### Input Requirements

- `adata.X`: RNA count matrix (cells × genes)
- `adata.obsm["atac"]`: ATAC count matrix (cells × peaks)
  or `adata.layers[atac_layer]`
- `adata.uns["atac_var_names"]`: ATAC feature names
- `adata.obs`: must contain `donor`, `condition`, `cell_type` columns

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `groupby` | required | Columns to group by (e.g. `["donor", "condition", "cell_type"]`) |
| `condition_col` | required | Condition column name |
| `donor_col` | required | Donor column name |
| `min_cells_per_group` | 20 | Minimum cells per pseudobulk group |
| `aggregation` | `"sum"` | Aggregation function (`sum` or `mean`) |
| `atac_layer` | `None` | ATAC layer name (if not in `obsm["atac"]`) |

### Output

- `data.rna`: pseudobulk RNA count matrix (pseudobulk_samples × genes)
- `data.atac`: pseudobulk ATAC count matrix (pseudobulk_samples × peaks)
- `data.obs`: pseudobulk metadata with `context` column (cell_type_condition)

## Important Notes

1. **Cells are not replicates.** Each pseudobulk group is one observation.
2. **Groups below min_cells_per_group are dropped** with a warning.
3. **Donor must appear in both conditions** for within-donor comparison.
   A donor appearing in only one condition creates confounding.
4. **External pseudobulk is recommended** for production use.
   Aggregate externally with established tools, then load with `from_matrices()`.

## Multi-Cell-Type Batch Analysis

For analyzing multiple cell types in one run:

```python
# Load with all cell types
data = MoDEData.from_pseudobulk(adata, groupby=["donor", "condition", "cell_type"], ...)

# Run MoDES — context column captures cell type
modes = MoDES(data=data, condition_col="condition")
result = modes.run()

# Filter by cell type
t_cell_events = result.event_table[result.event_table["context"].str.contains("T_cell")]
```

## Experimental Status

Single-cell pseudobulk support is experimental in v0.1.0. For production use:

1. Aggregate RNA and ATAC to `donor × condition × cell_type` pseudobulk externally
2. Load pseudobulk matrices with `MoDEData.from_matrices()`

Native cell-level mixed models are planned for a future release (v1.3+).
