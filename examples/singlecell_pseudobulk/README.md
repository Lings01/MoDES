# Single-Cell Pseudobulk Example

Demonstrates pseudobulk aggregation of single-cell multiome data for MoDES.

## Requirements

- `anndata` installed
- Single-cell multiome AnnData with RNA in `X` and ATAC in `obsm["atac"]`

## Run

```bash
python run_pseudobulk.py
```

## Notes

- Cells are aggregated by donor × condition × cell_type
- Minimum 20 cells per group (configurable)
- Each pseudobulk sample is treated as one independent observation
- Never treat individual cells as independent replicates
