"""Minimal MoDES example with real-world style gene symbols.

Usage:
    python examples/minimal_bulk/run_minimal.py

This demonstrates loading bulk RNA+ATAC data with plain gene symbols
(STAT1, GZMB, IL7R) and using external peak-gene links to build events.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import pandas as pd
from modes import MoDES, MoDEData

# Load data
base = os.path.dirname(__file__)
data = MoDEData.from_matrices(
    rna_counts=os.path.join(base, "rna_counts.tsv"),
    atac_counts=os.path.join(base, "atac_counts.tsv"),
    metadata=os.path.join(base, "metadata.tsv"),
    condition_col="condition",
    index_col=0,
)

# Load external peak-gene links
links = pd.read_csv(os.path.join(base, "peak_gene_links.tsv"), sep="\t")

modes = MoDES(
    data=data,
    condition_col="condition",
    external_links=links,
    fdr_threshold=0.5,
)

result = modes.run()
print(result.summary())

# Save outputs
out_dir = os.path.join(base, "output")
result.to_tsv(out_dir)
print(f"Output written to: {out_dir}")
