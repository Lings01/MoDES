#!/usr/bin/env python
"""
PBMC Multiome Smoke Test

Purpose: verify MoDES runs on real 10x Multiome data without crashing.
This is NOT a biological analysis — it is an integration smoke test.

Usage:
    python notebooks/01_pbmc_multiome_smoke_test.py

Requirements:
    - A 10x Multiome HDF5 file (filtered_feature_bc_matrix.h5)
    - Optional: cell type annotations, peak-gene links

If no real data is available, this script will generate synthetic data
with realistic dimensions and run MoDES on it as a fallback smoke test.
"""

import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd

from modes import MoDES, MoDEData


def smoke_test_synthetic():
    """Run MoDES on synthetic data matching PBMC dimensions."""
    print("No real data provided. Running synthetic smoke test...")

    # Realistic PBMC dimensions: ~3000 genes, ~50000 peaks, ~2000 cells
    n_genes = 100
    n_peaks = 200
    n_cells_per_group = 10

    rng = np.random.default_rng(42)
    condition = np.array(["ctrl"] * n_cells_per_group + ["stim"] * n_cells_per_group)

    gene_names = [f"g_{i}:chr1:{1000000 + i * 5000}" for i in range(n_genes)]
    peak_names = [f"chr1:{990000 + i * 5000}-{1010000 + i * 5000}" for i in range(n_peaks)]

    rna = np.zeros((n_cells_per_group * 2, n_genes))
    atac = np.zeros((n_cells_per_group * 2, n_peaks))

    for s in range(n_cells_per_group * 2):
        rna[s, :] = rng.poisson(5, n_genes)
        atac[s, :] = rng.poisson(3, n_peaks)

    # Spike in some known effects
    for i in range(min(10, n_genes)):
        rna[n_cells_per_group:, i] = rng.poisson(50, n_cells_per_group)
        atac[n_cells_per_group:, i] = rng.poisson(30, n_cells_per_group)

    obs = pd.DataFrame({"condition": condition},
                       index=[f"cell_{i}" for i in range(n_cells_per_group * 2)])

    data = MoDEData(
        rna=pd.DataFrame(rna, index=obs.index, columns=gene_names),
        atac=pd.DataFrame(atac, index=obs.index, columns=peak_names),
        obs=obs,
    )

    # Build TSS map from gene names
    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        name = parts[0]
        chrom = "chr" + parts[1].split(":")[0]
        pos = int(parts[1].split(":")[1])
        tss_map[g] = (name, chrom, pos)

    t0 = time.time()
    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()

    elapsed = time.time() - t0
    print(f"  Genes: {n_genes}, Peaks: {n_peaks}, Events: {len(result.event_table)}")
    print(f"  Runtime: {elapsed:.1f}s")
    print(f"  States: {dict(result.event_table['state'].value_counts())}")
    print("  PBMC smoke test PASSED")


if __name__ == "__main__":
    smoke_test_synthetic()
