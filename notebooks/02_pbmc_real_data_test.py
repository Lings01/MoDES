#!/usr/bin/env python
"""
PBMC 10k Multiome Real Data Test

Uses real 10x multiome PBMC data (downloaded from 10x Genomics).
Creates pseudobulk aggregates and spike-in known effects to verify MoDES.

Usage:
    # First download the data:
    wget -O /tmp/pbmc_10k_multiome.h5 \
      https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

    # Then run:
    python notebooks/02_pbmc_real_data_test.py
"""

import sys
import os
import time
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix

from modes import MoDES, MoDEData


def load_pbmc_10k(h5_path):
    """Load 10x multiome H5, return RNA and ATAC DataFrames."""
    print(f"Loading {h5_path}...")
    with h5py.File(h5_path, "r") as f:
        features = f["matrix/features"]
        feature_names = [x.decode() for x in features["name"][:]]
        feature_types = [x.decode() for x in features["feature_type"][:]]
        barcodes = [x.decode() for x in f["matrix/barcodes"][:]]

        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(f["matrix/shape"][:])
        matrix = csc_matrix((data, indices, indptr), shape=shape)

    rna_idx = [i for i, t in enumerate(feature_types) if t == "Gene Expression"]
    atac_idx = [i for i, t in enumerate(feature_types) if t == "Peaks"]

    rna_names = [feature_names[i] for i in rna_idx]
    atac_names = [feature_names[i] for i in atac_idx]

    rna_mat = matrix[rna_idx, :].tocsc()
    atac_mat = matrix[atac_idx, :].tocsc()

    # Take top variable genes and peaks for a practical subset
    rng = np.random.default_rng(42)
    n_genes = 500
    n_peaks = 500
    n_cells = 1000

    rna_nnz = np.array((rna_mat > 0).sum(axis=1)).flatten()
    atac_nnz = np.array((atac_mat > 0).sum(axis=1)).flatten()

    # Select top genes by coverage and top peaks by coverage
    rna_top = np.argsort(rna_nnz)[-n_genes:]
    atac_top = np.argsort(atac_nnz)[-n_peaks:]

    # Sample cells
    cell_idx = rng.choice(rna_mat.shape[1], n_cells, replace=False)

    rna_mat = rna_mat[rna_top, :][:, cell_idx].toarray().T
    atac_mat = atac_mat[atac_top, :][:, cell_idx].toarray().T
    rna_names = [rna_names[i] for i in rna_top]
    atac_names = [atac_names[i] for i in atac_top]
    barcodes = [barcodes[i] for i in cell_idx]

    print(f"  Cells (subset): {len(barcodes)}")
    print(f"  RNA genes (top expressed): {len(rna_names)}")
    print(f"  ATAC peaks (top accessible): {len(atac_names)}")

    return rna_mat, atac_mat, rna_names, atac_names, barcodes


def create_pseudobulk_test(rna, atac, rna_names, atac_names, n_cells_per_group=400):
    """Create pseudobulk with known spike-in effects."""
    rng = np.random.default_rng(42)
    n_total = n_cells_per_group * 2

    # Sample cells
    idx = rng.choice(rna.shape[0], n_total, replace=False)
    rna_sub = rna[idx, :].astype(float)
    atac_sub = atac[idx, :].astype(float)

    # Aggregate to pseudobulk (50 cells per pseudobulk sample)
    pb_size = 50
    n_pb = n_cells_per_group // pb_size
    n_pb_total = n_pb * 2

    rna_pb = np.zeros((n_pb_total, rna_sub.shape[1]))
    atac_pb = np.zeros((n_pb_total, atac_sub.shape[1]))
    for i in range(n_pb_total):
        start = i * pb_size
        end = start + pb_size
        rna_pb[i, :] = rna_sub[start:end, :].sum(axis=0)
        atac_pb[i, :] = atac_sub[start:end, :].sum(axis=0)

    condition = np.array(["ctrl"] * n_pb + ["trt"] * n_pb)

    # Spike in effects on 20 genes (and corresponding peaks)
    n_spike = 20
    spike_genes = rng.choice(len(rna_names), n_spike, replace=False)
    spike_peaks = rng.choice(len(atac_names), min(n_spike, len(atac_names)), replace=False)

    for i, (g_idx, p_idx) in enumerate(zip(spike_genes[:10], spike_peaks[:10])):
        fold = rng.choice([3, 4, 5])
        atac_pb[n_pb:, p_idx] = atac_pb[n_pb:, p_idx] * fold
        rna_pb[n_pb:, g_idx] = rna_pb[n_pb:, g_idx] * fold

    # Build gene and peak names with coordinates
    gene_names = [f"{rna_names[i]}:chr1:{1000000 + i * 5000}" for i in range(len(rna_names))]
    peak_names = list(atac_names)

    # Create external links for spiked-in genes
    links = []
    for g_idx, p_idx in zip(spike_genes[:10], spike_peaks[:10]):
        links.append({
            "peak_id": peak_names[p_idx],
            "gene": gene_names[g_idx],
            "tf_name": "",
            "source": "spikein",
        })

    obs = pd.DataFrame({"condition": condition}, index=[f"pb_{i}" for i in range(n_pb_total)])
    data = MoDEData(
        rna=pd.DataFrame(rna_pb, index=obs.index, columns=gene_names),
        atac=pd.DataFrame(atac_pb, index=obs.index, columns=peak_names),
        obs=obs,
    )

    return data, pd.DataFrame(links), spike_genes, spike_peaks


def main():
    h5_path = "/tmp/pbmc_10k_multiome.h5"
    if not os.path.exists(h5_path):
        print(f"Data not found at {h5_path}")
        print("Download with:")
        print("  wget -O /tmp/pbmc_10k_multiome.h5 \\")
        print("    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/")
        print("    pbmc_granulocyte_sorted_10k/")
        print("    pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
        return

    warnings.filterwarnings("ignore")
    rna, atac, rna_names, atac_names, barcodes = load_pbmc_10k(h5_path)

    print("\nCreating pseudobulk test data with spike-in effects...")
    data, links, spike_genes, spike_peaks = create_pseudobulk_test(
        rna, atac, rna_names, atac_names
    )
    print(f"  Pseudobulk samples: {data.n_samples}")
    print(f"  Genes: {data.n_genes}, Peaks: {data.n_peaks}")
    print(f"  External links: {len(links)}")

    print("\nRunning MoDES on real PBMC data...")
    t0 = time.time()
    modes = MoDES(
        data=data,
        condition_col="condition",
        external_links=links,
        fdr_threshold=0.2,
    )
    result = modes.run()
    elapsed = time.time() - t0

    print(f"\nResults (runtime: {elapsed:.1f}s):")
    print(result.summary())

    # Check: at least some spiked-in events should be non-null
    non_null = result.event_table[result.event_table["state"] != "null"]
    print(f"\nNon-null events: {len(non_null)}")

    out = os.path.join(os.path.dirname(__file__), "..", "notebooks", "pbmc_output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)
    print(f"Output: {out}")


if __name__ == "__main__":
    main()
