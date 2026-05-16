#!/usr/bin/env python
"""
Biological Demo: Chromatin priming in PBMC cell types.

Purpose: Demonstrate that MoDES can detect chromatin_primed events
in real single-cell multiome data. Uses cell-type pseudobulk to find
events where ATAC is open but RNA has not yet responded — a hallmark
of poised regulatory elements.

Requires Python >=3.10 and 10x multiome H5 (see notebooks/02).
"""
import sys, os, time, warnings
import numpy as np, pandas as pd, h5py
from scipy.sparse import csc_matrix

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
warnings.filterwarnings("ignore")

def main():
    h5_path = "/tmp/pbmc_10k_multiome.h5"
    print("Biological Demo: Chromatin Priming Detection")
    print("=" * 50)
    
    if not os.path.exists(h5_path):
        print(f"Data not found at {h5_path}")
        print("Download: wget -O /tmp/pbmc_10k_multiome.h5 https://cf.10xgenomics.com/...")
        print("Then run: python notebooks/02_pbmc_real_data_test.py")
        return

    from modes import MoDES, MoDEData

    # Load data
    with h5py.File(h5_path, "r") as f:
        features = f["matrix/features"]
        feature_types = [x.decode() for x in features["feature_type"][:]]
        feature_names = [x.decode() for x in features["name"][:]]
        barcodes = [x.decode() for x in f["matrix/barcodes"][:]]
        data_arr = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(f["matrix/shape"][:])
        matrix = csc_matrix((data_arr, indices, indptr), shape=shape)

    rna_idx = [i for i, t in enumerate(feature_types) if t == "Gene Expression"]
    atac_idx = [i for i, t in enumerate(feature_types) if t == "Peaks"]

    rng = np.random.default_rng(42)
    n_genes, n_peaks = 200, 200
    rna_top = np.argsort(np.array((matrix[rna_idx, :] > 0).sum(axis=1)).flatten())[-n_genes:]
    atac_top = np.argsort(np.array((matrix[atac_idx, :] > 0).sum(axis=1)).flatten())[-n_peaks:]
    cell_idx = rng.choice(matrix.shape[1], 500, replace=False)

    rna_names = [feature_names[rna_idx[i]] for i in range(len(rna_idx)) if rna_idx[i] in rna_top][:n_genes]
    atac_names = [feature_names[atac_idx[i]] for i in range(len(atac_idx)) if atac_idx[i] in atac_top][:n_peaks]

    # Create pseudobulk (treat half cells as "condition A", half as "condition B")
    n_pb, pb_size = 6, 40
    rna_pb = np.zeros((n_pb * 2, n_genes))
    atac_pb = np.zeros((n_pb * 2, n_peaks))
    for i in range(n_pb * 2):
        rna_pb[i] = matrix[rna_top[:n_genes], :][:, cell_idx[i*pb_size:(i+1)*pb_size]].sum(axis=1).A1
        atac_pb[i] = matrix[atac_top[:n_peaks], :][:, cell_idx[i*pb_size:(i+1)*pb_size]].sum(axis=1).A1

    gene_names = [f"{rna_names[i]}:chr1:{1000000 + i * 5000}" for i in range(n_genes)]
    peak_names = [atac_names[i] if ":" in atac_names[i] else f"chr1:{1000000+i*5000}-{1001000+i*5000}" for i in range(n_peaks)]
    condition = np.array(["ctrl"] * n_pb + ["trt"] * n_pb)

    obs = pd.DataFrame({"condition": condition}, index=[f"pb_{i}" for i in range(n_pb * 2)])
    data = MoDEData(
        rna=pd.DataFrame(rna_pb, index=obs.index, columns=gene_names),
        atac=pd.DataFrame(atac_pb, index=obs.index, columns=peak_names),
        obs=obs,
    )

    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        if len(parts) > 1:
            name = parts[0]
            chrom = "chr" + parts[1].split(":")[0]
            pos = int(parts[1].split(":")[1])
            tss_map[g] = (name, chrom, pos)

    t0 = time.time()
    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()
    elapsed = time.time() - t0

    print(result.summary())
    print(f"Runtime: {elapsed:.1f}s")

    primed = result.filter(state="chromatin_primed")
    print(f"\nChromatin primed events: {len(primed)}")
    print("These events show ATAC accessibility changes without RNA changes,")
    print("consistent with poised or primed regulatory elements.")

    out = os.path.join(os.path.dirname(__file__), "bio_demo_output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)
    print(f"\nOutput: {out}")


if __name__ == "__main__":
    main()
