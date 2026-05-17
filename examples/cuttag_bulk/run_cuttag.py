#!/usr/bin/env python
"""MoDES-Epi example: RNA + CUT&Tag (H3K27ac) multi-modal analysis.

Usage: python examples/cuttag_bulk/run_cuttag.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import numpy as np
import pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n = 12
    condition = np.array(["ctrl"] * (n // 2) + ["trt"] * (n // 2))
    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n)])

    # RNA: 10 genes on chr1
    gene_names = [f"g{i}:chr1:{1000000+i*500000}" for i in range(10)]
    rna = pd.DataFrame(rng.poisson(200, (n, 10)), index=obs.index, columns=gene_names)

    # ATAC: 5 peaks
    peak_names = [f"chr1:{990000+i*500000}-{1010000+i*500000}" for i in range(5)]
    atac = pd.DataFrame(rng.poisson(100, (n, 5)), index=obs.index, columns=peak_names)

    # H3K27ac CUT&Tag: 3 regions
    epi_names = [f"chr1:{990000+i*500000}-{1010000+i*500000}|H3K27ac" for i in range(3)]
    epi = pd.DataFrame(rng.poisson(60, (n, 3)), index=obs.index, columns=epi_names)

    # Spike H3K27ac effect in treatment
    epi.iloc[n // 2:, 0] = rng.poisson(300, n // 2)
    rna.iloc[n // 2:, 0] = rng.poisson(500, n // 2)

    epi_feat = pd.DataFrame({
        "feature_id": epi_names,
        "chr": ["chr1"] * 3, "start": [990000 + i * 500000 for i in range(3)],
        "end": [1010000 + i * 500000 for i in range(3)],
        "assay": ["CUTTAG"] * 3, "target": ["H3K27ac"] * 3,
    })

    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        tss_map[g] = (parts[0], "chr" + parts[1].split(":")[0],
                      int(parts[1].split(":")[1]))

    data = MoDEData.from_epigenomic_matrices(
        rna_counts=rna, epigenomic_counts=epi,
        epigenomic_features=epi_feat, metadata=obs,
        condition_col="condition", target="H3K27ac", atac_counts=atac,
    )

    print(f"Loaded: {data.n_samples} samples, {data.n_genes} genes, {data.n_peaks} peaks")
    print(f"Modalities: {list(data.modalities.keys())}")
    print(f"CUT&Tag target: {data.modality_specs['h3k27ac_cuttag'].target}")
    print(f"Regulatory role: {data.modality_specs['h3k27ac_cuttag'].regulatory_role}")

    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()
    print(result.summary())


if __name__ == "__main__":
    main()
