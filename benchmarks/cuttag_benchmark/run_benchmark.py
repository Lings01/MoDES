#!/usr/bin/env python3
"""CUT&Tag benchmark: RNA + H3K27ac CUT&Tag state recovery."""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np, pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n_per_group = 10
    n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    # 100 genes on chr1, 100 CUT&Tag peaks on chr1, 50 ATAC peaks
    n_genes, n_cuttag, n_atac = 100, 80, 50
    gene_names = [f"G{i}:chr1:{1000000+i*10000}" for i in range(n_genes)]
    cuttag_names = [f"chr1:{990000+i*10000}-{1010000+i*10000}|H3K27ac" for i in range(n_cuttag)]
    atac_names = [f"chr1:{990000+i*20000}-{1010000+i*20000}" for i in range(n_atac)]

    # Baseline counts
    rna = np.zeros((n_total, n_genes))
    cuttag = np.zeros((n_total, n_cuttag))
    atac = np.zeros((n_total, n_atac))

    for s in range(n_total):
        rna[s] = rng.poisson(200, n_genes)
        cuttag[s] = rng.poisson(100, n_cuttag)
        atac[s] = rng.poisson(80, n_atac)

    # Ground truth: 20 activating (H3K27ac↑ RNA↑), 20 primed (H3K27ac↑ RNA→),
    # 20 mark-only (H3K27ac↑ ATAC→ RNA→), 20 null
    truth = {}
    # Activating: H3K27ac↑ RNA↑
    rna[n_per_group:, :20] = rng.poisson(1000, (n_per_group, 20))
    cuttag[n_per_group:, :20] = rng.poisson(600, (n_per_group, 20))
    for i in range(20): truth[gene_names[i]] = "epigenomic_concordant"

    # Primed: H3K27ac↑ RNA null (paired)
    cuttag[n_per_group:, 20:40] = rng.poisson(600, (n_per_group, 20))
    for i in range(20, 40): truth[gene_names[i]] = "active_enhancer_primed"

    # Mark-only: H3K27ac↑, ATAC→, RNA→
    cuttag[n_per_group:, 40:60] = rng.poisson(600, (n_per_group, 20))
    for i in range(40, 60): truth[gene_names[i]] = "mark_only"

    # Null
    for i in range(60, n_genes): truth[gene_names[i]] = "null"

    # Filler to balance library sizes
    n_filler = 500
    rna_full = np.column_stack([rna, rng.poisson(300, (n_total, n_filler))])
    atac_full = np.column_stack([atac, rng.poisson(200, (n_total, n_filler))])

    gn_full = gene_names + [f"FG_{i}" for i in range(n_filler)]
    pn_full = atac_names + [f"FP_{i}" for i in range(n_filler)]
    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    epi_feat = pd.DataFrame({
        "feature_id": cuttag_names,
        "chr": ["chr1"] * n_cuttag,
        "start": [990000 + i * 10000 for i in range(n_cuttag)],
        "end": [1010000 + i * 10000 for i in range(n_cuttag)],
        "assay": ["CUTTAG"] * n_cuttag,
        "target": ["H3K27ac"] * n_cuttag,
    })

    print("=" * 60)
    print("CUT&Tag Benchmark: RNA + H3K27ac CUT&Tag")
    print(f"{n_total} samples, {n_genes} genes, {n_cuttag} CUT&Tag peaks")
    print(f"Ground truth: 20 activating + 20 primed + 20 mark-only + {n_genes-60} null")

    t0 = time.time()
    data = MoDEData.from_epigenomic_matrices(
        rna_counts=pd.DataFrame(rna_full, index=obs.index, columns=gn_full),
        epigenomic_counts=pd.DataFrame(cuttag, index=obs.index, columns=cuttag_names),
        epigenomic_features=epi_feat,
        metadata=obs, condition_col="condition",
        target="H3K27ac", atac_counts=pd.DataFrame(atac_full, index=obs.index, columns=pn_full),
    )
    tss = {g: (g.split(":")[0], "chr1", int(g.split(":")[2])) for g in gene_names}
    result = MoDES(data=data, condition_col="condition", tss_map=tss).run()
    t = time.time() - t0

    # Evaluate
    pred = {}
    for _, r in result.event_table.iterrows():
        g = str(r["gene"]).split(":")[0]
        if g not in pred: pred[g] = r["state"]

    correct = sum(1 for g, t in truth.items() if pred.get(g) == t or (t in ("epigenomic_concordant",) and pred.get(g) == "concordant"))
    print(f"\nAccuracy: {correct}/{len(truth)} ({correct/len(truth)*100:.0f}%)")
    print(f"Runtime: {t:.1f}s")
    sc = result.event_table["state"].value_counts()
    for s in ["concordant", "chromatin_primed", "rna_only", "null"]:
        print(f"  {s:20s}: {sc.get(s, 0)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)
    print(f"\nOutput: {out}")


if __name__ == "__main__":
    main()
