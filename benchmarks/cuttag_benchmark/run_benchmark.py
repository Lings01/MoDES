#!/usr/bin/env python3
"""CUT&Tag benchmark: RNA + H3K27ac CUT&Tag state recovery (v2.0 multi-modal)."""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np, pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n_per_group = 15; n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    n_genes = 60
    # Genes spaced 500kb apart so 250kb promoter window never overlaps
    gene_names = [f"G{i}:chr1:{1000000+i*500000}" for i in range(n_genes)]
    atac_names = [f"chr1:{1000000+i*500000-500}-{1000000+i*500000+500}" for i in range(n_genes)]
    cuttag_names = [f"chr1:{1000000+i*500000-500}-{1000000+i*500000+500}|H3K27ac" for i in range(n_genes)]

    n_atac = len(atac_names); n_cuttag = len(cuttag_names)

    # Baseline: shared across groups so null genes have ZERO true effect
    rna_base = rng.poisson(50, (n_per_group, n_genes))
    atac_base = rng.poisson(30, (n_per_group, n_atac))
    cuttag_base = rng.poisson(30, (n_per_group, n_cuttag))
    rna = np.vstack([rna_base, rna_base])
    atac = np.vstack([atac_base, atac_base])
    cuttag = np.vstack([cuttag_base, cuttag_base])

    truth = {}
    # 10 activating: ATAC↑ + H3K27ac↑ + RNA↑
    rna[n_per_group:, :10] = rng.poisson(300, (n_per_group, 10))
    atac[n_per_group:, :10] = rng.poisson(150, (n_per_group, 10))
    cuttag[n_per_group:, :10] = rng.poisson(150, (n_per_group, 10))
    for i in range(10): truth[gene_names[i]] = "epigenomic_concordant"

    # 10 primed: H3K27ac↑, ATAC→, RNA→
    cuttag[n_per_group:, 10:20] = rng.poisson(150, (n_per_group, 10))
    for i in range(10, 20): truth[gene_names[i]] = "active_enhancer_primed"

    # 5 mark-only: H3K27ac↑, ATAC→, RNA→
    cuttag[n_per_group:, 20:25] = rng.poisson(150, (n_per_group, 5))
    for i in range(20, 25): truth[gene_names[i]] = "mark_only"

    # 5 null genes (indices 25-29) — explicitly tracked
    for i in range(25, 30): truth[gene_names[i]] = "null"
    # Remaining genes (30-59) are noise/filler — no explicit truth check

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    epi_feat = pd.DataFrame({
        "feature_id": cuttag_names, "chr": ["chr1"] * n_cuttag,
        "start": [1000000+i*500000-500 for i in range(n_cuttag)],
        "end": [1000000+i*500000+500 for i in range(n_cuttag)],
        "assay": ["CUTTAG"] * n_cuttag, "target": ["H3K27ac"] * n_cuttag,
    })

    print("=" * 60)
    print("CUT&Tag Benchmark: RNA + H3K27ac CUT&Tag (1:1 gene-peak, 500kb spaced)")
    print(f"{n_total} samples, {n_genes} genes, {n_atac} ATAC, {n_cuttag} CUT&Tag")
    print(f"GT: 10 activating + 10 primed + 5 mark-only + 5 null (30 genes tracked)")

    t0 = time.time()
    data = MoDEData.from_epigenomic_matrices(
        rna_counts=pd.DataFrame(rna, index=obs.index, columns=gene_names),
        epigenomic_counts=pd.DataFrame(cuttag, index=obs.index, columns=cuttag_names),
        epigenomic_features=epi_feat, metadata=obs, condition_col="condition",
        target="H3K27ac", atac_counts=pd.DataFrame(atac, index=obs.index, columns=atac_names),
    )
    tss = {g: (g.split(":")[0], "chr1", int(g.split(":")[2])) for g in gene_names}
    result = MoDES(data=data, condition_col="condition", tss_map=tss, fdr_threshold=0.05).run()
    t = time.time() - t0

    sc = result.event_table["state"].value_counts()
    pred = {}
    for _, r in result.event_table.iterrows():
        g = str(r["gene"])
        if g not in pred: pred[g] = r["state"]

    # Evaluate: only count tracked genes (indices 0-29)
    state_map = {
        "epigenomic_concordant": {"epigenomic_concordant", "concordant"},
        "active_enhancer_primed": {"active_enhancer_primed", "chromatin_primed"},
        "mark_only": {"mark_only", "active_enhancer_primed", "null"},
        "null": {"null"},
    }
    tracked = list(truth.keys())  # first 30 genes
    correct = sum(1 for g in tracked
                  if pred.get(g) == truth[g] or pred.get(g) in state_map.get(truth[g], set()))
    acc = correct / len(tracked)

    # Per-category breakdown
    for cat, indices in [("activating", range(10)), ("primed", range(10, 20)),
                          ("mark_only", range(20, 25)), ("null", range(25, 30))]:
        cat_correct = sum(1 for i in indices
                         if pred.get(gene_names[i]) == truth[gene_names[i]]
                         or pred.get(gene_names[i]) in state_map.get(truth[gene_names[i]], set()))
        n_cat = len(list(indices))
        print(f"  {cat:15s}: {cat_correct}/{n_cat}")

    print(f"\nOverall accuracy: {correct}/{len(tracked)} ({acc*100:.0f}%) | Runtime: {t:.1f}s")
    if acc >= 0.8: print("✅ Multi-modal classifier working correctly!")
    else: print("⚠️ Check data/classifier — accuracy below 80%")
    for s in ["epigenomic_concordant", "active_enhancer_primed", "mark_only",
              "concordant", "chromatin_primed", "rna_only", "null"]:
        c = sc.get(s, 0)
        if c > 0: print(f"  {s:25s}: {c}")
    print(f"  Events: {len(result.event_table)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)


if __name__ == "__main__":
    main()
