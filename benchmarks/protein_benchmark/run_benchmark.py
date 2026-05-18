#!/usr/bin/env python3
"""Protein benchmark: RNA + ATAC + Protein state recovery (v2.0 1:1 gene-peak)."""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np, pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n_per_group = 15; n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    n_genes, n_atac, n_proteins = 60, 60, 60
    # 1:1 gene-peak spacing (500kb apart to avoid 250kb window overlap)
    gene_names = [f"G{i}:chr1:{1000000+i*500000}" for i in range(n_genes)]
    atac_names = [f"chr1:{1000000+i*500000-500}-{1000000+i*500000+500}" for i in range(n_atac)]
    prot_names = [f"PROT_{i}" for i in range(n_proteins)]

    # Shared baseline: ctrl and trt get same random baseline (zero true effect)
    rna_base = rng.poisson(50, (n_per_group, n_genes))
    atac_base = rng.poisson(30, (n_per_group, n_atac))
    protein_base = rng.poisson(100, (n_per_group, n_proteins))
    rna = np.vstack([rna_base, rna_base])
    atac = np.vstack([atac_base, atac_base])
    protein = np.vstack([protein_base, protein_base])

    truth = {}
    # Full activation (0-9): ATAC↑ RNA↑ Protein↑
    rna[n_per_group:, :10] = rng.poisson(300, (n_per_group, 10))
    atac[n_per_group:, :10] = rng.poisson(150, (n_per_group, 10))
    protein[n_per_group:, :10] = rng.poisson(500, (n_per_group, 10))
    for i in range(10): truth[gene_names[i]] = "full_activation"

    # Protein buffered (10-19): RNA↑ Protein→
    rna[n_per_group:, 10:20] = rng.poisson(300, (n_per_group, 10))
    for i in range(10, 20): truth[gene_names[i]] = "protein_buffered"

    # Protein memory (20-24): RNA→ Protein↑
    protein[n_per_group:, 20:25] = rng.poisson(500, (n_per_group, 5))
    for i in range(20, 25): truth[gene_names[i]] = "protein_memory"

    # Null (25-29) — explicitly tracked
    for i in range(25, 30): truth[gene_names[i]] = "null"

    prot_links = pd.DataFrame({
        "protein_id": [f"PROT_{i}" for i in range(30)],
        "gene": gene_names[:30],
    })

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])

    print("=" * 60)
    print("Protein Benchmark: RNA + ATAC + Protein (1:1 gene-peak, 500kb spaced)")
    print(f"{n_total} samples, {n_genes} genes, {n_atac} ATAC, {n_proteins} proteins")
    print(f"GT: 10 full_activation + 10 buffered + 5 memory + 5 null (30 tracked)")

    t0 = time.time()
    data = MoDEData.from_protein_matrices(
        rna_counts=pd.DataFrame(rna, index=obs.index, columns=gene_names),
        atac_counts=pd.DataFrame(atac, index=obs.index, columns=atac_names),
        protein_counts=pd.DataFrame(protein, index=obs.index, columns=prot_names),
        protein_gene_links=prot_links, metadata=obs, condition_col="condition",
    )
    tss = {g: (g.split(":")[0], "chr1", int(g.split(":")[2])) for g in gene_names}
    result = MoDES(data=data, condition_col="condition", tss_map=tss, fdr_threshold=0.1).run()
    t = time.time() - t0

    sc = result.event_table["state"].value_counts()

    pred = {}
    for _, r in result.event_table.iterrows():
        g = str(r["gene"])
        if g not in pred: pred[g] = r["state"]

    # Per-category evaluation
    for cat, indices in [("full_activation", range(10)), ("protein_buffered", range(10, 20)),
                          ("protein_memory", range(20, 25)), ("null", range(25, 30))]:
        cat_correct = sum(1 for i in indices if pred.get(gene_names[i]) == truth[gene_names[i]])
        n_cat = len(list(indices))
        print(f"  {cat:20s}: {cat_correct}/{n_cat}")

    tracked = list(truth.keys())
    correct = sum(1 for g in tracked if pred.get(g) == truth[g])
    acc = correct / len(tracked)
    print(f"\nOverall accuracy: {correct}/{len(tracked)} ({acc*100:.0f}%) | Runtime: {t:.1f}s")
    if acc >= 0.9: print("✅ Protein multi-modal classifier working correctly!")
    else: print("⚠️ Check data/classifier")

    for s in ["full_activation", "protein_buffered", "protein_memory",
              "concordant", "chromatin_primed", "rna_only", "null"]:
        c = sc.get(s, 0)
        if c > 0: print(f"  {s:20s}: {c}")
    print(f"  Events: {len(result.event_table)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)


if __name__ == "__main__":
    main()
