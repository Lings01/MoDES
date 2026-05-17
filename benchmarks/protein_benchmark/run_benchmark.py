#!/usr/bin/env python3
"""Protein benchmark: RNA + ATAC + Protein state recovery."""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np, pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n_per_group = 10
    n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    n_genes, n_atac, n_proteins = 100, 50, 30
    gene_names = [f"G{i}:chr1:{1000000+i*10000}" for i in range(n_genes)]
    atac_names = [f"chr1:{990000+i*20000}-{1010000+i*20000}" for i in range(n_atac)]
    prot_names = [f"PROT_{i}" for i in range(n_proteins)]

    rna = np.zeros((n_total, n_genes))
    atac = np.zeros((n_total, n_atac))
    protein = np.zeros((n_total, n_proteins))

    for s in range(n_total):
        rna[s] = rng.poisson(200, n_genes)
        atac[s] = rng.poisson(100, n_atac)
        protein[s] = rng.poisson(150, n_proteins)

    truth = {}
    # Full activation (0-9): ATAC↑ RNA↑ Protein↑
    rna[n_per_group:, :10] = rng.poisson(1000, (n_per_group, 10))
    atac[n_per_group:, :10] = rng.poisson(800, (n_per_group, 10))
    protein[n_per_group:, :10] = rng.poisson(1200, (n_per_group, 10))
    for i in range(10): truth[gene_names[i]] = "full_activation"

    # Protein buffered (10-19): RNA↑ Protein→
    rna[n_per_group:, 10:20] = rng.poisson(1000, (n_per_group, 10))
    for i in range(10, 20): truth[gene_names[i]] = "protein_buffered"

    # Protein memory (20-24): RNA→ Protein↑
    protein[n_per_group:, 20:25] = rng.poisson(1200, (n_per_group, 5))
    for i in range(20, 25): truth[gene_names[i]] = "protein_memory"

    # Null
    for i in range(25, n_genes): truth[gene_names[i]] = "null"

    prot_links = pd.DataFrame({
        "protein_id": prot_names[:25],
        "gene": gene_names[:25],
    })

    n_filler = 500
    rna_full = np.column_stack([rna, rng.poisson(300, (n_total, n_filler))])
    atac_full = np.column_stack([atac, rng.poisson(200, (n_total, n_filler))])
    gn_full = gene_names + [f"FG_{i}" for i in range(n_filler)]
    pn_full = atac_names + [f"FP_{i}" for i in range(n_filler)]

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])

    print("=" * 60)
    print("Protein Benchmark: RNA + ATAC + Protein")
    print(f"{n_total} samples, {n_genes} genes, {n_atac} ATAC peaks, {n_proteins} proteins")
    print(f"Ground truth: 10 full_activation + 10 buffered + 5 memory + 75 null")

    t0 = time.time()
    data = MoDEData.from_protein_matrices(
        rna_counts=pd.DataFrame(rna_full, index=obs.index, columns=gn_full),
        atac_counts=pd.DataFrame(atac_full, index=obs.index, columns=pn_full),
        protein_counts=pd.DataFrame(protein, index=obs.index, columns=prot_names),
        protein_gene_links=prot_links, metadata=obs, condition_col="condition",
    )
    tss = {g: (g.split(":")[0], "chr1", int(g.split(":")[2])) for g in gene_names}
    result = MoDES(data=data, condition_col="condition", tss_map=tss).run()
    t = time.time() - t0

    # Evaluate (v2.0: multi-modal protein states available)
    sc = result.event_table["state"].value_counts()
    print(f"\nRuntime: {t:.1f}s")
    print(f"Events: {len(result.event_table)}")
    for s in ["full_activation", "protein_buffered", "protein_memory",
              "concordant", "chromatin_primed", "rna_only", "null"]:
        c = sc.get(s, 0)
        if c > 0: print(f"  {s:20s}: {c}")

    # Check: ground-truth genes should be classified correctly
    pred = {}
    for _, r in result.event_table.iterrows():
        g = str(r["gene"])
        if g not in pred: pred[g] = r["state"]

    full_act_correct = sum(1 for i in range(10)
                          if pred.get(gene_names[i]) == "full_activation")
    buffered_correct = sum(1 for i in range(10, 20)
                          if pred.get(gene_names[i]) == "protein_buffered")
    memory_correct = sum(1 for i in range(20, 25)
                        if pred.get(gene_names[i]) == "protein_memory")
    print(f"\nFull activation → full_activation: {full_act_correct}/10 detected")
    print(f"Protein buffered → protein_buffered: {buffered_correct}/10 detected")
    print(f"Protein memory → protein_memory: {memory_correct}/5 detected")

    # Check modality evidence output
    if result.modality_evidence is not None and len(result.modality_evidence) > 0:
        mods = result.modality_evidence["modality"].unique()
        print(f"Modality evidence layers: {list(mods)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)
    print(f"\nOutput: {out}")


if __name__ == "__main__":
    main()
