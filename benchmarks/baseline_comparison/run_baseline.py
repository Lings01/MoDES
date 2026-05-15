#!/usr/bin/env python
"""Baseline comparison: naive overlap vs MoDES state classification.

The naive baseline classifies events as:
  - concordant: ATAC significant AND RNA significant AND same direction
  - primed: ATAC significant only
  - rna_only: RNA significant only
  - null: neither significant

This is compared against MoDES full pipeline output on the same data.

Usage:
    python benchmarks/baseline_comparison/run_baseline.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import numpy as np
import pandas as pd

from modes import MoDEData, MoDES


def naive_classify(atac_fdr, rna_fdr, atac_dir, rna_dir, threshold=0.1):
    atac_sig = atac_fdr < threshold
    rna_sig = rna_fdr < threshold
    same_dir = atac_dir == rna_dir

    if atac_sig and rna_sig and same_dir:
        return "concordant"
    elif atac_sig and rna_sig and not same_dir:
        return "discordant_opposite"
    elif atac_sig and not rna_sig:
        return "chromatin_primed"
    elif not atac_sig and rna_sig:
        return "rna_only"
    return "null"


def main():
    rng = np.random.default_rng(42)
    n, n_events, n_filler = 10, 30, 200
    n_total = n * 2
    condition = np.array(["ctrl"] * n + ["trt"] * n)

    chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"]
    states_config = ["concordant"] * 6 + ["chromatin_primed"] * 6 + ["rna_only"] * 6 + \
                    ["discordant_opposite"] * 6 + ["null"] * 6
    gene_names, peak_names, true_states = [], [], []
    for i in range(n_events):
        ch = chromosomes[i % 10]
        pos = 500000 + (i // 10) * 2000000
        gene_names.append(f"g_{i}:{ch}:{pos}")
        peak_names.append(f"{ch}:{pos-100}-{pos+100}")
        true_states.append(states_config[i])

    all_gn = gene_names + [f"fg_{i}" for i in range(n_filler)]
    all_pn = peak_names + [f"fp_{i}" for i in range(n_filler)]

    rna = np.zeros((n_total, n_events + n_filler))
    atac = np.zeros((n_total, n_events + n_filler))
    for s in range(n_total):
        rna[s, :] = rng.poisson(300, n_events + n_filler)
        atac[s, :] = rng.poisson(200, n_events + n_filler)

    effect_map = {"concordant": (5.0, 4.0), "chromatin_primed": (5.0, 0.0),
                  "rna_only": (0.0, 4.0), "discordant_opposite": (4.0, -3.0),
                  "null": (0.0, 0.0)}
    for i in range(n_events):
        a_eff, r_eff = effect_map[true_states[i]]
        rna[n:, i] = rng.poisson(np.maximum(1, rna[n:, i] * np.exp(r_eff)))
        atac[n:, i] = rng.poisson(np.maximum(1, atac[n:, i] * np.exp(a_eff)))

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    data = MoDEData(rna=pd.DataFrame(rna, index=obs.index, columns=all_gn),
                    atac=pd.DataFrame(atac, index=obs.index, columns=all_pn), obs=obs)

    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        name, chrom = parts[0], "chr" + parts[1].split(":")[0]
        pos = int(parts[1].split(":")[1])
        tss_map[g] = (name, chrom, pos)

    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()

    # Naive baseline using MoDES effect estimates
    et = result.event_table
    naive = [naive_classify(r["atac_fdr"], r["rna_fdr"], r["atac_direction"], r["rna_direction"])
             for _, r in et.iterrows()]

    pred_map = dict(zip(et["gene"], et["state"]))
    naive_map = dict(zip(et["gene"], naive))

    modes_correct = sum(pred_map.get(g) == t for g, t in zip(gene_names, true_states))
    naive_correct = sum(naive_map.get(g) == t for g, t in zip(gene_names, true_states))

    print("Baseline Comparison Benchmark")
    print("=" * 50)
    print(f"MoDES accuracy: {modes_correct / n_events:.2f}")
    print(f"Naive  accuracy: {naive_correct / n_events:.2f}")
    print("MoDES uses conditional decomposition + empirical Bayes refinement.")
    print("Naive uses only per-modality significance + direction matching.")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)


if __name__ == "__main__":
    main()
