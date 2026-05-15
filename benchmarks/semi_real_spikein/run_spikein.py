#!/usr/bin/env python
"""Semi-real spike-in benchmark: real count distributions with controlled effects.

Generates count matrices using realistic library sizes and sparsity patterns,
spikes in known event-state effects, evaluates MoDES recovery.

Usage:
    python benchmarks/semi_real_spikein/run_spikein.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import numpy as np
import pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n_per_group, n_events, n_filler = 12, 20, 300
    n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"]
    states = ["concordant", "chromatin_primed", "rna_only", "discordant_opposite", "null"] * 4
    gene_names, peak_names, true_states = [], [], []
    for i in range(n_events):
        ch = chromosomes[i % 10]
        pos = 500000 + (i // 10) * 2000000
        gene_names.append(f"g_{i}:{ch}:{pos}")
        peak_names.append(f"{ch}:{pos-100}-{pos+100}")
        true_states.append(states[i])

    all_gn = gene_names + [f"fg_{i}" for i in range(n_filler)]
    all_pn = peak_names + [f"fp_{i}" for i in range(n_filler)]

    # Realistic library size distribution (log-normal)
    lib_sizes = rng.lognormal(8.0, 0.5, n_total)
    # Realistic sparsity: most genes low, few high
    gene_means = rng.exponential(50, n_events + n_filler)

    rna = np.zeros((n_total, n_events + n_filler))
    atac = np.zeros((n_total, n_events + n_filler))

    for s in range(n_total):
        expected = gene_means * lib_sizes[s] / np.exp(8.0)
        rna[s, :] = rng.poisson(np.maximum(1, expected))
        atac[s, :] = rng.poisson(np.maximum(1, expected * 0.6))

    # Spike in effects for trt samples
    effect_map = {"concordant": (4.0, 3.0), "chromatin_primed": (4.0, 0.0),
                  "rna_only": (0.0, 3.0), "discordant_opposite": (3.0, -2.0),
                  "null": (0.0, 0.0)}
    for i in range(n_events):
        a_eff, r_eff = effect_map[true_states[i]]
        rna[n_per_group:, i] = rng.poisson(np.maximum(1, rna[n_per_group:, i] * np.exp(r_eff)))
        atac[n_per_group:, i] = rng.poisson(np.maximum(1, atac[n_per_group:, i] * np.exp(a_eff)))

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

    pred_map = dict(zip(result.event_table["gene"], result.event_table["state"]))
    predicted = [pred_map.get(g, "null") for g in gene_names]
    truth = true_states

    correct = sum(p == t for p, t in zip(predicted, truth))
    acc = correct / len(truth)

    print("Semi-Real Spike-in Benchmark")
    print("=" * 50)
    print(f"Realistic library sizes + sparsity, {n_events} spiked-in events")
    print(f"Accuracy: {acc:.2f}")
    for s in ["concordant", "chromatin_primed", "rna_only", "discordant_opposite", "null"]:
        idx = [i for i, t in enumerate(truth) if t == s]
        corr = sum(predicted[i] == s for i in idx)
        print(f"  {s:20s}: {corr}/{len(idx)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)


if __name__ == "__main__":
    main()
