#!/usr/bin/env python
"""Synthetic benchmark: evaluate MoDES state recovery against ground truth.

Generates simulated RNA+ATAC data with known event states, runs MoDES,
and computes per-state precision, recall, F1, and confusion matrix.

Usage:
    python benchmarks/simulated_event_states/run_benchmark.py
"""

import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import numpy as np
import pandas as pd

from modes import MoDES, MoDEData


def generate_benchmark_data(
    n_per_group=10,
    n_events_per_state=50,
    seed=42,
):
    """Generate data with known ground truth states."""
    rng = np.random.default_rng(seed)
    n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)

    # Define states and their effect patterns
    state_configs = [
        ("concordant", (3.0, 2.5)),
        ("chromatin_primed", (3.0, 0.0)),
        ("rna_only", (0.0, 2.5)),
        ("discordant_opposite", (3.0, -1.5)),
        ("null", (0.0, 0.0)),
    ]

    # Add filler to balance library sizes
    n_filler = 200
    total_events = n_events_per_state * len(state_configs)
    gene_names = []
    peak_names = []
    true_states = []
    atac_effects = []
    rna_effects = []

    for state, (a_eff, r_eff) in state_configs:
        for i in range(n_events_per_state):
            gene_names.append(f"g_{state}_{i}:chr1:{1000+i*2000}")
            peak_names.append(f"chr1:{900+i*2000}-{1100+i*2000}")
            true_states.append(state)
            atac_effects.append(a_eff)
            rna_effects.append(r_eff)

    n_real = len(gene_names)
    all_gene_names = gene_names + [f"filler_g_{i}" for i in range(n_filler)]
    all_peak_names = peak_names + [f"filler_p_{i}" for i in range(n_filler)]

    rna = np.zeros((n_total, n_real + n_filler))
    atac = np.zeros((n_total, n_real + n_filler))

    baseline_rna = rng.poisson(300, n_real + n_filler)
    baseline_atac = rng.poisson(200, n_real + n_filler)

    for s in range(n_total):
        is_trt = condition[s] == "trt"
        rna[s, :] = rng.poisson(baseline_rna)
        atac[s, :] = rng.poisson(baseline_atac)

    for i in range(n_real):
        if condition[0] == "ctrl":  # trt indices: n_per_group to end
            atac[n_per_group:, i] = rng.poisson(
                np.exp(atac_effects[i]) * baseline_atac[i], n_per_group
            ).astype(float)
            rna[n_per_group:, i] = rng.poisson(
                np.exp(rna_effects[i]) * baseline_rna[i], n_per_group
            ).astype(float)

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    rna_df = pd.DataFrame(rna, index=obs.index, columns=all_gene_names)
    atac_df = pd.DataFrame(atac, index=obs.index, columns=all_peak_names)

    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        name = parts[0]
        coord_parts = parts[1].split(":")
        chrom = "chr" + coord_parts[0]
        pos = int(coord_parts[1])
        tss_map[g] = (name, chrom, pos)

    gt = pd.DataFrame({
        "gene": gene_names,
        "peak_id": peak_names,
        "true_state": true_states,
    })

    return MoDEData(rna=rna_df, atac=atac_df, obs=obs), gt, tss_map


def evaluate(truth, predicted):
    """Compute per-state metrics."""
    states = ["concordant", "chromatin_primed", "rna_only",
              "discordant_opposite", "null"]
    metrics = []
    for s in states:
        tp = ((truth == s) & (predicted == s)).sum()
        fp = ((truth != s) & (predicted == s)).sum()
        fn = ((truth == s) & (predicted != s)).sum()
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        metrics.append({
            "state": s, "n_true": int((truth == s).sum()),
            "precision": round(precision, 3), "recall": round(recall, 3),
            "f1": round(f1, 3),
        })

    accuracy = (truth == predicted).mean()
    return pd.DataFrame(metrics), round(accuracy, 3)


def main():
    print("=" * 60)
    print("MoDES Benchmark: Simulated Event State Recovery")
    print("=" * 60)

    t0 = time.time()
    print("\n[1/4] Generating benchmark data...")
    data, gt, tss_map = generate_benchmark_data(n_per_group=15, n_events_per_state=40)
    print(f"      {len(data.gene_names)} genes, {data.n_samples} samples, {len(gt)} events")

    print("\n[2/4] Running MoDES...")
    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()

    print("\n[3/4] Evaluating...")
    pred_map = dict(zip(result.event_table["gene"], result.event_table["state"]))
    predicted = np.array([pred_map.get(g, "null") for g in gt["gene"]])
    truth = gt["true_state"].values

    metrics_df, accuracy = evaluate(truth, predicted)

    print(f"\n[4/4] Results:")
    print(f"      Overall accuracy: {accuracy}")
    print(f"      Runtime: {time.time() - t0:.1f}s")
    print()
    print(metrics_df.to_string(index=False))

    out_dir = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out_dir, exist_ok=True)
    gt.to_csv(os.path.join(out_dir, "truth.tsv"), sep="\t", index=False)
    result.to_tsv(os.path.join(out_dir, "moDES_output"))
    metrics_df.to_csv(os.path.join(out_dir, "metrics.tsv"), sep="\t", index=False)
    print(f"\nOutputs saved to: {out_dir}")


if __name__ == "__main__":
    main()
