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
        ("concordant", (2.0, 2.0)),
        ("chromatin_primed", (2.0, 0.0)),
        ("rna_only", (0.0, 2.0)),
        ("discordant_opposite", (2.0, -1.5)),
        ("null", (0.0, 0.0)),
    ]

    # Add filler to balance library sizes
    n_filler = 2000
    total_events = n_events_per_state * len(state_configs)
    gene_names = []
    peak_names = []
    true_states = []
    atac_effects = []
    rna_effects = []

    # Spread events across chromosomes to avoid cross-linking
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10"]
    event_idx = 0
    for state, (a_eff, r_eff) in state_configs:
        for i in range(n_events_per_state):
            ch = chromosomes[event_idx % len(chromosomes)]
            pos = 500000 + (event_idx // len(chromosomes)) * 1000000
            gene_names.append(f"g_{state}_{i}:{ch}:{pos}")
            peak_names.append(f"{ch}:{pos-100}-{pos+100}")
            true_states.append(state)
            atac_effects.append(a_eff)
            rna_effects.append(r_eff)
            event_idx += 1

    n_real = len(gene_names)
    all_gene_names = gene_names + [f"filler_g_{i}" for i in range(n_filler)]
    all_peak_names = peak_names + [f"filler_p_{i}" for i in range(n_filler)]

    rna = np.zeros((n_total, n_real + n_filler), dtype=float)
    atac = np.zeros((n_total, n_real + n_filler), dtype=float)

    baseline_rna_ctrl = np.concatenate([
        rng.poisson(300, n_real),
        rng.poisson(800, n_filler),
    ]).astype(float)
    baseline_atac_ctrl = np.concatenate([
        rng.poisson(200, n_real),
        rng.poisson(400, n_filler),
    ]).astype(float)

    # For null events: generate paired null values (one pair per ctrl-trt pair).
    # Each ctrl[i] and trt[i] share the same baseline for null genes,
    # ensuring zero condition effect for those features.
    null_rna_pairs = rng.poisson(baseline_rna_ctrl[:n_real], (n_per_group, n_real)).astype(float)
    null_atac_pairs = rng.poisson(baseline_atac_ctrl[:n_real], (n_per_group, n_real)).astype(float)

    for p in range(n_per_group):
        ctrl_idx = p
        trt_idx = n_per_group + p

        # Both ctrl and trt use the pre-generated pair values for null genes
        # and independent Poisson for signal genes
        for i in range(n_real):
            rna_eff = rna_effects[i]
            atac_eff = atac_effects[i]

            if rna_eff == 0.0:
                # Paired null: same value for ctrl and trt → coef=0
                val = null_rna_pairs[p, i]
                rna[ctrl_idx, i] = val
                rna[trt_idx, i] = val
            else:
                rna[ctrl_idx, i] = float(rng.poisson(baseline_rna_ctrl[i]))
                rna[trt_idx, i] = float(rng.poisson(np.exp(rna_eff) * baseline_rna_ctrl[i]))

            if atac_eff == 0.0:
                val = null_atac_pairs[p, i]
                atac[ctrl_idx, i] = val
                atac[trt_idx, i] = val
            else:
                atac[ctrl_idx, i] = float(rng.poisson(baseline_atac_ctrl[i]))
                atac[trt_idx, i] = float(rng.poisson(np.exp(atac_eff) * baseline_atac_ctrl[i]))

        # Filler genes: same baseline for both
        filler_rna = rng.poisson(baseline_rna_ctrl[n_real:])
        filler_atac = rng.poisson(baseline_atac_ctrl[n_real:])
        rna[ctrl_idx, n_real:] = filler_rna
        rna[trt_idx, n_real:] = filler_rna  # identical filler → no library size bias
        atac[ctrl_idx, n_real:] = filler_atac
        atac[trt_idx, n_real:] = filler_atac

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    rna_df = pd.DataFrame(rna, index=obs.index, columns=all_gene_names)
    atac_df = pd.DataFrame(atac, index=obs.index, columns=all_peak_names)

    tss_map = {}
    for g in gene_names:
        # Parse "g_CONC_0:chr1:500000" format
        parts = g.split(":chr")
        name = parts[0]
        rest = parts[1]
        chrom = "chr" + rest.split(":")[0]
        pos = int(rest.split(":")[1])
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
    data, gt, tss_map = generate_benchmark_data(n_per_group=8, n_events_per_state=10)
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
