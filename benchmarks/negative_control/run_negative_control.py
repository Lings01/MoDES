#!/usr/bin/env python
"""Negative-control benchmark: shuffled labels should produce mostly null events.

Usage:
    python benchmarks/negative_control/run_negative_control.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import numpy as np
import pandas as pd
from modes import MoDES, MoDEData


def main():
    rng = np.random.default_rng(42)
    n, n_filler = 10, 200
    n_total = n * 2

    gene_names = [f"g_{i}:chr1:{500000+i*1000000}" for i in range(10)]
    peak_names = [f"chr1:{499900+i*1000000}-{500100+i*1000000}" for i in range(10)]
    all_gene_names = gene_names + [f"fg_{i}" for i in range(n_filler)]
    all_peak_names = peak_names + [f"fp_{i}" for i in range(n_filler)]

    rna = np.zeros((n_total, 10 + n_filler))
    atac = np.zeros((n_total, 10 + n_filler))

    for s in range(n_total):
        rna[s, :] = rng.poisson(300, 10 + n_filler)
        atac[s, :] = rng.poisson(200, 10 + n_filler)

    # True condition
    true_cond = np.array(["ctrl"] * n + ["trt"] * n)
    # Shuffled condition (breaks any real signal)
    shuffled_cond = true_cond.copy()
    rng.shuffle(shuffled_cond)

    obs = pd.DataFrame({"condition": shuffled_cond},
                       index=[f"s{i}" for i in range(n_total)])
    data = MoDEData(
        rna=pd.DataFrame(rna, index=obs.index, columns=all_gene_names),
        atac=pd.DataFrame(atac, index=obs.index, columns=all_peak_names),
        obs=obs,
    )
    tss_map = {}
    for g in gene_names:
        parts = g.split(":chr")
        name, chrom, pos = parts[0], "chr1", int(parts[1].split(":")[1])
        tss_map[g] = (name, chrom, pos)

    modes = MoDES(data=data, condition_col="condition", tss_map=tss_map)
    result = modes.run()

    state_counts = result.event_table["state"].value_counts()
    null_frac = state_counts.get("null", 0) / len(result.event_table)

    print("Negative Control Benchmark")
    print("=" * 40)
    print(f"Condition labels: shuffled (no true signal)")
    print(f"Events: {len(result.event_table)}")
    print(f"Null fraction: {null_frac:.2f} (expected: high)")
    print(f"State distribution: {dict(state_counts)}")

    out = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out, exist_ok=True)
    result.to_tsv(out)

    # Success: most events should be null
    assert null_frac > 0.5, f"Expected >50% null, got {null_frac:.1%}"
    print("PASS: shuffled labels produce mostly null events")


if __name__ == "__main__":
    main()
