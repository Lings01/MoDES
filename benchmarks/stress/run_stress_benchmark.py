#!/usr/bin/env python3
"""MoDES v2.0 Stress Benchmarks: null calibration, link-noise, batch-confounded, ablation.

Tests the robustness of the grammar-based multi-modal classifier under adversarial
conditions. Outputs benchmark_metrics.tsv with per-test metrics.
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np, pandas as pd
from modes import MoDES, MoDEData


def _build_data(rng, n_per_group=15, n_genes=40, with_effect=True):
    """Build RNA+ATAC test data with 1:1 gene-peak pairing at 500kb spacing."""
    n_total = n_per_group * 2
    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)
    gene_names = [f"G{i}:chr1:{1000000+i*500000}" for i in range(n_genes)]
    peak_names = [f"chr1:{1000000+i*500000-500}-{1000000+i*500000+500}" for i in range(n_genes)]

    rna_base = rng.poisson(50, (n_per_group, n_genes))
    atac_base = rng.poisson(30, (n_per_group, n_genes))
    rna = np.vstack([rna_base, rna_base]).astype(float)
    atac = np.vstack([atac_base, atac_base]).astype(float)

    if with_effect:
        # 10 concordant (ATAC↑+RNA↑), 10 chromatin_primed (ATAC↑), 20 null
        atac[n_per_group:, :20] = rng.poisson(200, (n_per_group, 20))
        rna[n_per_group:, :10] = rng.poisson(300, (n_per_group, 10))

    obs = pd.DataFrame({"condition": condition}, index=[f"s{i}" for i in range(n_total)])
    tss = {g: (g.split(":")[0], "chr1", int(g.split(":")[2])) for g in gene_names}
    return (pd.DataFrame(rna, index=obs.index, columns=gene_names),
            pd.DataFrame(atac, index=obs.index, columns=peak_names),
            obs, tss, gene_names)


def run_moDES(rna, atac, obs, tss, links=None, fdr=0.1):
    data = MoDEData(rna=rna, atac=atac, obs=obs)
    return MoDES(data=data, condition_col="condition", tss_map=tss,
                 external_links=links, fdr_threshold=fdr).run()


def evaluate(result, true_conc_genes, true_primed_genes, true_null_genes):
    """Compute benchmark metrics."""
    et = result.event_table
    n_events = len(et)
    if n_events == 0:
        return {"n_events": 0, "non_null_rate": 0, "false_concordant_rate": 0,
                "true_concordant_rate": 0, "concordant_count": 0, "primed_count": 0}

    pred = {}
    for _, r in et.iterrows():
        g = str(r["gene"])
        if g not in pred:
            pred[g] = r["state"]

    n_conc = sum(1 for g in true_conc_genes if pred.get(g) in ("concordant", "epigenomic_concordant"))
    n_primed = sum(1 for g in true_primed_genes if pred.get(g) in ("chromatin_primed", "active_enhancer_primed"))
    n_null_correct = sum(1 for g in true_null_genes if pred.get(g) in ("null", "unresolved"))

    non_null = sum(1 for v in pred.values() if v not in ("null", "unresolved"))
    non_null_rate = non_null / max(len(pred), 1)
    false_conc = sum(1 for g in true_null_genes if pred.get(g) in ("concordant", "epigenomic_concordant"))
    false_conc_rate = false_conc / max(len(true_null_genes), 1)

    return {
        "n_events": n_events,
        "non_null_rate": round(non_null_rate, 3),
        "false_concordant_rate": round(false_conc_rate, 3),
        "true_concordant_rate": round(n_conc / max(len(true_conc_genes), 1), 3),
        "true_primed_rate": round(n_primed / max(len(true_primed_genes), 1), 3),
    }


def main():
    rng = np.random.default_rng(42)
    all_metrics = []

    # --- 1. Null calibration: shuffle condition labels ---
    print("=" * 60)
    print("MoDES v2.0 Stress Benchmarks")
    print("=" * 60)
    print("[1/4] Null calibration (shuffled condition labels)...")
    rna_df, atac_df, obs, tss, genes = _build_data(rng, with_effect=True)
    true_conc = genes[:10]; true_primed = genes[10:20]; true_null = genes[20:]
    obs_shuffled = obs.copy()
    obs_shuffled["condition"] = rng.permutation(obs["condition"].values)
    data_null = MoDEData(rna=rna_df, atac=atac_df, obs=obs_shuffled)
    result_null = MoDES(data=data_null, condition_col="condition", tss_map=tss).run()
    m = evaluate(result_null, true_conc, true_primed, true_null)
    m["test"] = "null_calibration"
    all_metrics.append(m)
    print(f"    non_null_rate={m['non_null_rate']:.3f} false_concordant_rate={m['false_concordant_rate']:.3f}")

    # --- 2. Link-noise: coordinate-based vs random links vs true links ---
    print("[2/4] Link-noise (coordinate-based vs random vs true links)...")
    rna_df, atac_df, obs, tss, genes = _build_data(rng, with_effect=True)
    true_conc = genes[:10]; true_primed = genes[10:20]; true_null = genes[20:]

    # Coordinate-based (no external links) — baseline
    result_coord = run_moDES(rna_df, atac_df, obs, tss)
    m_coord = evaluate(result_coord, true_conc, true_primed, true_null)
    m_coord["test"] = "link_coordinate_based"
    all_metrics.append(m_coord)
    print(f"    coord: non_null_rate={m_coord['non_null_rate']:.3f} true_conc_rate={m_coord['true_concordant_rate']:.3f}")

    # True external links for signal genes
    true_links = pd.DataFrame({
        "peak_id": [f"chr1:{1000000+i*500000-500}-{1000000+i*500000+500}" for i in range(20)],
        "gene": genes[:20],
    })
    result_true = run_moDES(rna_df, atac_df, obs, tss, links=true_links)
    m_true = evaluate(result_true, true_conc, true_primed, true_null)
    m_true["test"] = "link_external_true"
    all_metrics.append(m_true)
    print(f"    external_true: true_conc_rate={m_true['true_concordant_rate']:.3f} n_events={m_true['n_events']}")

    # Random links added to coordinate-based
    random_links = pd.DataFrame({
        "peak_id": rng.choice(atac_df.columns, 10),
        "gene": rng.choice(rna_df.columns, 10),
    })
    result_rand = run_moDES(rna_df, atac_df, obs, tss, links=random_links)
    m_rand = evaluate(result_rand, true_conc, true_primed, true_null)
    m_rand["test"] = "link_random_extra"
    all_metrics.append(m_rand)
    print(f"    random_extra: false_conc_rate={m_rand['false_concordant_rate']:.3f}")

    # --- 3. Batch-confounded ---
    print("[3/4] Batch-confounded (partial vs full)...")
    for confound_level in ["partial", "full"]:
        rna_df, atac_df, obs, tss, genes = _build_data(rng, with_effect=True)
        true_conc = genes[:10]; true_primed = genes[10:20]; true_null = genes[20:]
        n = len(obs)
        if confound_level == "full":
            obs["batch"] = np.where(obs["condition"] == "ctrl", "A", "B")
        else:
            half = n // 2
            obs["batch"] = ["A"] * half + ["B"] * (n - half)
        data_batch = MoDEData(rna=rna_df, atac=atac_df, obs=obs)
        try:
            result_batch = MoDES(data=data_batch, condition_col="condition",
                                batch_col="batch", tss_map=tss).run()
            m = evaluate(result_batch, true_conc, true_primed, true_null)
            m["test"] = f"batch_{confound_level}"
            all_metrics.append(m)
            print(f"    {confound_level}: non_null_rate={m['non_null_rate']:.3f}")
        except ValueError as e:
            print(f"    {confound_level}: rank deficiency detected (expected)")

    # --- 4. Ablation: no conditional decomp, no quality filter ---
    print("[4/4] Ablation tests...")
    rna_df, atac_df, obs, tss, genes = _build_data(rng, with_effect=True)
    true_conc = genes[:10]; true_primed = genes[10:20]; true_null = genes[20:]

    # Full MoDES v2.0
    result_full = run_moDES(rna_df, atac_df, obs, tss)
    m_full = evaluate(result_full, true_conc, true_primed, true_null)
    m_full["test"] = "ablation_full"
    all_metrics.append(m_full)

    # Without conditional decomposition (set decomposition to NaN for all events)
    data_abl = MoDEData(rna=rna_df, atac=atac_df, obs=obs)
    modes_abl = MoDES(data=data_abl, condition_col="condition", tss_map=tss, fdr_threshold=0.1)
    modes_abl.build_events()
    modes_abl.estimate_effects()
    # Create dummy conditional effects with required columns
    n_ev = len(modes_abl.events)
    modes_abl.conditional_effects = pd.DataFrame({
        "event_id": modes_abl.events["event_id"],
        "rna_after_atac_coef": [np.nan] * n_ev,
        "rna_after_atac_se": [np.nan] * n_ev,
        "rna_after_atac_z": [np.nan] * n_ev,
        "rna_after_atac_pval": [1.0] * n_ev,
        "rna_after_atac_fdr": [1.0] * n_ev,
        "convergence": [False] * n_ev,
        "attenuation": [np.nan] * n_ev,
    })
    modes_abl.build_evidence()
    modes_abl.classify_states()
    modes_abl._build_modality_evidence()
    result_no_decomp = modes_abl._assemble_results()
    m_nd = evaluate(result_no_decomp, true_conc, true_primed, true_null)
    m_nd["test"] = "ablation_no_conditional"
    all_metrics.append(m_nd)

    # Naive DE+DA overlap
    et = result_full.event_table
    naive_conc = len(et[(et["atac_fdr"] < 0.1) & (et["rna_fdr"] < 0.1)])
    naive_total = len(et)
    m_naive = {"test": "ablation_naive_overlap",
               "non_null_rate": round(naive_conc / max(naive_total, 1), 3),
               "false_concordant_rate": 0.0, "true_concordant_rate": 0.0,
               "n_events": naive_total}
    all_metrics.append(m_naive)

    out_dir = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out_dir, exist_ok=True)

    # --- 5. Calibration: assignment_score vs empirical accuracy ---
    print("[5/5] Calibration analysis...")
    et_full = result_full.event_table
    known_genes = dict(zip(genes[:10], ["concordant"] * 10))
    known_genes.update(dict(zip(genes[10:20], ["chromatin_primed"] * 10)))

    cal_rows = []
    for _, row in et_full.iterrows():
        g = row["gene"]
        if g in known_genes:
            cal_rows.append({
                "gene": g,
                "true_state": known_genes[g],
                "pred_state": row["state"],
                "assignment_score": row.get("state_assignment_score", np.nan),
            })
    cal_df = pd.DataFrame(cal_rows)
    if len(cal_df) > 0 and cal_df["assignment_score"].notna().sum() >= 4:
        cal_df["correct"] = cal_df["true_state"] == cal_df["pred_state"]
        cal_df["score_bin"] = pd.cut(
            cal_df["assignment_score"].fillna(0),
            bins=[0, 0.2, 0.5, 0.8, 1.0, 100],
            labels=["0-0.2", "0.2-0.5", "0.5-0.8", "0.8-1.0", ">1.0"],
        )
        cal_summary = cal_df.groupby("score_bin", observed=False).agg(
            n=("correct", "count"),
            n_correct=("correct", "sum"),
            mean_score=("assignment_score", "mean"),
        ).reset_index()
        cal_summary["empirical_accuracy"] = cal_summary["n_correct"] / cal_summary["n"]
        cal_summary["calibration_gap"] = cal_summary["empirical_accuracy"] - cal_summary["mean_score"].clip(0, 1)
        print(cal_summary.to_string(index=False))
        cal_summary.to_csv(
            os.path.join(os.path.dirname(__file__), "output", "calibration.tsv"),
            sep="\t", index=False,
        )

    # --- Output ---
    metrics_df = pd.DataFrame(all_metrics)
    print(f"\nBenchmark Metrics:")
    print(metrics_df.to_string(index=False))
    metrics_df.to_csv(os.path.join(out_dir, "benchmark_metrics.tsv"), sep="\t", index=False)
    print(f"\nOutput: {out_dir}")


if __name__ == "__main__":
    main()
