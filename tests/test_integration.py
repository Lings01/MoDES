"""Integration tests for the full MoDES pipeline."""

import numpy as np
import pandas as pd

from modes import MoDES


def _build_known_data():
    """
    Build deterministic data with precisely controlled ground truth.

    Design: null features use fixed baseline + symmetric Normal(0,sigma) noise
    so that E[ctrl] = E[trt] by construction. Signal features use large
    multiplicative fold changes for reliable detection even at small n.

    Returns MoDEData, ground truth DataFrame, and tss_map.
    """
    from modes.data import MoDEData

    # Only 4 genes/peaks with coordinates → only 4 real events
    gene_names = [
        "g_concordant:chr1:500000",
        "g_primed:chr2:500000",
        "g_rna_only:chr3:500000",
        "g_null:chr4:500000",
    ]
    peak_names = [
        "chr1:499000-501000",
        "chr2:499000-501000",
        "chr3:499000-501000",
        "chr4:499000-501000",
    ]
    # Filler genes/peaks to balance library sizes (no coordinates → no events)
    n_filler = 100
    filler_genes = [f"filler_gene_{i}" for i in range(n_filler)]
    filler_peaks = [f"filler_peak_{i}" for i in range(n_filler)]
    all_gene_names = gene_names + filler_genes
    all_peak_names = peak_names + filler_peaks

    # 10 ctrl + 10 trt = 20 samples. Enough power for 30x signal, low power for null.
    n = 10
    n_total = n * 2
    condition = np.array(["ctrl"] * n + ["trt"] * n)

    rng = np.random.default_rng(42)

    # For null modalities: trt uses the SAME n values as ctrl (reversed order).
    # This guarantees identical group means → coef=0, SE>0, p≈1.0.
    atac_null_ctrl = rng.poisson(200, n).astype(float)
    rna_null_ctrl = rng.poisson(300, n).astype(float)
    atac_null_trt = atac_null_ctrl[::-1].copy()  # same multiset, reversed
    rna_null_trt = rna_null_ctrl[::-1].copy()

    atac = np.zeros((n_total, 4), dtype=float)
    rna = np.zeros((n_total, 4), dtype=float)

    # rna_only ATAC null: identical multisets → p=1.0
    atac[:n, 2] = atac_null_ctrl
    atac[n:, 2] = atac_null_trt
    # primed RNA null: identical multisets → p=1.0
    rna[:n, 1] = rna_null_ctrl
    rna[n:, 1] = rna_null_trt
    # null state: both modalities null (independent Poisson, will be non-sig with small n)
    atac[:, 3] = rng.poisson(200, n_total).astype(float)
    rna[:, 3] = rng.poisson(300, n_total).astype(float)

    # --- Signal modalities: 30x fold change in trt ---
    # Add filler columns to balance library sizes between ctrl and trt.
    # Without balancing, trt's signal genes inflate library size offsets,
    # causing spurious significance in null genes.
    # Use n_filler defined above to balance library sizes
    atac_full = np.zeros((n_total, 4 + n_filler), dtype=float)
    rna_full = np.zeros((n_total, 4 + n_filler), dtype=float)

    # Filler: identical null values for all samples (dominates library size)
    for fi in range(n_filler):
        rna_full[:, 4 + fi] = rng.poisson(500, n_total)
        atac_full[:, 4 + fi] = rng.poisson(400, n_total)

    # Concordant (col 0): ATAC↑ RNA↑
    atac_full[:n, 0] = rng.poisson(100, n)
    rna_full[:n, 0] = rng.poisson(150, n)
    atac_full[n:, 0] = rng.poisson(3000, n)
    rna_full[n:, 0] = rng.poisson(4500, n)

    # Primed (col 1): ATAC↑ (30x), RNA null
    atac_full[:n, 1] = rng.poisson(100, n)
    atac_full[n:, 1] = rng.poisson(3000, n)

    # RNA_only (col 2): RNA↑ (30x), ATAC null
    rna_full[:n, 2] = rng.poisson(150, n)
    rna_full[n:, 2] = rng.poisson(4500, n)

    # Copy null modalities into full arrays
    atac_full[:n, 2] = atac_null_ctrl   # rna_only ATAC null
    atac_full[n:, 2] = atac_null_trt
    rna_full[:n, 1] = rna_null_ctrl     # primed RNA null
    rna_full[n:, 1] = rna_null_trt
    atac_full[:, 3] = rng.poisson(200, n_total).astype(float)  # null ATAC
    rna_full[:, 3] = rng.poisson(300, n_total).astype(float)   # null RNA

    atac = atac_full
    rna = rna_full

    true_states = ["concordant", "chromatin_primed", "rna_only", "null"]

    obs = pd.DataFrame({
        "condition": condition,
    }, index=[f"s{i}" for i in range(n_total)])

    rna_df = pd.DataFrame(rna, index=obs.index, columns=all_gene_names)
    atac_df = pd.DataFrame(atac, index=obs.index, columns=all_peak_names)

    gt = pd.DataFrame({
        "gene": gene_names,
        "peak_id": peak_names,
        "true_state": true_states,
    })

    tss_map = {
        "g_concordant:chr1:500000": ("g_concordant", "chr1", 500000),
        "g_primed:chr2:500000": ("g_primed", "chr2", 500000),
        "g_rna_only:chr3:500000": ("g_rna_only", "chr3", 500000),
        "g_null:chr4:500000": ("g_null", "chr4", 500000),
    }

    return MoDEData(rna=rna_df, atac=atac_df, obs=obs), gt, tss_map


class TestIntegration:
    def test_state_recovery_accuracy(self):
        """Test that MoDES recovers known ground truth states."""
        data, gt, tss_map = _build_known_data()

        modes = MoDES(
            data=data,
            tss_map=tss_map,
            condition_col="condition",
            fdr_threshold=0.1,
        )
        result = modes.run()

        # Map ground truth to events
        recovered = []
        for _, gr in gt.iterrows():
            mask = (
                (result.event_table["gene"] == gr["gene"]) &
                (result.event_table["peak_id"] == gr["peak_id"])
            )
            if mask.sum() > 0:
                row = result.event_table[mask].iloc[0]
                recovered.append({
                    "gene": gr["gene"],
                    "peak_id": gr["peak_id"],
                    "true_state": gr["true_state"],
                    "predicted_state": row["state"],
                })

        rec_df = pd.DataFrame(recovered)
        assert len(rec_df) == 4, "Should have all 4 ground truth events match"

        # Concordant should be recovered as concordant
        conc_mask = rec_df["true_state"] == "concordant"
        assert conc_mask.sum() == 1
        conc_state = rec_df.loc[conc_mask, "predicted_state"].iloc[0]
        assert conc_state == "concordant", f"Expected concordant, got {conc_state}"

        # Chromatin primed should be recovered
        primed_mask = rec_df["true_state"] == "chromatin_primed"
        assert primed_mask.sum() == 1
        primed_state = rec_df.loc[primed_mask, "predicted_state"].iloc[0]
        assert primed_state == "chromatin_primed", \
            f"Expected chromatin_primed, got {primed_state}"

        rna_mask = rec_df["true_state"] == "rna_only"
        assert rna_mask.sum() == 1
        rna_state = rec_df.loc[rna_mask, "predicted_state"].iloc[0]
        assert rna_state == "rna_only", \
            f"Expected rna_only, got {rna_state}"

    def test_full_pipeline_outputs(self, synthetic_bulk_data_small):
        """All output files are generated correctly."""
        import tempfile
        import os

        data, gt, tss_map = synthetic_bulk_data_small
        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
        result = modes.run()

        with tempfile.TemporaryDirectory() as d:
            result.to_tsv(d)
            result.to_graphml(os.path.join(d, "net.graphml"))
            result.to_report(os.path.join(d, "report.html"))

            files = os.listdir(d)
            assert "event_table.tsv" in files
            assert "event_state_confidence.tsv" in files
            assert "event_layer_effects.tsv" in files
            assert "model_diagnostics.tsv" in files
            assert "run_params.tsv" in files
            assert "net.graphml" in files
            assert "report.html" in files

    def test_state_distribution_nonempty(self, synthetic_bulk_data_small):
        """Pipeline should run without error and produce output rows."""
        data, gt, tss_map = synthetic_bulk_data_small
        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition", fdr_threshold=0.2)
        result = modes.run()

        assert len(result.event_table) > 0, "No events in output"
        assert "state" in result.event_table.columns
        assert "state_confidence" in result.event_table.columns
        assert "artifact_risk" in result.event_table.columns
        assert "event_fdr" in result.event_table.columns
        unique_states = set(result.event_table["state"])
        assert len(unique_states) >= 1, "No states found in output"
        # All states must be biological states
        from modes.states import StateClassifier
        assert unique_states.issubset(StateClassifier.BIOLOGICAL_STATES)
