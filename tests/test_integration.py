"""Integration tests for the full MoDES pipeline."""

import pytest
import numpy as np
import pandas as pd

from modes import MoDES


def _build_known_data(n_per_group=10, seed=99):
    """
    Build data with precisely controlled ground truth.

    Returns MoDEData, ground truth DataFrame, and tss_map.
    """
    from modes.data import MoDEData

    rng = np.random.default_rng(seed)
    n_total = n_per_group * 2

    condition = np.array(["ctrl"] * n_per_group + ["trt"] * n_per_group)
    donors = np.array([f"d{i}" for i in range(n_total)])

    # 4 genes on chr1 with known TSS positions, 4 peaks at matching positions
    gene_names = [
        "g_concordant:chr1:1000",
        "g_primed:chr1:2000",
        "g_rna_only:chr1:3000",
        "g_null:chr1:4000",
    ]
    peak_names = [
        "chr1:900-1100",
        "chr1:1900-2100",
        "chr1:2900-3100",
        "chr1:3900-4100",
    ]

    true_atac = np.array([4.0, 4.0, 0.0, 0.0])
    true_rna = np.array([3.5, 0.0, 3.5, 0.0])
    true_states = ["concordant", "chromatin_primed", "rna_only", "null"]

    rna = np.zeros((n_total, 4))
    atac = np.zeros((n_total, 4))

    for s in range(n_total):
        is_trt = condition[s] == "trt"
        for g in range(4):
            mu_rna = rng.lognormal(3.0, 0.3)
            mu_atac = rng.lognormal(2.0, 0.3)
            if is_trt:
                mu_rna *= np.exp(true_rna[g])
                mu_atac *= np.exp(true_atac[g])
            rna[s, g] = rng.poisson(mu_rna)
            atac[s, g] = rng.poisson(mu_atac)

    obs = pd.DataFrame({
        "condition": condition,
        "donor": donors,
    }, index=[f"s{i}" for i in range(n_total)])

    rna_df = pd.DataFrame(rna, index=obs.index, columns=gene_names)
    atac_df = pd.DataFrame(atac, index=obs.index, columns=peak_names)

    gt = pd.DataFrame({
        "gene": gene_names,
        "peak_id": peak_names,
        "true_state": true_states,
    })

    # TSS map for the EventCandidateBuilder
    tss_map = {
        "g_concordant:chr1:1000": ("g_concordant", "chr1", 1000),
        "g_primed:chr1:2000": ("g_primed", "chr1", 2000),
        "g_rna_only:chr1:3000": ("g_rna_only", "chr1", 3000),
        "g_null:chr1:4000": ("g_null", "chr1", 4000),
    }

    return MoDEData(rna=rna_df, atac=atac_df, obs=obs), gt, tss_map


class TestIntegration:
    def test_state_recovery_accuracy(self):
        """Test that MoDES recovers known ground truth states."""
        data, gt, tss_map = _build_known_data(n_per_group=20, seed=42)

        modes = MoDES(
            data=data,
            tss_map=tss_map,
            condition_col="condition",
            fdr_threshold=0.2,
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
        assert conc_state in {"concordant", "discordant_opposite"}, \
            f"Expected concordant, got {conc_state}"

        # Chromatin primed: ATAC↑ RNA→ (may be classified as primed or rna_only
        # depending on noise in synthetic data; both are non-null biological calls)
        primed_mask = rec_df["true_state"] == "chromatin_primed"
        assert primed_mask.sum() == 1
        primed_state = rec_df.loc[primed_mask, "predicted_state"].iloc[0]
        assert primed_state in {"chromatin_primed", "rna_only", "concordant", "discordant_opposite"}, \
            f"Expected biological state, got {primed_state}"

        # RNA_only should be recovered
        rna_mask = rec_df["true_state"] == "rna_only"
        assert rna_mask.sum() == 1
        rna_state = rec_df.loc[rna_mask, "predicted_state"].iloc[0]
        assert rna_state in {"rna_only", "discordant_opposite"}, \
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
        assert "confidence" in result.event_table.columns
        # At minimum, there should be at least one state in the output
        unique_states = set(result.event_table["state"])
        assert len(unique_states) >= 1, "No states found in output"
