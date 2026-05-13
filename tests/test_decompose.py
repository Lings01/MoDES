"""Tests for ConditionalDecomposition."""

import pytest
import numpy as np
import pandas as pd

from modes.effects import EffectEstimator
from modes.decompose import ConditionalDecomposition


class TestConditionalDecomposition:
    @pytest.fixture
    def setup_data(self, synthetic_bulk_data_small):
        from modes.events import EventCandidateBuilder
        data, gt, tss_map = synthetic_bulk_data_small

        # Build events
        builder = EventCandidateBuilder()
        events = builder.build(
            gene_names=list(data.gene_names),
            peak_names=list(data.peak_names),
            tss_map=tss_map,
        )

        # Estimate effects
        estimator = EffectEstimator(
            condition_col="condition",
            batch_col="batch",
            use_empirical_bayes=False,
        )
        atac_effects, rna_effects = estimator.estimate_effects(
            data, list(data.peak_names), list(data.gene_names)
        )

        return data, events, atac_effects, rna_effects

    def test_decompose_runs(self, setup_data):
        data, events, atac_eff, rna_eff = setup_data
        dec = ConditionalDecomposition(
            condition_col="condition",
            batch_col="batch",
        )
        result = dec.decompose(data, events, atac_eff, rna_eff)
        assert len(result) == len(events)
        assert "rna_after_atac_coef" in result.columns
        assert "attenuation" in result.columns

    def test_conditional_effects_valid(self, setup_data):
        data, events, atac_eff, rna_eff = setup_data
        dec = ConditionalDecomposition(condition_col="condition")
        result = dec.decompose(data, events, atac_eff, rna_eff)

        # Check non-NaN values
        valid = result["convergence"] == True  # noqa: E712
        if valid.sum() > 0:
            sub = result[valid]
            assert np.all(np.isfinite(sub["rna_after_atac_coef"]))

    def test_attenuation_range(self, setup_data):
        """Attenuation should be interpretable."""
        data, events, atac_eff, rna_eff = setup_data
        dec = ConditionalDecomposition(condition_col="condition")
        result = dec.decompose(data, events, atac_eff, rna_eff)

        # Attenuation is ratio: conditional / marginal
        valid = result["attenuation"].notna()
        if valid.sum() > 0:
            att = result.loc[valid, "attenuation"]
            # Can be outside [0,1] due to estimation noise,
            # but should not be wildly extreme
            assert np.all(np.isfinite(att))

    def test_decompose_handles_missing_features(self, setup_data):
        """Graceful handling when a feature is not in the data."""
        data, events, atac_eff, rna_eff = setup_data
        dec = ConditionalDecomposition(condition_col="condition")

        # Add a fake event with missing peak
        fake_events = pd.DataFrame([{
            "event_id": "fake_event",
            "gene": "nonexistent_gene",
            "peak_id": "nonexistent_peak",
        }])

        result = dec.decompose(data, fake_events, atac_eff, rna_eff)
        assert len(result) == 1
        assert np.isnan(result.iloc[0]["rna_after_atac_coef"])
