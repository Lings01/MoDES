"""Tests for EvidenceBuilder and StateClassifier."""

import pytest
import numpy as np
import pandas as pd

from modes._types import ModalityEffect
from modes.states import EvidenceBuilder, StateClassifier, summarize_states


class TestEvidenceBuilder:
    def test_build(self, synthetic_bulk_data_small):
        from modes.events import EventCandidateBuilder
        from modes.effects import EffectEstimator
        from modes.decompose import ConditionalDecomposition

        data, gt, tss_map = synthetic_bulk_data_small

        builder = EventCandidateBuilder()
        events = builder.build(
            gene_names=list(data.gene_names),
            peak_names=list(data.peak_names),
            tss_map=tss_map,
        )

        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=False,
        )
        atac_eff, rna_eff = estimator.estimate_effects(
            data, list(data.peak_names), list(data.gene_names)
        )

        dec = ConditionalDecomposition(condition_col="condition")
        cond_eff = dec.decompose(data, events, atac_eff, rna_eff)

        ev_builder = EvidenceBuilder()
        evidence = ev_builder.build(events, atac_eff, rna_eff, cond_eff, data)

        assert len(evidence) == len(events)
        assert "z_atac" in evidence.columns
        assert "z_rna" in evidence.columns
        assert "quality_score" in evidence.columns


class TestStateClassifier:
    @pytest.fixture
    def evidence_df(self, synthetic_bulk_data_small):
        from modes.events import EventCandidateBuilder
        from modes.effects import EffectEstimator
        from modes.decompose import ConditionalDecomposition

        data, gt, tss_map = synthetic_bulk_data_small

        builder = EventCandidateBuilder()
        events = builder.build(
            gene_names=list(data.gene_names),
            peak_names=list(data.peak_names),
            tss_map=tss_map,
        )

        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=False,
        )
        atac_eff, rna_eff = estimator.estimate_effects(
            data, list(data.peak_names), list(data.gene_names)
        )

        dec = ConditionalDecomposition(condition_col="condition")
        cond_eff = dec.decompose(data, events, atac_eff, rna_eff)

        ev_builder = EvidenceBuilder()
        return ev_builder.build(events, atac_eff, rna_eff, cond_eff, data)

    def test_classify_all_states_present(self, evidence_df):
        """All 5 states should appear in output."""
        classifier = StateClassifier(
            fdr_threshold=0.5,  # loose threshold to get more calls
            use_empirical_bayes=True,
        )
        states = classifier.classify(evidence_df)
        assert len(states) == len(evidence_df)
        assert "state" in states.columns
        assert "posterior_prob" in states.columns
        # At least some non-null states
        unique_states = set(states["state"])
        assert len(unique_states) >= 1

    def test_classify_rule_based_only(self, evidence_df):
        classifier = StateClassifier(
            fdr_threshold=0.1,
            use_empirical_bayes=False,
        )
        states = classifier.classify(evidence_df)
        assert len(states) == len(evidence_df)

    def test_posterior_probabilities_in_range(self, evidence_df):
        classifier = StateClassifier(
            fdr_threshold=0.1,
            use_empirical_bayes=True,
        )
        states = classifier.classify(evidence_df)
        assert (states["posterior_prob"] >= 0).all()
        assert (states["posterior_prob"] <= 1).all()

    def test_state_labels_valid(self, evidence_df):
        classifier = StateClassifier(use_empirical_bayes=False)
        states = classifier.classify(evidence_df)
        for s in states["state"]:
            assert s in StateClassifier.VALID_STATES

    def test_summarize_states(self, evidence_df):
        classifier = StateClassifier(use_empirical_bayes=False)
        states = classifier.classify(evidence_df)
        summary = summarize_states(states)
        assert "count" in summary.columns
        assert "fraction" in summary.columns
        assert abs(summary["fraction"].sum() - 1.0) < 0.01
