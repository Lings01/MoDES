"""Tests for EvidenceBuilder and StateClassifier (v2.0 grammar-based)."""

import pandas as pd
import numpy as np
import pytest

from modes.states import EvidenceBuilder, StateClassifier, summarize_states


class TestEvidenceBuilder:
    def test_build(self, synthetic_bulk_data_small):
        from modes.decompose import ConditionalDecomposition
        from modes.effects import EffectEstimator
        from modes.events import EventCandidateBuilder

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
        from modes.decompose import ConditionalDecomposition
        from modes.effects import EffectEstimator
        from modes.events import EventCandidateBuilder

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
        """Output should have new v2.0 columns."""
        classifier = StateClassifier(fdr_threshold=0.5)
        states = classifier.classify(evidence_df)
        assert len(states) == len(evidence_df)
        assert "state" in states.columns
        assert "state_assignment_score" in states.columns
        assert "state_support_score" in states.columns
        unique_states = set(states["state"])
        assert len(unique_states) >= 1

    def test_classify_rule_based_only(self, evidence_df):
        classifier = StateClassifier(fdr_threshold=0.1)
        states = classifier.classify(evidence_df)
        assert len(states) == len(evidence_df)

    def test_state_assignment_score_in_range(self, evidence_df):
        classifier = StateClassifier(fdr_threshold=0.1)
        states = classifier.classify(evidence_df)
        scores = states["state_assignment_score"].dropna()
        assert (scores >= 0).all()

    def test_state_labels_valid(self, evidence_df):
        classifier = StateClassifier(fdr_threshold=0.1)
        states = classifier.classify(evidence_df)
        from modes.modalities.state_rules import ALL_STATE_NAMES
        for s in states["state"]:
            assert s in ALL_STATE_NAMES

    def test_summarize_states(self, evidence_df):
        classifier = StateClassifier(fdr_threshold=0.1)
        states = classifier.classify(evidence_df)
        summary = summarize_states(states)
        assert "count" in summary.columns
        assert "fraction" in summary.columns
        assert abs(summary["fraction"].sum() - 1.0) < 0.01


def test_low_quality_significant_event_gets_artifact_risk():
    """Low-quality single-modality signal should get high artifact_risk."""
    evidence = pd.DataFrame({
        "event_id": ["e1"],
        "z_atac": [5.0],
        "z_rna": [0.1],
        "z_rna_given_atac": [0.1],
        "atac_fdr": [1e-6],
        "rna_fdr": [1.0],
        "atac_direction": [1],
        "rna_direction": [0],
        "atac_pval": [1e-6],
        "rna_pval": [1.0],
        "atac_coef": [2.0],
        "rna_coef": [0.0],
        "quality_score": [0.05],
    })

    classifier = StateClassifier(
        fdr_threshold=0.1,
        quality_threshold=0.3,
    )
    states = classifier.classify(evidence)
    # Low quality single-modality signal: ATAC-only is chromatin_primed
    assert states.loc[0, "state"] == "chromatin_open_primed"
    assert states.loc[0, "artifact_risk"] == "high"
    assert "single_modality_low_quality" in states.loc[0, "artifact_reason"]


def test_artifact_risk_low_medium_high():
    """artifact_risk should be low/medium/high based on quality_score."""
    evidence = pd.DataFrame({
        "event_id": [f"e{i}" for i in range(3)],
        "z_atac": [0.1, 0.1, 0.1],
        "z_rna": [0.1, 0.1, 0.1],
        "z_rna_given_atac": [0.1, 0.1, 0.1],
        "atac_fdr": [1.0, 1.0, 1.0],
        "rna_fdr": [1.0, 1.0, 1.0],
        "atac_direction": [0, 0, 0],
        "rna_direction": [0, 0, 0],
        "atac_pval": [1.0, 1.0, 1.0],
        "rna_pval": [1.0, 1.0, 1.0],
        "atac_coef": [0.0, 0.0, 0.0],
        "rna_coef": [0.0, 0.0, 0.0],
        "quality_score": [0.9, 0.5, 0.1],
    })

    classifier = StateClassifier(quality_threshold=0.3)
    states = classifier.classify(evidence)
    assert states.loc[0, "artifact_risk"] == "low"
    assert states.loc[1, "artifact_risk"] == "medium"
    # quality=0.1 < 0.3: borderline/low_quality_score → medium
    assert states.loc[2, "artifact_risk"] == "medium"


def test_rule_based_core_states_exact():
    """Direct evidence-to-state test: grammar rules classify correctly."""
    evidence = pd.DataFrame({
        "event_id": ["e_conc", "e_primed", "e_rna", "e_null"],
        "z_atac": [5.0, 5.0, 0.1, 0.1],
        "z_rna": [4.0, 0.1, 4.5, 0.1],
        "z_rna_given_atac": [0.2, 0.1, 4.0, 0.1],
        "atac_fdr": [1e-6, 1e-6, 1.0, 1.0],
        "rna_fdr": [1e-6, 1.0, 1e-6, 1.0],
        "atac_direction": [1, 1, 0, 0],
        "rna_direction": [1, 0, 1, 0],
        "atac_pval": [1e-6, 1e-6, 1.0, 1.0],
        "rna_pval": [1e-6, 1.0, 1e-6, 1.0],
        "atac_coef": [2.0, 2.0, 0.0, 0.0],
        "rna_coef": [1.5, 0.0, 1.5, 0.0],
        "quality_score": [1.0, 1.0, 1.0, 1.0],
    })
    clf = StateClassifier(fdr_threshold=0.1)
    states = clf.classify(evidence).set_index("event_id")
    assert states.loc["e_conc", "state"] == "concordant_activation"
    assert states.loc["e_primed", "state"] == "chromatin_open_primed"
    assert states.loc["e_rna", "state"] == "rna_up_only"
    assert states.loc["e_null", "state"] == "null"


# ── P1: Full-direction RNA+ATAC state tests ───────────────────────────

def _make_evidence(event_id, atac_fdr, atac_dir, atac_coef, rna_fdr, rna_dir, rna_coef, quality=1.0):
    return pd.DataFrame({
        "event_id": [event_id],
        "z_atac": [8.0 if atac_fdr < 0.01 else 0.1],
        "z_rna": [8.0 if rna_fdr < 0.01 else 0.1],
        "z_rna_given_atac": [0.1],
        "atac_fdr": [atac_fdr], "rna_fdr": [rna_fdr],
        "atac_direction": [atac_dir], "rna_direction": [rna_dir],
        "atac_pval": [atac_fdr], "rna_pval": [rna_fdr],
        "atac_coef": [atac_coef], "rna_coef": [rna_coef],
        "quality_score": [quality],
    })


def test_concordant_activation():
    ev = _make_evidence("e1", 1e-6, 1, 1.5, 1e-6, 1, 2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "concordant_activation"


def test_concordant_repression():
    ev = _make_evidence("e1", 1e-6, -1, -1.5, 1e-6, -1, -2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "concordant_repression"


def test_discordant_opening_repression():
    ev = _make_evidence("e1", 1e-6, 1, 1.5, 1e-6, -1, -2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "discordant_opening_repression"


def test_discordant_closing_activation():
    ev = _make_evidence("e1", 1e-6, -1, -1.5, 1e-6, 1, 2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "discordant_closing_activation"


def test_chromatin_open_primed():
    ev = _make_evidence("e1", 1e-6, 1, 1.5, 1.0, 0, 0.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "chromatin_open_primed"


def test_chromatin_closed_primed():
    ev = _make_evidence("e1", 1e-6, -1, -1.5, 1.0, 0, 0.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "chromatin_closed_primed"


def test_rna_up_only():
    ev = _make_evidence("e1", 1.0, 0, 0.0, 1e-6, 1, 2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "rna_up_only"


def test_rna_down_only():
    ev = _make_evidence("e1", 1.0, 0, 0.0, 1e-6, -1, -2.0)
    clf = StateClassifier(fdr_threshold=0.1)
    out = clf.classify(ev)
    assert out.loc[0, "state"] == "rna_down_only"


# ── P1: Filter tests ─────────────────────────────────────────────────

def test_filter_uses_min_state_support_adjusted_score():
    from modes import MoDESResult
    et = pd.DataFrame({
        "event_id": ["e1", "e2", "e3"],
        "state": ["null", "null", "null"],
        "state_support_adjusted_score": [0.5, 2.0, 5.0],
    })
    result = MoDESResult(event_table=et)
    out = result.filter(min_state_support_adjusted_score=2.0)
    assert set(out["state_support_adjusted_score"].values) == {2.0, 5.0}


def test_filter_state_family():
    from modes import MoDESResult
    et = pd.DataFrame({
        "event_id": ["e1", "e2", "e3"],
        "state": ["concordant_activation", "concordant_repression", "rna_up_only"],
        "state_family": ["concordant", "concordant", "rna_only"],
    })
    result = MoDESResult(event_table=et)
    out = result.filter(state_family="concordant")
    assert len(out) == 2
    assert set(out["state"]) == {"concordant_activation", "concordant_repression"}


# ── P1: CUT&Tag mark_only test ────────────────────────────────────────

def test_mark_only_wins_when_mark_up_rna_atac_absent():
    """mark_only should win over active_enhancer_primed when both RNA and ATAC absent."""
    evidence = pd.DataFrame({
        "event_id": ["e1"],
        "atac_fdr": [1.0], "atac_direction": [0], "atac_pval": [1.0], "atac_coef": [0.0],
        "rna_fdr": [1.0], "rna_direction": [0], "rna_pval": [1.0], "rna_coef": [0.0],
        "z_atac": [0.1], "z_rna": [0.1], "z_rna_given_atac": [0.1],
        "h3k27ac_z": [8.0], "h3k27ac_fdr": [1e-8], "h3k27ac_pval": [1e-8],
        "h3k27ac_direction": [1], "h3k27ac_coef": [2.0],
        "h3k27ac_available": [True], "h3k27ac_matched": [True],
        "h3k27ac_measured": [True], "h3k27ac_missing_reason": [""],
        "quality_score": [1.0],
    })
    from modes.modalities.base import ModalitySpec
    clf = StateClassifier(
        fdr_threshold=0.1,
        modality_specs={
            "h3k27ac": ModalitySpec(
                name="h3k27ac", assay="CUTTAG", feature_type="region",
                target="H3K27ac", regulatory_role="activating_mark",
                expected_rna_direction=1,
            )
        },
    )
    out = clf.classify(evidence)
    assert out.loc[0, "state"] == "mark_only", f"Expected mark_only, got {out.loc[0, 'state']}"


# ── P1: Protein absent/missing tests ──────────────────────────────────

def test_protein_buffered_requires_measured_protein_absent():
    """protein_buffered_up should trigger when protein matched+measured+not sig."""
    evidence = pd.DataFrame({
        "event_id": ["e1"],
        "atac_fdr": [1.0], "atac_direction": [0], "atac_pval": [1.0], "atac_coef": [0.0],
        "rna_fdr": [1e-6], "rna_direction": [1], "rna_pval": [1e-6], "rna_coef": [2.0],
        "z_atac": [0.1], "z_rna": [8.0], "z_rna_given_atac": [0.1],
        "protein_z": [0.1], "protein_fdr": [0.5], "protein_pval": [0.5],
        "protein_direction": [0], "protein_coef": [0.0],
        "protein_available": [True], "protein_matched": [True],
        "protein_measured": [True], "protein_missing_reason": [""],
        "quality_score": [1.0],
    })
    from modes.modalities.base import ModalitySpec
    clf = StateClassifier(
        fdr_threshold=0.1,
        modality_specs={
            "protein": ModalitySpec(name="protein", assay="PROTEIN",
                                     feature_type="protein", regulatory_role="protein_output"),
        },
    )
    out = clf.classify(evidence)
    assert out.loc[0, "state"] == "protein_buffered_up", f"Got {out.loc[0, 'state']}"


def test_protein_missing_does_not_trigger_buffered():
    """When protein modality is unmatched, protein_buffered should NOT trigger."""
    evidence = pd.DataFrame({
        "event_id": ["e1"],
        "atac_fdr": [1.0], "atac_direction": [0], "atac_pval": [1.0], "atac_coef": [0.0],
        "rna_fdr": [1e-6], "rna_direction": [1], "rna_pval": [1e-6], "rna_coef": [2.0],
        "z_atac": [0.1], "z_rna": [8.0], "z_rna_given_atac": [0.1],
        "protein_z": [0.0], "protein_fdr": [float("nan")], "protein_pval": [float("nan")],
        "protein_direction": [0], "protein_coef": [float("nan")],
        "protein_available": [True], "protein_matched": [False],
        "protein_measured": [False], "protein_missing_reason": ["no_feature_match"],
        "quality_score": [1.0],
    })
    from modes.modalities.base import ModalitySpec
    clf = StateClassifier(
        fdr_threshold=0.1,
        modality_specs={
            "protein": ModalitySpec(name="protein", assay="PROTEIN",
                                     feature_type="protein", regulatory_role="protein_output"),
        },
    )
    out = clf.classify(evidence)
    assert out.loc[0, "state"] != "protein_buffered_up", \
        f"Should NOT be protein_buffered_up when protein unmatched, got {out.loc[0, 'state']}"
