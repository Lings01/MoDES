"""Tests for MoDES orchestrator and MoDESResult."""

import os

import numpy as np
import pandas as pd
import pytest

from modes import MoDES


class TestMoDES:
    def test_full_pipeline_runs(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small

        modes = MoDES(
            data=data,
            tss_map=tss_map,
            condition_col="condition",
            batch_col="batch",
            fdr_threshold=0.2,
        )

        result = modes.run()
        assert result is not None
        assert len(result.event_table) > 0

    def test_pipeline_outputs_have_expected_columns(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small

        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
        result = modes.run()

        expected_cols = [
            "event_id", "gene", "peak_id", "state", "state_confidence",
            "atac_coef", "rna_coef", "atac_fdr", "rna_fdr",
        ]
        for col in expected_cols:
            assert col in result.event_table.columns, f"Missing column: {col}"

    def test_pipeline_with_covariates(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        data.obs["age"] = np.random.uniform(30, 70, data.n_samples)

        modes = MoDES(
            data=data,
            tss_map=tss_map,
            condition_col="condition",
            covariate_cols=["age"],
            batch_col="batch",
        )
        result = modes.run()
        assert result is not None

    def test_pipeline_with_external_links(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small

        # Create a simple external link
        external = pd.DataFrame({
            "peak_id": [data.peak_names[0]],
            "gene": [data.gene_names[0]],
            "tf_name": ["STAT1"],
            "source": ["scenic"],
        })

        modes = MoDES(
            data=data,
            tss_map=tss_map,
            condition_col="condition",
            external_links=external,
        )
        result = modes.run()
        assert result is not None

    def test_step_by_step(self, synthetic_bulk_data_small):
        """Test that each pipeline step can be called independently."""
        data, gt, tss_map = synthetic_bulk_data_small
        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")

        events = modes.build_events()
        assert events is not None
        assert len(events) > 0

        atac, rna = modes.estimate_effects()
        assert len(atac) > 0
        assert len(rna) > 0

        cond = modes.decompose()
        assert len(cond) == len(events)

        evidence = modes.build_evidence()
        assert len(evidence) == len(events)

        states = modes.classify_states()
        assert len(states) == len(events)

    def test_raises_without_prerequisite(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")

        with pytest.raises(RuntimeError):
            modes.estimate_effects()  # events not built

        modes.build_events()
        with pytest.raises(RuntimeError):
            modes.decompose()  # effects not estimated


class TestMoDESResult:
    @pytest.fixture
    def result(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
        return modes.run()

    def test_summary(self, result):
        text = result.summary()
        assert "MoDES Results Summary" in text
        assert "Total events" in text

    def test_filter_by_state(self, result):
        filtered = result.filter(state="concordant")
        assert all(s == "concordant" for s in filtered["state"])

    def test_filter_by_confidence(self, result):
        filtered = result.filter(min_confidence=0.5)
        assert len(filtered) <= len(result.event_table)

    def test_filter_by_fdr(self, result):
        filtered = result.filter(fdr_threshold=1.0)  # keep all
        assert len(filtered) == len(result.event_table)

        filtered_strict = result.filter(fdr_threshold=0.01)
        assert len(filtered_strict) <= len(result.event_table)

    def test_to_tsv(self, result, tmp_path):
        output_dir = str(tmp_path / "output")
        result.to_tsv(output_dir)

        assert os.path.exists(os.path.join(output_dir, "event_table.tsv"))
        assert os.path.exists(os.path.join(output_dir, "event_state_confidence.tsv"))
        assert os.path.exists(os.path.join(output_dir, "event_layer_effects.tsv"))

        # Verify files are readable
        et = pd.read_csv(os.path.join(output_dir, "event_table.tsv"), sep="\t")
        assert len(et) == len(result.event_table)

    def test_to_graphml(self, result, tmp_path):
        path = str(tmp_path / "network.graphml")
        result.to_graphml(path)
        assert os.path.exists(path)

    def test_to_report(self, result, tmp_path):
        path = str(tmp_path / "report.html")
        result.to_report(path)
        assert os.path.exists(path)

        # Check it's valid HTML
        with open(path) as f:
            content = f.read()
        assert "<html" in content
        assert "MoDES" in content


def test_event_table_contains_artifact_risk(synthetic_bulk_data_small):
    """Event table should contain artifact_risk, artifact_reason, event_pval, event_fdr."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()
    assert "artifact_risk" in result.event_table.columns
    assert "artifact_reason" in result.event_table.columns
    assert "event_pval" in result.event_table.columns
    assert "event_fdr" in result.event_table.columns
    assert result.event_table["event_fdr"].between(0, 1).all()
    # artifact_risk values should be from the allowed set
    valid_risks = {"low", "medium", "high"}
    assert set(result.event_table["artifact_risk"].dropna().unique()).issubset(valid_risks)


def test_filter_exclude_high_artifact(synthetic_bulk_data_small):
    """Filter should support exclude_high_artifact."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()
    filtered = result.filter(exclude_high_artifact=True)
    if "artifact_risk" in filtered.columns:
        assert "high" not in set(filtered["artifact_risk"])


def test_to_tsv_writes_model_diagnostics(synthetic_bulk_data_small, tmp_path):
    """to_tsv should write model_diagnostics.tsv."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()
    result.to_tsv(str(tmp_path))
    assert os.path.exists(os.path.join(str(tmp_path), "model_diagnostics.tsv"))
    diag = pd.read_csv(os.path.join(str(tmp_path), "model_diagnostics.tsv"), sep="\t")
    assert "model_used" in diag.columns
    assert "family" in diag.columns
    assert "converged" in diag.columns
    assert "dropped_covariates" in diag.columns


def test_no_annotation_raises_clear_error(synthetic_bulk_data_small):
    """Plain gene symbols without annotation should raise clear error."""
    data, _, tss_map = synthetic_bulk_data_small

    # Replace gene names with plain symbols (no coordinates)
    data.rna.columns = ["STAT1", "GZMB", "IL7R"] + list(data.rna.columns[3:])

    _modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")  # noqa: F841

    # tss_map only has entries for original gene names, not STAT1/GZMB/IL7R
    # Since some genes still match and some don't, we need to check that
    # external_links or a proper tss_map covering ALL genes is needed.
    # For this test, remove tss_map to trigger the error
    modes2 = MoDES(data=data, condition_col="condition")
    with pytest.raises(ValueError, match="No candidate events"):
        modes2.build_events()


def test_external_links_dataframe_not_boolean_evaluated(synthetic_bulk_data_small):
    """DataFrame external_links should not trigger boolean ambiguity."""
    data, _, tss_map = synthetic_bulk_data_small
    links = pd.DataFrame({
        "peak_id": [data.peak_names[0]],
        "gene": [data.gene_names[0]],
    })
    modes = MoDES(
        data=data,
        tss_map=tss_map,
        condition_col="condition",
        external_links=links,
    )
    events = modes.build_events()
    assert len(events) > 0


def test_report_escapes_html(synthetic_bulk_data_small, tmp_path):
    """HTML report should escape script tags in gene names."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    result.event_table.loc[0, "gene"] = "<script>alert(1)</script>"
    path = str(tmp_path / "report.html")
    result.to_report(path)

    with open(path) as f:
        html_text = f.read()
    assert "<script>" not in html_text
    assert "&lt;script&gt;" in html_text


def test_event_table_schema_exact(synthetic_bulk_data_small):
    """Output schema must match the frozen v1.1 column order."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    expected = [
        "event_id", "tf_name", "peak_id", "gene", "context",
        "atac_coef", "atac_se", "atac_pval", "atac_fdr", "atac_direction",
        "rna_coef", "rna_se", "rna_pval", "rna_fdr", "rna_direction",
        "rna_after_atac_coef", "rna_after_atac_se",
        "rna_after_atac_pval", "rna_after_atac_fdr",
        "state", "state_confidence", "quality_score",
        "artifact_risk", "artifact_reason",
        "event_pval", "event_fdr",
    ]
    assert list(result.event_table.columns) == expected


def test_save_load_roundtrip(synthetic_bulk_data_small, tmp_path):
    """save() + load() should produce identical event_tables."""
    from modes import MoDESResult
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    out = str(tmp_path / "results")
    result.save(out)
    loaded = MoDESResult.load(out)

    assert list(result.event_table.columns) == list(loaded.event_table.columns)
    assert len(result.event_table) == len(loaded.event_table)
    assert result.params == loaded.params


def test_filter_states(synthetic_bulk_data_small):
    """Filter by states list should work."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    filtered = result.filter(states=["concordant", "chromatin_primed"])
    assert set(filtered["state"].unique()).issubset({"concordant", "chromatin_primed"})


def test_filter_genes(synthetic_bulk_data_small):
    """Filter by gene list should work."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    genes = result.event_table["gene"].iloc[:2].tolist()
    filtered = result.filter(genes=genes)
    assert set(filtered["gene"].unique()).issubset(set(genes))


def test_filter_peaks(synthetic_bulk_data_small):
    """Filter by peak list should work."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    peaks = result.event_table["peak_id"].iloc[:2].tolist()
    filtered = result.filter(peaks=peaks)
    assert set(filtered["peak_id"].unique()).issubset(set(peaks))


def test_run_params_has_versions(synthetic_bulk_data_small):
    """run_params should include version info."""
    data, _, tss_map = synthetic_bulk_data_small
    modes = MoDES(data=data, tss_map=tss_map, condition_col="condition")
    result = modes.run()

    assert "modes_version" in result.params
    assert "python_version" in result.params
    assert "numpy_version" in result.params
    assert "n_external_links" in result.params
