"""Tests for MoDES orchestrator and MoDESResult."""

import pytest
import os
import numpy as np
import pandas as pd

from modes import MoDES, MoDESResult, MoDEData


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
            "event_id", "gene", "peak_id", "state", "posterior",
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

    def test_filter_by_posterior(self, result):
        filtered = result.filter(min_posterior=0.5)
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
        assert os.path.exists(os.path.join(output_dir, "event_state_probability.tsv"))
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
