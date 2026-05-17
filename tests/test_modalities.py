"""Tests for the v2.0 modality abstraction layer."""

import pytest
import numpy as np
import pandas as pd

from modes.modalities import (
    ModalitySpec, CUTTAG_REGISTRY,
    get_cuttag_target_info, validate_cuttag_features,
    make_cuttag_spec,
)
from modes.modalities.base import RNA_SPEC, ATAC_SPEC


class TestModalitySpec:
    def test_rna_spec(self):
        assert RNA_SPEC.name == "rna"
        assert RNA_SPEC.assay == "RNA"
        assert RNA_SPEC.feature_type == "gene"

    def test_atac_spec(self):
        assert ATAC_SPEC.name == "atac"
        assert ATAC_SPEC.assay == "ATAC"
        assert ATAC_SPEC.expected_rna_direction == 1

    def test_cuttag_spec_from_target(self):
        spec = make_cuttag_spec("h3k27ac_cuttag", "H3K27ac", "CUTTAG")
        assert spec.assay == "CUTTAG"
        assert spec.target == "H3K27ac"
        assert spec.regulatory_role == "activating_enhancer"
        assert spec.expected_rna_direction == 1

    def test_repressive_mark_spec(self):
        spec = make_cuttag_spec("h3k27me3_cuttag", "H3K27me3")
        assert spec.expected_rna_direction == -1
        assert spec.regulatory_role == "polycomb_repression"

    def test_ctcf_has_no_direction(self):
        spec = make_cuttag_spec("ctcf_cutrun", "CTCF", "CUTRUN")
        assert spec.expected_rna_direction is None

    def test_is_epigenomic(self):
        spec = make_cuttag_spec("h3k27ac_cuttag", "H3K27ac")
        assert spec.is_epigenomic()

    def test_is_activating(self):
        spec = make_cuttag_spec("h3k27ac_cuttag", "H3K27ac")
        assert spec.is_activating()

    def test_is_repressive(self):
        spec = make_cuttag_spec("h3k27me3_cuttag", "H3K27me3")
        assert spec.is_repressive()

    def test_unknown_target_fallback(self):
        info = get_cuttag_target_info("UNKNOWN_MARK")
        assert info["role"] == "unknown"
        assert info["expected_rna_direction"] is None


class TestCUTTAGRegistry:
    def test_all_keys_have_role(self):
        for target, info in CUTTAG_REGISTRY.items():
            assert "role" in info
            assert "expected_rna_direction" in info

    def test_activating_marks_have_positive_direction(self):
        activating = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3"]
        for t in activating:
            assert CUTTAG_REGISTRY[t]["expected_rna_direction"] == 1

    def test_repressive_marks_have_negative_direction(self):
        repressive = ["H3K27me3", "H3K9me3"]
        for t in repressive:
            assert CUTTAG_REGISTRY[t]["expected_rna_direction"] == -1


class TestValidateCUTTAGFeatures:
    def test_valid_features(self):
        df = pd.DataFrame({
            "feature_id": ["chr1:100-200|H3K27ac"],
            "chr": ["chr1"], "start": [100], "end": [200],
            "assay": ["CUTTAG"], "target": ["H3K27ac"],
        })
        issues = validate_cuttag_features(df)
        assert len(issues) == 0

    def test_missing_columns(self):
        df = pd.DataFrame({"feature_id": ["x"]})
        issues = validate_cuttag_features(df)
        assert len(issues) > 0


class TestMultiModalityData:
    def test_epigenomic_matrices_constructor(self):
        """from_epigenomic_matrices should work with RNA + CUT&Tag."""
        from modes.data import MoDEData

        n = 8
        rna = pd.DataFrame(
            np.random.poisson(100, (n, 5)),
            index=[f"s{i}" for i in range(n)],
            columns=["G1", "G2", "G3", "G4", "G5"],
        )
        epi = pd.DataFrame(
            np.random.poisson(50, (n, 3)),
            index=[f"s{i}" for i in range(n)],
            columns=[
                "chr1:100-200|H3K27ac",
                "chr1:300-500|H3K27ac",
                "chr2:100-300|H3K27ac",
            ],
        )
        epi_feat = pd.DataFrame({
            "feature_id": [
                "chr1:100-200|H3K27ac",
                "chr1:300-500|H3K27ac",
                "chr2:100-300|H3K27ac",
            ],
            "chr": ["chr1", "chr1", "chr2"],
            "start": [100, 300, 100],
            "end": [200, 500, 300],
            "assay": ["CUTTAG"] * 3,
            "target": ["H3K27ac"] * 3,
        })
        meta = pd.DataFrame(
            {"condition": ["ctrl"] * 4 + ["trt"] * 4},
            index=[f"s{i}" for i in range(n)],
        )
        data = MoDEData.from_epigenomic_matrices(
            rna_counts=rna,
            epigenomic_counts=epi,
            epigenomic_features=epi_feat,
            metadata=meta,
            condition_col="condition",
            target="H3K27ac",
        )
        assert data.n_samples == n
        assert data.n_genes == 5
        assert "h3k27ac_cuttag" in data.modalities
        assert "rna" in data.modality_specs
        assert "h3k27ac_cuttag" in data.modality_specs
