"""Tests for multi-condition and pseudotime engines."""

import numpy as np
from modes.modalities.dynamic import (
    build_contrast_matrix,
    estimate_pseudotime_lag,
)


class TestContrastMatrix:
    def test_build_binary(self):
        cond = np.array(["ctrl", "ctrl", "trt", "trt"])
        X, names = build_contrast_matrix(cond, [("ctrl", "trt")])
        assert X.shape == (4, 2)
        assert "trt_vs_ctrl" in names
        assert X[0, 0] == 1.0  # intercept
        assert X[2, 1] == 1.0  # trt sample

    def test_build_multi(self):
        cond = np.array(["ctrl", "drugA", "drugB", "ctrl", "drugA", "drugB"])
        contrasts = [("ctrl", "drugA"), ("ctrl", "drugB"), ("drugA", "drugB")]
        X, names = build_contrast_matrix(cond, contrasts)
        assert X.shape == (6, 4)
        assert "drugA_vs_ctrl" in names
        assert "drugB_vs_ctrl" in names
        assert "drugB_vs_drugA" in names

    def test_invalid_contrast_raises(self):
        cond = np.array(["ctrl", "trt"])
        try:
            build_contrast_matrix(cond, [("unknown", "trt")])
            assert False, "Should raise"
        except ValueError:
            pass


class TestPseudotimeLag:
    def test_estimate_lag(self):
        rng = np.random.default_rng(42)
        n = 50
        pseudotime = np.sort(rng.uniform(0, 10, n))
        # ATAC peaks early, RNA peaks later (lag=3)
        atac = np.exp(-((pseudotime - 2) ** 2) / 4) + rng.normal(0, 0.1, n)
        rna = np.exp(-((pseudotime - 5) ** 2) / 4) + rng.normal(0, 0.1, n)
        result = estimate_pseudotime_lag(atac, rna, pseudotime, max_lag=8)
        assert "best_lag" in result
        assert "lag_correlation" in result
        assert "lag_profile" in result
        assert isinstance(result["best_lag"], int)
        # ATAC peaks before RNA → lag should be positive
        assert result["best_lag"] > 0

    def test_no_lag(self):
        rng = np.random.default_rng(42)
        n = 30
        pseudotime = np.sort(rng.uniform(0, 10, n))
        atac = np.exp(-((pseudotime - 5) ** 2) / 4) + rng.normal(0, 0.1, n)
        rna = np.exp(-((pseudotime - 5) ** 2) / 4) + rng.normal(0, 0.1, n)
        result = estimate_pseudotime_lag(atac, rna, pseudotime, max_lag=5)
        assert abs(result["best_lag"]) <= 5
