"""Tests for statistical utilities."""

import numpy as np
from modes.utils import (
    benjamini_hochberg,
    empirical_bayes_moderate,
    compute_quality_score,
)


class TestBenjaminiHochberg:
    def test_basic(self):
        pvalues = np.array([0.001, 0.01, 0.05, 0.5, 0.9])
        qvalues = benjamini_hochberg(pvalues)
        assert len(qvalues) == len(pvalues)
        assert qvalues[0] <= qvalues[-1]
        assert (qvalues >= 0).all()
        assert (qvalues <= 1).all()

    def test_all_significant(self):
        pvalues = np.array([0.001, 0.002, 0.003])
        qvalues = benjamini_hochberg(pvalues)
        assert qvalues[0] <= 0.05

    def test_all_null(self):
        pvalues = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        qvalues = benjamini_hochberg(pvalues)
        assert (qvalues > 0.1).all()

    def test_empty(self):
        qvalues = benjamini_hochberg(np.array([]))
        assert len(qvalues) == 0


class TestEmpiricalBayes:
    def test_moderate_shrinks_large_ses(self):
        coefs = np.array([1.0] * 20 + [0.0] * 80)
        ses = np.concatenate([
            np.random.uniform(0.1, 0.3, 20),
            np.random.uniform(0.5, 2.0, 80),
        ])
        mod_ses = empirical_bayes_moderate(coefs, ses, df=3)
        assert len(mod_ses) == len(ses)
        # Moderated SEs should be less extreme than raw
        assert np.std(mod_ses) <= np.std(ses) * 1.1

    def test_moderate_with_nans(self):
        coefs = np.array([1.0, 0.5, np.nan, 0.2, 0.0])
        ses = np.array([0.2, 0.3, 0.4, 0.5, np.inf])
        mod_ses = empirical_bayes_moderate(coefs, ses)
        assert len(mod_ses) == 5


class TestQualityScore:
    def test_good_quality(self):
        counts = np.random.poisson(50, 100)
        score = compute_quality_score(counts, None)
        assert 0 <= score <= 1
        assert score > 0.5  # high counts should have good quality

    def test_poor_quality(self):
        counts = np.array([0] * 90 + [1] * 10)
        score = compute_quality_score(counts, None)
        assert score < 0.5

    def test_empty(self):
        score = compute_quality_score(np.array([]), None)
        assert score == 0.0
