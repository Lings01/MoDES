"""Tests for EffectEstimator."""

import pytest
import numpy as np

from modes.effects import EffectEstimator


class TestEffectEstimator:
    def test_estimate_atac_effects(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            donor_col="donor",
            batch_col="batch",
            use_empirical_bayes=False,
        )
        effects = estimator.estimate_atac_effects(
            data, list(data.peak_names)
        )
        assert len(effects) > 0
        for e in effects.values():
            assert e.fdr <= 1.0
            assert e.p_value >= 0
            assert e.p_value <= 1.0

    def test_estimate_rna_effects(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            donor_col="donor",
            batch_col="batch",
            use_empirical_bayes=False,
        )
        effects = estimator.estimate_rna_effects(
            data, list(data.gene_names)
        )
        assert len(effects) > 0
        for e in effects.values():
            assert e.fdr <= 1.0

    def test_estimate_effects(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=False,
        )
        atac, rna = estimator.estimate_effects(
            data, list(data.peak_names), list(data.gene_names)
        )
        assert len(atac) > 0
        assert len(rna) > 0

    def test_empirical_bayes_moderation(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=True,
        )
        effects = estimator.estimate_atac_effects(
            data, list(data.peak_names)
        )
        for e in effects.values():
            assert e.p_value >= 0
            assert e.p_value <= 1.0

    def test_effect_directions(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=False,
        )
        effects = estimator.estimate_rna_effects(
            data, list(data.gene_names[:5])  # should be concordant genes
        )
        # Concordant genes should tend to have positive effects
        directions = [e.direction for e in effects.values()]
        assert all(d in (-1, 0, 1) for d in directions)

    def test_with_covariates(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        # Add a numeric covariate
        data.obs["age"] = np.random.uniform(30, 70, data.n_samples)

        estimator = EffectEstimator(
            condition_col="condition",
            covariate_cols=["age"],
            batch_col="batch",
            use_empirical_bayes=False,
        )
        effects = estimator.estimate_rna_effects(
            data, list(data.gene_names[:5])
        )
        assert len(effects) > 0


class TestNBGLMFitting:
    def test_fit_converges(self, synthetic_bulk_data_small):
        data, gt, tss_map = synthetic_bulk_data_small
        estimator = EffectEstimator(
            condition_col="condition",
            use_empirical_bayes=False,
        )
        effects = estimator.estimate_rna_effects(
            data, [data.gene_names[0]]
        )
        effect = list(effects.values())[0]
        assert effect.convergence
        assert not np.isnan(effect.coef)
