"""Effect size estimation for ATAC and RNA modalities."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
from scipy import stats as scipy_stats

from modes._types import ModalityEffect
from modes.utils import benjamini_hochberg, empirical_bayes_moderate


class EffectEstimator:
    """
    Estimate condition effects for ATAC and RNA modalities using NB GLM.

    Parameters
    ----------
    condition_col : str
        Column in data.obs specifying the condition of interest.
    covariate_cols : list of str, optional
        Additional covariates to include in GLMs.
    donor_col : str, optional
        Column for donor/replicate. Treated as fixed effect with EB shrinkage.
    batch_col : str, optional
        Column for batch. Treated as fixed effect.
    use_empirical_bayes : bool
        Apply empirical Bayes variance moderation. Default True.
    """

    def __init__(
        self,
        condition_col: str,
        covariate_cols: list[str] | None = None,
        donor_col: str | None = None,
        batch_col: str | None = None,
        use_empirical_bayes: bool = True,
        contrast: tuple | None = None,
        allow_poisson_fallback: bool = True,
        allow_simplified_fallback: bool = False,
        cov_type: str = "nonrobust",
    ):
        self.condition_col = condition_col
        self.covariate_cols = covariate_cols or []
        self.donor_col = donor_col
        self.batch_col = batch_col
        self.use_empirical_bayes = use_empirical_bayes
        self.contrast = contrast
        self.allow_poisson_fallback = allow_poisson_fallback
        self.allow_simplified_fallback = allow_simplified_fallback
        self.cov_type = cov_type

    def estimate_atac_effects(
        self,
        data,
        peak_names: list[str],
    ) -> dict[str, ModalityEffect]:
        """
        Estimate condition effects for ATAC peaks.

        Parameters
        ----------
        data : MoDEData
        peak_names : list of str
            Peak names to estimate effects for.

        Returns
        -------
        dict mapping peak_id -> ModalityEffect
        """
        effects = {}
        _, atac_ls = data.get_library_sizes()
        X_base = self._build_design_matrix(data, offset_col=None)

        for peak in peak_names:
            if peak not in data.atac.columns:
                continue
            y = data.atac[peak].values.astype(float)
            # P0 opt: skip low-information features
            if _skip_feature(y):
                effects[peak] = ModalityEffect(
                    coef=0.0, se=1e6, z_score=0.0, p_value=1.0,
                    convergence=False,
                    model_summary={
                        "model_used": "skipped_low_count",
                        "warning": "Below detection threshold",
                        "converged": False,
                        "dropped_covariates": False,
                    },
                )
                continue
            effect = self._fit_nb_glm(
                y, X_base, offset=atac_ls, feature_name=peak
            )
            effects[peak] = effect

        # Apply BH correction across all ATAC effects
        pvals = np.array([e.p_value for e in effects.values()])
        fdrs = benjamini_hochberg(pvals)
        for (key, e), fdr in zip(effects.items(), fdrs):
            e.fdr = float(fdr)
            e.direction = int(np.sign(e.coef)) if np.isfinite(e.coef) and e.coef != 0 else 0

        # EB moderation
        if self.use_empirical_bayes:
            effects = self._apply_eb_moderation_atac(effects, data)

        return effects

    def estimate_rna_effects(
        self,
        data,
        gene_names: list[str],
    ) -> dict[str, ModalityEffect]:
        """
        Estimate condition effects for RNA genes.

        Parameters
        ----------
        data : MoDEData
        gene_names : list of str

        Returns
        -------
        dict mapping gene_name -> ModalityEffect
        """
        effects = {}
        rna_ls, _ = data.get_library_sizes()
        X_base = self._build_design_matrix(data, offset_col=None)

        for gene in gene_names:
            if gene not in data.rna.columns:
                continue
            y = data.rna[gene].values.astype(float)
            if _skip_feature(y):
                effects[gene] = ModalityEffect(
                    coef=0.0, se=1e6, z_score=0.0, p_value=1.0,
                    convergence=False,
                    model_summary={
                        "model_used": "skipped_low_count",
                        "warning": "Below detection threshold",
                        "converged": False,
                        "dropped_covariates": False,
                    },
                )
                continue
            effect = self._fit_nb_glm(
                y, X_base, offset=rna_ls, feature_name=gene
            )
            effects[gene] = effect

        # BH correction
        pvals = np.array([e.p_value for e in effects.values()])
        fdrs = benjamini_hochberg(pvals)
        for (key, e), fdr in zip(effects.items(), fdrs):
            e.fdr = float(fdr)
            e.direction = int(np.sign(e.coef)) if np.isfinite(e.coef) and e.coef != 0 else 0

        # EB moderation
        if self.use_empirical_bayes:
            effects = self._apply_eb_moderation_rna(effects, data)

        return effects

    def estimate_effects(
        self,
        data,
        peak_names: list[str],
        gene_names: list[str],
    ) -> tuple[dict[str, ModalityEffect], dict[str, ModalityEffect]]:
        """Run both ATAC and RNA effect estimation."""
        atac_effects = self.estimate_atac_effects(data, peak_names)
        rna_effects = self.estimate_rna_effects(data, gene_names)
        return atac_effects, rna_effects

    def _build_design_matrix(
        self,
        data,
        offset_col: str | None = None,
    ) -> np.ndarray:
        """
        Build design matrix from data.obs.

        Always includes intercept + condition column.
        Adds covariates, donor, and batch as specified.
        """
        obs = data.obs
        n = data.n_samples

        # Condition
        cond = obs[self.condition_col].values
        if np.issubdtype(cond.dtype, np.number):
            cond_col = cond.astype(float)
            ref_label = "0"
            tgt_label = "1"
        else:
            categories = sorted(set(cond))
            if len(categories) != 2:
                raise NotImplementedError(
                    "MoDES v0.1 supports only binary condition. "
                    "Please run pairwise contrasts or implement a contrast matrix."
                )
            # Determine reference and target levels
            if self.contrast is not None:
                ref_label, tgt_label = self.contrast[0], self.contrast[1]
                if ref_label not in categories or tgt_label not in categories:
                    raise ValueError(
                        f"Contrast levels {self.contrast} not found in condition. "
                        f"Available: {categories}"
                    )
            else:
                ref_label = categories[0]
                tgt_label = categories[1]
            cond_col = (cond == tgt_label).astype(float)

        # Start with intercept + condition
        X_cols = [np.ones(n), cond_col]
        contrast_label = f"{self.condition_col}[{tgt_label}_vs_{ref_label}]"
        col_names = ["intercept", contrast_label]

        # Covariates
        for cov in self.covariate_cols:
            if cov in obs.columns:
                vals = obs[cov].values
                if np.issubdtype(vals.dtype, np.number):
                    X_cols.append(vals.astype(float))
                else:
                    # One-hot encode categorical
                    for cat in sorted(set(vals))[1:]:  # skip first as reference
                        X_cols.append((vals == cat).astype(float))
                        col_names.append(f"{cov}_{cat}")
                    continue
                col_names.append(cov)

        # Batch
        if self.batch_col and self.batch_col in obs.columns:
            batches = obs[self.batch_col].values
            for b in sorted(set(batches))[1:]:
                X_cols.append((batches == b).astype(float))
                col_names.append(f"batch_{b}")

        # Donor
        if self.donor_col and self.donor_col in obs.columns:
            donors = obs[self.donor_col].values
            for d in sorted(set(donors))[1:]:
                X_cols.append((donors == d).astype(float))
                col_names.append(f"donor_{d}")

        X = np.column_stack(X_cols)

        # Rank check
        rank = np.linalg.matrix_rank(X)
        if rank < X.shape[1]:
            msg = (
                f"Design matrix is rank deficient: rank={rank}, "
                f"n_cols={X.shape[1]}.\n"
                "Possible causes:\n"
                "  - condition fully confounded with donor\n"
                "  - condition fully confounded with batch\n"
                "  - too many covariates for the sample size\n"
                "  - donor appears in only one condition\n"
                "Ensure each donor/batch has observations in both conditions."
            )
            raise ValueError(msg)

        return X

    def _fit_nb_glm(
        self,
        y: np.ndarray,
        X: np.ndarray,
        offset: np.ndarray = None,
        feature_name: str = "",
    ) -> ModalityEffect:
        """Fit NB GLM for a single feature."""
        result = _safe_fit_nb_glm(
            y, X, offset,
            allow_poisson=self.allow_poisson_fallback,
            allow_simplified=self.allow_simplified_fallback,
            cov_type=self.cov_type,
        )

        def _build_summary(res) -> dict:
            """Extract model diagnostics from result."""
            return {
                "family": getattr(res, "_modes_family", "unknown"),
                "model_used": getattr(res, "_modes_model_used", "unknown"),
                "alpha": getattr(res, "_modes_alpha", None),
                "alpha_estimated": getattr(res, "_modes_alpha_estimated", False),
                "converged": bool(getattr(res, "converged", False)),
                "dropped_covariates": getattr(res, "_modes_dropped_covariates", False),
                "warning": getattr(res, "_modes_warning", ""),
            }

        if result is None:
            return ModalityEffect(
                coef=np.nan, se=np.nan, z_score=0.0, p_value=1.0,
                convergence=False,
                model_summary={
                    "family": "unknown",
                    "model_used": "failed",
                    "alpha": None,
                    "alpha_estimated": False,
                    "converged": False,
                    "dropped_covariates": False,
                    "warning": "GLM did not converge at any level",
                },
            )

        # Condition coefficient is at index 1 (after intercept)
        coef = result.params[1] if len(result.params) > 1 else np.nan
        se = result.bse[1] if len(result.bse) > 1 else np.nan

        if np.isnan(coef) or np.isnan(se) or se <= 0:
            s = _build_summary(result)
            s["converged"] = False
            s["warning"] = "Invalid coefficient or standard error."
            return ModalityEffect(
                coef=float(coef) if not np.isnan(coef) else 0.0,
                se=float(se) if not np.isnan(se) else 1e6,
                z_score=0.0,
                p_value=1.0,
                convergence=False,
                model_summary=s,
            )

        z = coef / se
        pval = 2 * scipy_stats.norm.sf(abs(z))

        return ModalityEffect(
            coef=float(coef),
            se=float(se),
            z_score=float(z),
            p_value=float(pval),
            convergence=bool(getattr(result, "converged", False)),
            model_summary=_build_summary(result),
        )

    def _apply_eb_moderation_atac(
        self,
        effects: dict[str, ModalityEffect],
        data,
    ) -> dict[str, ModalityEffect]:
        """Apply EB variance moderation to ATAC effects."""
        keys = list(effects.keys())
        coefs = np.array([effects[k].coef for k in keys])
        ses = np.array([effects[k].se for k in keys])

        df = max(data.n_samples - min(data.n_samples, self._n_params(data)), 1)
        mod_ses = empirical_bayes_moderate(coefs, ses, df=df)

        for i, key in enumerate(keys):
            if np.isfinite(mod_ses[i]) and mod_ses[i] > 0:
                e = effects[key]
                e.se = float(mod_ses[i])
                e.z_score = float(e.coef / mod_ses[i])
                # Recompute p-value using t-distribution with moderated df
                t_stat = e.z_score
                pval = 2 * scipy_stats.t.sf(abs(t_stat), df + 3)
                e.p_value = float(pval)

        # Re-BH correction
        pvals = np.array([effects[k].p_value for k in keys])
        fdrs = benjamini_hochberg(pvals)
        for key, fdr in zip(keys, fdrs):
            effects[key].fdr = float(fdr)
            e = effects[key]
            e.direction = int(np.sign(e.coef)) if np.isfinite(e.coef) and e.coef != 0 else 0

        return effects

    def _apply_eb_moderation_rna(
        self,
        effects: dict[str, ModalityEffect],
        data,
    ) -> dict[str, ModalityEffect]:
        """Apply EB variance moderation to RNA effects."""
        # Same logic as ATAC
        return self._apply_eb_moderation_atac(effects, data)

    def _n_params(self, data) -> int:
        """Number of parameters in the design matrix."""
        return self._build_design_matrix(data).shape[1]


def _skip_feature(y: np.ndarray, min_nonzero: int = 3, min_total: float = 10.0) -> bool:
    """Return True if feature has too few nonzero samples or too low total count."""
    return (y > 0).sum() < min_nonzero or y.sum() < min_total


def _safe_fit_nb_glm(
    y: np.ndarray,
    X: np.ndarray,
    offset: np.ndarray = None,
    allow_poisson: bool = True,
    allow_simplified: bool = False,
    cov_type: str = "nonrobust",
) -> Any | None:
    """
    Fit NB GLM with fallback handling.

    Returns statsmodels result or None on failure.
    Each successful result is tagged with _modes_* diagnostic attributes.
    """
    try:
        import statsmodels.api as sm

        # Scale offset for numerical stability
        if offset is not None:
            offset_adj = offset - np.mean(offset)
        else:
            offset_adj = np.zeros(len(y))

        # Try NB GLM with fixed/default alpha
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = sm.GLM(
                endog=y,
                exog=X,
                offset=offset_adj,
                family=sm.families.NegativeBinomial(),
            )
            result = model.fit(maxiter=100, disp=0, cov_type=cov_type)

        if getattr(result, "converged", False):
            result._modes_model_used = "nb_default_alpha"
            result._modes_family = "negative_binomial"
            result._modes_alpha = None
            result._modes_alpha_estimated = False
            result._modes_dropped_covariates = False
            return result

        # Fallback 1: NB2 with fixed alpha=1
        model2 = sm.GLM(
            endog=y,
            exog=X,
            offset=offset_adj,
            family=sm.families.NegativeBinomial(alpha=1.0),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result2 = model2.fit(maxiter=100, disp=0, cov_type=cov_type)

        if getattr(result2, "converged", False):
            result2._modes_model_used = "nb_fixed_alpha"
            result2._modes_family = "negative_binomial"
            result2._modes_alpha = 1.0
            result2._modes_alpha_estimated = False
            result2._modes_dropped_covariates = False
            return result2

        # Fallback 2: Poisson (if allowed)
        if allow_poisson:
            model3 = sm.GLM(
                endog=y,
                exog=X,
                offset=offset_adj,
                family=sm.families.Poisson(),
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result3 = model3.fit(maxiter=200, disp=0, cov_type=cov_type)
            result3._modes_model_used = "poisson_fallback"
            result3._modes_family = "poisson"
            result3._modes_alpha = None
            result3._modes_alpha_estimated = False
            result3._modes_dropped_covariates = False
            if not getattr(result3, "converged", False):
                result3._modes_warning = (
                    "Poisson fallback did not converge; "
                    "coefficients are returned with caution."
                )
            return result3
        return None

    except Exception:
        # Final fallback: simplified model (if allowed)
        if not allow_simplified:
            return None
        try:
            import statsmodels.api as sm

            X_simple = X[:, :2]
            model = sm.GLM(
                endog=y,
                exog=X_simple,
                offset=offset - np.mean(offset) if offset is not None else None,
                family=sm.families.NegativeBinomial(alpha=1.0),
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(maxiter=100, disp=0, cov_type=cov_type)
            result._modes_model_used = "nb_simple_fallback"
            result._modes_family = "negative_binomial"
            result._modes_alpha = 1.0
            result._modes_alpha_estimated = False
            result._modes_dropped_covariates = True
            return result
        except Exception:
            return None
