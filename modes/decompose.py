"""Conditional decomposition: determine if RNA effect is explained by local ATAC."""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from modes._types import ModalityEffect
from modes.utils import benjamini_hochberg


class ConditionalDecomposition:
    """
    Determine whether RNA condition effects are explained by local chromatin.

    For each event e = (peak_c, gene_g):
      1. Fit RNA_g ~ Condition + ATAC_cis + Covariates
      2. Compare the condition coefficient to the marginal RNA effect.
      3. Compute attenuation: how much the condition effect shrinks
         after controlling for ATAC.

    Parameters
    ----------
    condition_col : str
    covariate_cols : list of str, optional
    donor_col : str, optional
    batch_col : str, optional
    """

    def __init__(
        self,
        condition_col: str,
        covariate_cols: Optional[List[str]] = None,
        donor_col: Optional[str] = None,
        batch_col: Optional[str] = None,
    ):
        self.condition_col = condition_col
        self.covariate_cols = covariate_cols or []
        self.donor_col = donor_col
        self.batch_col = batch_col

    def decompose(
        self,
        data,
        events: pd.DataFrame,
        atac_effects: Dict[str, ModalityEffect],
        rna_effects: Dict[str, ModalityEffect],
    ) -> pd.DataFrame:
        """
        Compute conditional RNA effects for each event.

        Parameters
        ----------
        data : MoDEData
        events : DataFrame of EventCandidate
        atac_effects : dict peak_id -> ModalityEffect
        rna_effects : dict gene_name -> ModalityEffect

        Returns
        -------
        DataFrame with columns:
            event_id, rna_after_atac_coef, rna_after_atac_se,
            rna_after_atac_z, rna_after_atac_pval, rna_after_atac_fdr,
            attenuation,
            convergence
        """
        records = []
        rna_ls, _ = data.get_library_sizes()

        for _, event in events.iterrows():
            gene = event["gene"]
            peak = event["peak_id"]

            if gene not in data.rna.columns or peak not in data.atac.columns:
                records.append(
                    self._null_record(event["event_id"])
                )
                continue

            result = self._fit_conditional(
                data=data,
                gene=gene,
                peak=peak,
                rna_offset=rna_ls,
                rna_marginal_effect=rna_effects.get(gene),
            )
            result["event_id"] = event["event_id"]
            records.append(result)

        df = pd.DataFrame(records)

        # BH correction on conditional p-values
        if "rna_after_atac_pval" not in df.columns or df["rna_after_atac_pval"].isna().all():
            df["rna_after_atac_fdr"] = 1.0
        else:
            pvals = df["rna_after_atac_pval"].fillna(1.0).values
            fdrs = benjamini_hochberg(pvals)
            df["rna_after_atac_fdr"] = fdrs

        # Compute attenuation
        # attenuation = beta_conditional / beta_marginal
        # If near 0 -> ATAC explains most RNA effect
        # If near 1 -> ATAC doesn't explain RNA effect
        attenuations = []
        for _, row in df.iterrows():
            rid = row["event_id"]
            for _, ev in events.iterrows():
                if ev["event_id"] == rid:
                    rna_eff = rna_effects.get(ev["gene"])
                    break
            else:
                rna_eff = None

            if rna_eff is not None and abs(rna_eff.coef) > 1e-6:
                att = row["rna_after_atac_coef"] / rna_eff.coef
            else:
                att = np.nan
            attenuations.append(att)
        df["attenuation"] = attenuations

        return df

    def _fit_conditional(
        self,
        data,
        gene: str,
        peak: str,
        rna_offset: np.ndarray,
        rna_marginal_effect: Optional[ModalityEffect] = None,
    ) -> dict:
        """
        Fit RNA_g ~ Condition + ATAC_peak + Covariates.

        Returns dict of conditional effect estimates.
        """
        y = data.rna[gene].values.astype(float)
        # Build design matrix outside try to let validation errors propagate
        X_base = self._build_design_matrix(data)

        # Add ATAC as continuous covariate
        atac_vals = data.atac[peak].values.astype(float)
        atac_mean = np.mean(atac_vals)
        if atac_mean > 0:
            atac_vals = atac_vals / atac_mean
        X = np.column_stack([X_base, atac_vals])

        offset_adj = rna_offset - np.mean(rna_offset)

        try:
            import statsmodels.api as sm

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = sm.GLM(
                    endog=y,
                    exog=X,
                    offset=offset_adj,
                    family=sm.families.NegativeBinomial(),
                )
                result = model.fit(maxiter=100, disp=0)

            if not getattr(result, "converged", False):
                model2 = sm.GLM(
                    endog=y,
                    exog=X,
                    offset=offset_adj,
                    family=sm.families.Poisson(),
                )
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = model2.fit(maxiter=200, disp=0)

            converged = getattr(result, "converged", False)

            if converged and len(result.params) > 1:
                coef = result.params[1]
                se = result.bse[1]
                z = coef / se if se > 0 else 0.0
                pval = 2 * scipy_stats.norm.sf(abs(z))

                return {
                    "rna_after_atac_coef": float(coef),
                    "rna_after_atac_se": float(se),
                    "rna_after_atac_z": float(z),
                    "rna_after_atac_pval": float(pval),
                    "convergence": True,
                    "atac_coef": float(result.params[-1]),
                    "model_aic": float(result.aic),
                }
            else:
                return self._null_conditional()

        except (NotImplementedError, ValueError):
            raise
        except Exception:
            return self._null_conditional()

    def _null_record(self, event_id: str) -> dict:
        """Return null result for an event that couldn't be fitted."""
        return {
            "event_id": event_id,
            "rna_after_atac_coef": np.nan,
            "rna_after_atac_se": np.nan,
            "rna_after_atac_z": np.nan,
            "rna_after_atac_pval": 1.0,
            "rna_after_atac_fdr": 1.0,
            "attenuation": np.nan,
            "convergence": False,
            "atac_coef": np.nan,
            "model_aic": np.nan,
        }

    def _null_conditional(self) -> dict:
        """Return null result for failed fit."""
        return {
            "rna_after_atac_coef": np.nan,
            "rna_after_atac_se": np.nan,
            "rna_after_atac_z": np.nan,
            "rna_after_atac_pval": 1.0,
            "convergence": False,
            "atac_coef": np.nan,
            "model_aic": np.nan,
        }

    def _build_design_matrix(self, data) -> np.ndarray:
        """Build design matrix matching EffectEstimator format."""
        obs = data.obs
        n = data.n_samples

        cond = obs[self.condition_col].values
        if np.issubdtype(cond.dtype, np.number):
            cond_col = cond.astype(float)
        else:
            categories = sorted(set(cond))
            if len(categories) != 2:
                raise NotImplementedError(
                    "MoDES v0.1 supports only binary condition. "
                    "Please run pairwise contrasts or implement a contrast matrix."
                )
            cond_col = (cond == categories[1]).astype(float)

        X_cols = [np.ones(n), cond_col]

        for cov in self.covariate_cols:
            if cov in obs.columns:
                vals = obs[cov].values
                if np.issubdtype(vals.dtype, np.number):
                    X_cols.append(vals.astype(float))
                else:
                    for cat in sorted(set(vals))[1:]:
                        X_cols.append((vals == cat).astype(float))
                    continue

        if self.batch_col and self.batch_col in obs.columns:
            batches = obs[self.batch_col].values
            for b in sorted(set(batches))[1:]:
                X_cols.append((batches == b).astype(float))

        if self.donor_col and self.donor_col in obs.columns:
            donors = obs[self.donor_col].values
            for d in sorted(set(donors))[1:]:
                X_cols.append((donors == d).astype(float))

        X = np.column_stack(X_cols)

        # Rank check
        rank = np.linalg.matrix_rank(X)
        if rank < X.shape[1]:
            raise ValueError(
                f"Design matrix is rank deficient: rank={rank}, "
                f"n_cols={X.shape[1]}. "
                "Check confounding among condition, donor, batch, and covariates."
            )

        return X
