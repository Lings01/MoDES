"""Event state classification: rule-based + empirical Bayes."""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from modes._types import EventEvidence, EventState, ModalityEffect
from modes.utils import compute_quality_score


class EvidenceBuilder:
    """
    Construct evidence vectors D_e from effect estimates.

    D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality_score]
    """

    def build(
        self,
        events: pd.DataFrame,
        atac_effects: Dict[str, ModalityEffect],
        rna_effects: Dict[str, ModalityEffect],
        conditional_effects: pd.DataFrame,
        data,
    ) -> pd.DataFrame:
        """
        Build evidence vectors for all events.

        Parameters
        ----------
        events : DataFrame of EventCandidate
        atac_effects : dict peak_id -> ModalityEffect
        rna_effects : dict gene_name -> ModalityEffect
        conditional_effects : DataFrame from ConditionalDecomposition
        data : MoDEData

        Returns
        -------
        DataFrame with evidence vector columns.
        """
        records = []
        if "event_id" not in conditional_effects.columns:
            # No valid conditional effects; create empty evidence
            return pd.DataFrame(columns=[
                "event_id", "context", "z_atac", "z_rna",
                "z_rna_given_atac", "quality_score",
                "atac_fdr", "rna_fdr", "atac_direction",
                "rna_direction", "atac_pval", "rna_pval",
                "rna_after_atac_pval",
            ])

        cond_map = dict(zip(
            conditional_effects["event_id"].values,
            conditional_effects.index,
        ))

        for _, event in events.iterrows():
            eid = event["event_id"]
            gene = event["gene"]
            peak = event["peak_id"]

            atac_eff = atac_effects.get(peak)
            rna_eff = rna_effects.get(gene)

            if atac_eff is None:
                atac_eff = ModalityEffect(coef=0, se=1e6, z_score=0, p_value=1.0)
            if rna_eff is None:
                rna_eff = ModalityEffect(coef=0, se=1e6, z_score=0, p_value=1.0)

            # Conditional effect
            cond_row = None
            if eid in cond_map:
                cond_row = conditional_effects.iloc[cond_map[eid]]
            else:
                matching = conditional_effects[
                    conditional_effects["event_id"] == eid
                ]
                if len(matching) > 0:
                    cond_row = matching.iloc[0]

            if cond_row is not None:
                rna_after = ModalityEffect(
                    coef=cond_row["rna_after_atac_coef"],
                    se=cond_row["rna_after_atac_se"],
                    z_score=cond_row["rna_after_atac_z"],
                    p_value=cond_row["rna_after_atac_pval"],
                    fdr=cond_row["rna_after_atac_fdr"],
                    convergence=bool(cond_row.get("convergence", False)),
                )
            else:
                rna_after = ModalityEffect(
                    coef=0, se=1e6, z_score=0, p_value=1.0
                )

            # Quality score
            if peak in data.atac.columns:
                atac_counts = data.atac[peak].values
                quality = compute_quality_score(atac_counts, None)
            else:
                quality = 0.5

            evidence = EventEvidence(
                event_id=eid,
                context=event.get("context", ""),
                effect_atac=atac_eff,
                effect_rna=rna_eff,
                effect_rna_given_atac=rna_after,
                quality_score=quality,
            )

            records.append({
                "event_id": eid,
                "context": evidence.context,
                "z_atac": evidence.z_atac,
                "z_rna": evidence.z_rna,
                "z_rna_given_atac": evidence.z_rna_given_atac,
                "quality_score": evidence.quality_score,
                "atac_fdr": atac_eff.fdr,
                "rna_fdr": rna_eff.fdr,
                "atac_direction": atac_eff.direction,
                "rna_direction": rna_eff.direction,
                "atac_pval": atac_eff.p_value,
                "rna_pval": rna_eff.p_value,
                "rna_after_atac_pval": rna_after.p_value,
            })

        return pd.DataFrame(records)


class StateClassifier:
    """
    Classify events into regulatory states.

    Stage 1: Rule-based deterministic classification.
    Stage 2: Empirical Bayes refinement with posterior probabilities.

    Parameters
    ----------
    fdr_threshold : float
        FDR threshold for significance calls (default 0.1).
    z_threshold : float
        Minimum |z| for concordant direction check (default 1.0).
    quality_threshold : float
        Quality score below which events are flagged as artifact-like (default 0.3).
    use_empirical_bayes : bool
        Apply EB refinement. Default True.
    """

    # Phase 1 states
    VALID_STATES = [
        "null",
        "concordant",
        "chromatin_primed",
        "rna_only",
        "artifact_like",
    ]

    def __init__(
        self,
        fdr_threshold: float = 0.1,
        z_threshold: float = 1.0,
        quality_threshold: float = 0.3,
        use_empirical_bayes: bool = True,
    ):
        self.fdr_threshold = fdr_threshold
        self.z_threshold = z_threshold
        self.quality_threshold = quality_threshold
        self.use_empirical_bayes = use_empirical_bayes

    def classify(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """
        Classify events into regulatory states.

        Parameters
        ----------
        evidence : DataFrame
            Output of EvidenceBuilder.build().

        Returns
        -------
        DataFrame with columns: event_id, state, posterior_prob
        """
        if evidence.empty:
            return pd.DataFrame(columns=["event_id", "state", "posterior_prob"])

        # Stage 1: Rule-based
        states = self._rule_based_classify(evidence)

        # Stage 2: Empirical Bayes refinement
        if self.use_empirical_bayes:
            states = self._empirical_bayes_classify(evidence, states)

        return states

    def _rule_based_classify(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """First-pass rule-based classification."""
        results = []

        for _, row in evidence.iterrows():
            atac_sig = row["atac_fdr"] < self.fdr_threshold
            rna_sig = row["rna_fdr"] < self.fdr_threshold
            same_dir = row["atac_direction"] == row["rna_direction"]
            low_qual = row["quality_score"] < self.quality_threshold

            if atac_sig and rna_sig and same_dir:
                state = "concordant"
            elif atac_sig and not rna_sig:
                state = "chromatin_primed"
            elif not atac_sig and rna_sig:
                state = "rna_only"
            elif (atac_sig or rna_sig) and low_qual:
                state = "artifact_like"
            else:
                state = "null"

            # Also check direction mismatch cases
            if atac_sig and rna_sig and not same_dir:
                # Opposite directions -> likely artifact
                if low_qual:
                    state = "artifact_like"
                else:
                    # Keep as RNA_only since it could be trans-regulation
                    state = "rna_only"

            results.append({
                "event_id": row["event_id"],
                "state": state,
                "posterior_prob": 1.0,
            })

        return pd.DataFrame(results)

    def _empirical_bayes_classify(
        self,
        evidence: pd.DataFrame,
        states: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Refine classifications using empirical Bayes.

        Fits per-state Gaussian distributions to evidence vectors
        and computes posterior probabilities.
        """
        # Extract evidence vectors (z_atac, z_rna, z_rna_given_atac)
        # Exclude quality_score from Gaussian model since it's non-Gaussian
        ev_cols = ["z_atac", "z_rna", "z_rna_given_atac"]

        # Only use converged events for prior estimation
        valid = evidence["quality_score"] > 0
        valid &= np.isfinite(evidence[ev_cols].values).all(axis=1)

        if valid.sum() < 10:
            return states

        ev = evidence.loc[valid, ev_cols].values
        state_labels = states.loc[valid, "state"].values

        # Fit per-state means and covariances
        unique_states = [s for s in self.VALID_STATES if s in set(state_labels)]
        priors = {}
        means = {}
        covs = {}

        for s in unique_states:
            mask = state_labels == s
            n_s = mask.sum()
            priors[s] = n_s / len(state_labels)

            if n_s >= 3:
                means[s] = ev[mask].mean(axis=0)
                # Diagonal covariance for robustness
                cov = np.diag(np.var(ev[mask], axis=0, ddof=1).clip(min=0.01))
                covs[s] = cov
            else:
                means[s] = np.zeros(ev.shape[1])
                covs[s] = np.eye(ev.shape[1])

        # Compute posteriors for all events
        posteriors_all = []
        for i, (_, row) in enumerate(evidence.iterrows()):
            if not valid.iloc[i]:
                posteriors_all.append({
                    "event_id": row["event_id"],
                    "posterior_prob": 1.0,
                    "state": states.iloc[i]["state"],
                })
                continue

            x = row[ev_cols].values.astype(float)
            if np.isnan(x).any():
                posteriors_all.append({
                    "event_id": row["event_id"],
                    "posterior_prob": 1.0,
                    "state": states.iloc[i]["state"],
                })
                continue

            # Compute log likelihoods
            log_probs = {}
            for s in unique_states:
                try:
                    mvn = scipy_stats.multivariate_normal(
                        mean=means[s], cov=covs[s], allow_singular=True
                    )
                    ll = mvn.logpdf(x)
                    lp = np.log(max(priors[s], 1e-10))
                    log_probs[s] = ll + lp
                except Exception:
                    log_probs[s] = -np.inf

            # Numerical stability: subtract max log prob
            max_lp = max(log_probs.values())
            probs = {
                s: np.exp(lp - max_lp) for s, lp in log_probs.items()
            }
            total = sum(probs.values())
            probs = {s: p / total for s, p in probs.items()}

            # Apply quality penalty for artifact-like state
            quality = row["quality_score"]
            if quality < self.quality_threshold:
                penalty = quality / self.quality_threshold
                if "artifact_like" in probs:
                    probs["artifact_like"] = max(
                        probs.get("artifact_like", 0), 1 - penalty
                    )
                else:
                    probs["artifact_like"] = 1 - penalty
                # Renormalize
                total = sum(probs.values())
                probs = {s: p / total for s, p in probs.items()}

            best_state = max(probs, key=probs.get)
            best_prob = probs[best_state]

            posteriors_all.append({
                "event_id": row["event_id"],
                "posterior_prob": float(best_prob),
                "state": best_state,
            })

        result = pd.DataFrame(posteriors_all)
        return result


def summarize_states(states: pd.DataFrame) -> pd.DataFrame:
    """Return per-state count summary."""
    summary = states.groupby("state").size().reset_index(name="count")
    summary["fraction"] = summary["count"] / summary["count"].sum()
    return summary.sort_values("count", ascending=False)
