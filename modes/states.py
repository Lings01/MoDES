"""Event state classification: rule-based + empirical Bayes."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from modes._types import EventEvidence, ModalityEffect
from modes.utils import compute_quality_score


class EvidenceBuilder:
    """
    Construct evidence vectors D_e from effect estimates.

    D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality_score]
    """

    def __init__(self, batch_col: str | None = None):
        self.batch_col = batch_col

    def build(
        self,
        events: pd.DataFrame,
        atac_effects: dict[str, ModalityEffect],
        rna_effects: dict[str, ModalityEffect],
        conditional_effects: pd.DataFrame,
        data,
    ) -> pd.DataFrame:
        records = []
        if "event_id" not in conditional_effects.columns:
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

            # Quality scores with components (P0 opt: use new dict API + batch labels)
            bc = self.batch_col or "batch"
            batch_labels = data.obs[bc].values if bc in data.obs.columns else None
            quality = 0.5
            if peak in data.atac.columns:
                atac_qc = compute_quality_score(data.atac[peak].values, batch_labels)
                quality = 0.6 * atac_qc.get("quality_score", 0.5)
            if gene in data.rna.columns:
                rna_qc = compute_quality_score(data.rna[gene].values, batch_labels)
                quality = 0.6 * quality + 0.4 * rna_qc.get("quality_score", 0.5)

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
    Classify events into regulatory states with artifact risk assessment.

    Stage 1: Rule-based classification into biological states.
    Stage 2: Empirical Bayes refinement for confidence scores.

    Biological states: null, concordant, chromatin_primed, rna_only, discordant_opposite
    Artifact risk: low, medium, high

    Parameters
    ----------
    fdr_threshold : float
        FDR threshold for significance calls (default 0.1).
    z_threshold : float
        Minimum |z| for concordant direction check (default 1.0).
    quality_threshold : float
        Below this, artifact_risk becomes high (default 0.3).
    use_empirical_bayes : bool
        Apply EB refinement. Default True.
    """

    BIOLOGICAL_STATES = {
        "null",
        "concordant",
        "chromatin_primed",
        "rna_only",
        "discordant_opposite",
    }

    ARTIFACT_RISK_LEVELS = ["low", "medium", "high"]

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
        """Classify events into biological states with artifact risk."""
        if evidence.empty:
            return pd.DataFrame(columns=[
                "event_id", "state", "state_confidence",
                "artifact_risk", "artifact_reason",
            ])

        # Stage 1: Rule-based
        states = self._rule_based_classify(evidence)

        # Stage 2: Empirical Bayes refinement
        if self.use_empirical_bayes:
            states = self._empirical_bayes_classify(evidence, states)

        # Validate: no artifact_like as primary state
        for s in states["state"]:
            if s not in self.BIOLOGICAL_STATES:
                raise ValueError(f"Invalid biological state returned: {s}")

        return states

    def _rule_based_classify(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """First-pass rule-based classification with artifact risk."""
        results = []

        for _, row in evidence.iterrows():
            atac_sig = row["atac_fdr"] < self.fdr_threshold
            rna_sig = row["rna_fdr"] < self.fdr_threshold
            same_dir = row["atac_direction"] == row["rna_direction"]

            # Determine biological state only (NO artifact_like as state)
            if atac_sig and rna_sig and same_dir:
                state = "concordant"
            elif atac_sig and rna_sig and not same_dir:
                state = "discordant_opposite"
            elif atac_sig and not rna_sig:
                state = "chromatin_primed"
            elif not atac_sig and rna_sig:
                state = "rna_only"
            else:
                state = "null"

            # Compute artifact risk separately
            artifact_risk, artifact_reason = self._compute_artifact_risk(row)

            results.append({
                "event_id": row["event_id"],
                "state": state,
                "state_confidence": 1.0,
                "artifact_risk": artifact_risk,
                "artifact_reason": artifact_reason,
            })

        return pd.DataFrame(results)

    def _compute_artifact_risk(self, row: pd.Series) -> tuple:
        """Compute artifact risk level and reason string."""
        reasons = []
        quality_score = float(row.get("quality_score", 1.0))
        atac_sig = float(row.get("atac_fdr", 1.0)) < self.fdr_threshold
        rna_sig = float(row.get("rna_fdr", 1.0)) < self.fdr_threshold
        z_atac = abs(float(row.get("z_atac", 0.0)))
        z_rna = abs(float(row.get("z_rna", 0.0)))

        # Quality-based
        if quality_score < self.quality_threshold:
            reasons.append("low_quality_score")
        elif quality_score < self.quality_threshold * 2:
            reasons.append("borderline_quality")

        # Depth-based (low z-score with signal suggests weak evidence)
        if z_atac < 0.5 and atac_sig:
            reasons.append("low_atac_depth")
        if z_rna < 0.5 and rna_sig:
            reasons.append("low_rna_depth")

        # Single-modality with low quality
        if quality_score < self.quality_threshold * 2 and (atac_sig ^ rna_sig):
            reasons.append("single_modality_low_quality")

        if not reasons:
            return "low", ""
        if "single_modality_low_quality" in reasons:
            return "high", ";".join(reasons)
        return "medium", ";".join(reasons)

    def _empirical_bayes_classify(
        self,
        evidence: pd.DataFrame,
        states: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Refine classifications using empirical Bayes.

        Fits per-state Gaussian distributions to evidence vectors
        and computes state_confidence scores.
        """
        ev_cols = ["z_atac", "z_rna", "z_rna_given_atac"]

        valid = evidence["quality_score"] > 0
        valid &= np.isfinite(evidence[ev_cols].values).all(axis=1)

        if valid.sum() < 10:
            return states

        ev = evidence.loc[valid, ev_cols].values
        state_labels = states.loc[valid, "state"].values

        # EB states: biological categories only, excluding artifact_like
        eb_states = [s for s in self.BIOLOGICAL_STATES if s in set(state_labels)]
        if len(eb_states) < 2:
            return states

        priors = {}
        means = {}
        covs = {}

        for s in eb_states:
            mask = state_labels == s
            n_s = mask.sum()
            priors[s] = max(n_s / len(state_labels), 0.01)

            if n_s >= 3:
                means[s] = ev[mask].mean(axis=0)
                cov = np.diag(np.var(ev[mask], axis=0, ddof=1).clip(min=0.01))
                covs[s] = cov
            else:
                means[s] = np.zeros(ev.shape[1])
                covs[s] = np.eye(ev.shape[1])

        confidences_all = []
        for i, (_, row) in enumerate(evidence.iterrows()):
            current_state = states.iloc[i]["state"]
            current_risk = states.iloc[i]["artifact_risk"]
            current_reason = states.iloc[i].get("artifact_reason", "")

            if not valid.iloc[i]:
                confidences_all.append({
                    "event_id": row["event_id"],
                    "state_confidence": 1.0,
                    "state": current_state,
                    "artifact_risk": current_risk,
                    "artifact_reason": current_reason,
                })
                continue

            x = row[ev_cols].values.astype(float)
            if np.isnan(x).any():
                confidences_all.append({
                    "event_id": row["event_id"],
                    "state_confidence": 1.0,
                    "state": current_state,
                    "artifact_risk": current_risk,
                    "artifact_reason": current_reason,
                })
                continue

            log_probs = {}
            for s in eb_states:
                try:
                    mvn = scipy_stats.multivariate_normal(
                        mean=means[s], cov=covs[s], allow_singular=True
                    )
                    ll = mvn.logpdf(x)
                    lp = np.log(priors[s])
                    log_probs[s] = ll + lp
                except Exception:
                    log_probs[s] = -np.inf

            max_lp = max(log_probs.values())
            probs = {
                s: np.exp(lp - max_lp) for s, lp in log_probs.items()
            }
            total = sum(probs.values())
            probs = {s: p / max(total, 1e-10) for s, p in probs.items()}

            best_state = max(probs, key=probs.get)
            best_prob = probs[best_state]

            confidences_all.append({
                "event_id": row["event_id"],
                "state_confidence": float(best_prob),
                "state": best_state,
                "artifact_risk": current_risk,
                    "artifact_reason": current_reason,
            })

        result = pd.DataFrame(confidences_all)
        return result


def summarize_states(states: pd.DataFrame) -> pd.DataFrame:
    """Return per-state count summary."""
    summary = states.groupby("state").size().reset_index(name="count")
    summary["fraction"] = summary["count"] / summary["count"].sum()
    return summary.sort_values("count", ascending=False)
