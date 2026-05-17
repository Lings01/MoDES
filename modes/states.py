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

    D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality_score, ...extra modalities]
    """

    def __init__(self, batch_col: str | None = None,
                 extra_modality_effects: dict[str, dict] | None = None):
        self.batch_col = batch_col
        self.extra_effects = extra_modality_effects or {}

    def build(
        self,
        events: pd.DataFrame,
        atac_effects: dict[str, ModalityEffect],
        rna_effects: dict[str, ModalityEffect],
        conditional_effects: pd.DataFrame,
        data,
    ) -> pd.DataFrame:
        records = []
        # Build base + extra modality column names
        extra_cols = []
        for mod_name in self.extra_effects:
            for suffix in ["_z", "_fdr", "_pval", "_direction"]:
                extra_cols.append(f"{mod_name}{suffix}")

        if "event_id" not in conditional_effects.columns:
            return pd.DataFrame(columns=[
                "event_id", "context", "z_atac", "z_rna",
                "z_rna_given_atac", "quality_score",
                "atac_fdr", "rna_fdr", "atac_direction",
                "rna_direction", "atac_pval", "rna_pval",
                "rna_after_atac_pval",
            ] + extra_cols)

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
                atac_qs = atac_qc.get("quality_score", 0.5)
            else:
                atac_qs = 0.5
            if gene in data.rna.columns:
                rna_qc = compute_quality_score(data.rna[gene].values, batch_labels)
                rna_qs = rna_qc.get("quality_score", 0.5)
            else:
                rna_qs = 0.5
            quality = 0.6 * atac_qs + 0.4 * rna_qs

            evidence = EventEvidence(
                event_id=eid,
                context=event.get("context", ""),
                effect_atac=atac_eff,
                effect_rna=rna_eff,
                effect_rna_given_atac=rna_after,
                quality_score=quality,
            )

            record = {
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
            }
            # v2.0: add extra modality evidence columns
            for mod_name, eff_dict in self.extra_effects.items():
                spec = data.modality_specs.get(mod_name)
                mod_eff = None

                # Determine feature match strategy
                if spec and spec.assay == "PROTEIN":
                    # Use protein_gene_links to map gene → protein feature
                    links = getattr(data, 'protein_gene_links', None)
                    if links is not None and len(links) > 0:
                        # Try exact gene match first, then short-name match
                        matched = links[links["gene"] == gene]
                        if len(matched) == 0:
                            gene_short = str(gene).split(":")[0]
                            matched = links[links["gene"].astype(str) == gene_short]
                        for pid in matched["protein_id"].values:
                            mod_eff = eff_dict.get(str(pid))
                            if mod_eff:
                                break
                    # Fallback: try fuzzy match across all effect keys
                    if mod_eff is None:
                        for k, v in eff_dict.items():
                            if str(k) == str(gene).split(":")[0]:
                                mod_eff = v
                                break
                elif spec and spec.feature_type == "region":
                    feature = peak
                    mod_eff = eff_dict.get(feature)
                    # Try fuzzy match: strip |suffix from feature keys (CUT&Tag naming)
                    if mod_eff is None:
                        for k, v in eff_dict.items():
                            base = k.split("|")[0] if "|" in str(k) else str(k)
                            if base == str(feature).split("|")[0]:
                                mod_eff = v
                                break
                else:
                    # Gene-like modality: match by gene name
                    feature = gene
                    mod_eff = eff_dict.get(feature)
                    if mod_eff is None:
                        for k, v in eff_dict.items():
                            if str(k) == str(feature).split(":")[0]:
                                mod_eff = v
                                break

                if mod_eff and np.isfinite(mod_eff.coef):
                    record[f"{mod_name}_z"] = mod_eff.z_score
                    record[f"{mod_name}_fdr"] = mod_eff.fdr
                    record[f"{mod_name}_pval"] = mod_eff.p_value
                    record[f"{mod_name}_direction"] = mod_eff.direction
                else:
                    record[f"{mod_name}_z"] = 0.0
                    record[f"{mod_name}_fdr"] = 1.0
                    record[f"{mod_name}_pval"] = 1.0
                    record[f"{mod_name}_direction"] = 0
            records.append(record)

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

    # v2.0: Core RNA+ATAC states (always available)
    BIOLOGICAL_STATES = {
        "null",
        "concordant",
        "chromatin_primed",
        "rna_only",
        "discordant_opposite",
    }

    # v2.0: Extended states for multi-modal
    EPIGENOMIC_STATES = {
        "epigenomic_concordant", "active_enhancer_primed", "mark_only",
        "repressive_concordant", "derepression", "repressive_primed",
    }
    PROTEIN_STATES = {
        "full_activation", "protein_buffered", "protein_memory", "protein_opposite",
    }
    SPATIAL_STATES = {
        "spatial_region_specific", "spatial_niche_driven",
        "cell_intrinsic", "spatial_edge_artifact",
    }

    ARTIFACT_RISK_LEVELS = ["low", "medium", "high"]

    @classmethod
    def get_applicable_states(cls, modality_specs: dict | None = None) -> set:
        """Return all applicable biological states for the given modalities.

        Delegates to the grammar module (modes.modalities.grammar) when available,
        with static state sets as fallback.
        """
        states = set(cls.BIOLOGICAL_STATES)
        if modality_specs:
            try:
                from modes.modalities.grammar import get_all_states
                grammar_states = get_all_states(
                    include_epigenomic=any(
                        hasattr(s, 'is_epigenomic') and s.is_epigenomic()
                        for s in modality_specs.values()
                    ),
                    include_protein=any(
                        s.assay == "PROTEIN" for s in modality_specs.values()
                    ),
                    include_spatial=any(
                        s.assay == "SPATIAL" for s in modality_specs.values()
                    ),
                )
                return set(grammar_states)
            except ImportError:
                pass
            # Fallback: static state sets
            for spec in modality_specs.values():
                if hasattr(spec, 'is_epigenomic') and spec.is_epigenomic():
                    states |= cls.EPIGENOMIC_STATES
                if spec.assay == "PROTEIN":
                    states |= cls.PROTEIN_STATES
        return states

    def __init__(
        self,
        fdr_threshold: float = 0.1,
        z_threshold: float = 1.0,
        quality_threshold: float = 0.3,
        use_empirical_bayes: bool = True,
        modality_specs: dict | None = None,
    ):
        self.fdr_threshold = fdr_threshold
        self.z_threshold = z_threshold
        self.quality_threshold = quality_threshold
        self.use_empirical_bayes = use_empirical_bayes
        self.modality_specs = modality_specs or {}

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

        # Validate: all states must be in the applicable set
        valid_states = self.get_applicable_states(self.modality_specs)
        for s in states["state"]:
            if s not in valid_states:
                raise ValueError(f"Invalid biological state returned: {s}")

        return states

    def _rule_based_classify(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """First-pass rule-based classification with artifact risk (v2.0 multi-modal)."""
        results = []

        # Detect available extra modalities from evidence columns
        extra_mods = [c.replace("_z", "") for c in evidence.columns
                      if c.endswith("_z") and c not in ("z_atac", "z_rna", "z_rna_given_atac")]

        for _, row in evidence.iterrows():
            atac_sig = row["atac_fdr"] < self.fdr_threshold
            rna_sig = row["rna_fdr"] < self.fdr_threshold
            same_dir = row.get("atac_direction", 0) == row.get("rna_direction", 0)

            # v2.0: Check epigenomic activating marks (H3K27ac, H3K4me3, etc.)
            state = None  # defer; epigenomic states take priority
            for mod_name in extra_mods:
                epi_sig = row.get(f"{mod_name}_fdr", 1.0) < self.fdr_threshold
                epi_dir = row.get(f"{mod_name}_direction", 0)
                if not epi_sig:
                    continue

                # Determine if activating or repressive from modality_specs
                spec = self.modality_specs.get(mod_name)
                is_activating = spec and getattr(spec, 'expected_rna_direction', None) == 1
                is_repressive = spec and getattr(spec, 'expected_rna_direction', None) == -1

                if is_activating:
                    rna_same = epi_dir == row.get("rna_direction", 0) if rna_sig else False
                    if epi_sig and rna_sig and rna_same:
                        state = "epigenomic_concordant"
                    elif epi_sig and not rna_sig:
                        state = "active_enhancer_primed"
                    break
                elif is_repressive:
                    if epi_sig and rna_sig and epi_dir != row.get("rna_direction", 0):
                        state = "repressive_concordant"
                    elif epi_sig and rna_sig and epi_dir == row.get("rna_direction", 0):
                        state = "derepression"
                    elif epi_sig and not rna_sig:
                        state = "repressive_primed"
                    break

            # v2.0: Check protein modalities (if no epigenomic state assigned)
            if state is None:
                for mod_name in extra_mods:
                    spec = self.modality_specs.get(mod_name)
                    if spec is None or spec.assay != "PROTEIN":
                        continue
                    prot_sig = row.get(f"{mod_name}_fdr", 1.0) < self.fdr_threshold
                    prot_dir = row.get(f"{mod_name}_direction", 0)
                    if not prot_sig and not rna_sig and not atac_sig:
                        continue
                    rna_dir = row.get("rna_direction", 0)
                    atac_dir = row.get("atac_direction", 0)
                    if atac_sig and rna_sig and prot_sig and rna_dir == atac_dir and prot_dir == rna_dir:
                        state = "full_activation"
                    elif rna_sig and prot_sig and prot_dir != rna_dir:
                        state = "protein_opposite"
                    elif rna_sig and not prot_sig:
                        state = "protein_buffered"
                    elif not rna_sig and prot_sig:
                        state = "protein_memory"
                    break
                # If no protein state assigned, check for mark_only (epigenomic mark without RNA/ATAC)
                if state is None:
                    for mod_name in extra_mods:
                        spec = self.modality_specs.get(mod_name)
                        is_epi = spec and spec.is_epigenomic() if hasattr(spec, 'is_epigenomic') else False
                        if is_epi and not rna_sig and not atac_sig:
                            epi_sig = row.get(f"{mod_name}_fdr", 1.0) < self.fdr_threshold
                            if epi_sig:
                                state = "mark_only"
                                break

            # v2.0: Check spatial modalities (if no state assigned from molecular layers)
            if state is None:
                for mod_name in extra_mods:
                    spec = self.modality_specs.get(mod_name)
                    if spec is None or spec.assay != "SPATIAL":
                        continue
                    # Spatial evidence columns: neighbor_z, moran_z, edge_score
                    neighbor_z = abs(float(row.get(f"{mod_name}_neighbor_z", 0.0)))
                    moran_z = abs(float(row.get(f"{mod_name}_moran_z", 0.0)))
                    edge_score = float(row.get(f"{mod_name}_edge_score", 0.0))
                    if edge_score > 0.8:
                        state = "spatial_edge_artifact"
                    elif moran_z > 2.0 and neighbor_z > 2.0:
                        state = "spatial_niche_driven"
                    elif moran_z > 2.0 and neighbor_z <= 2.0:
                        state = "spatial_region_specific"
                    elif moran_z <= 2.0:
                        state = "cell_intrinsic"
                    break

            # Fall back to RNA+ATAC rules if no state was assigned from extra modalities
            if state is None:
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
        # v2.0: include extra modality z-scores in evidence vector space
        extra_z_cols = [c for c in evidence.columns
                       if c.endswith("_z") and c not in ev_cols]
        ev_cols = ev_cols + extra_z_cols

        valid = evidence["quality_score"] > 0
        valid &= np.isfinite(evidence[ev_cols].values).all(axis=1)

        if valid.sum() < 10:
            return states

        ev = evidence.loc[valid, ev_cols].values
        state_labels = states.loc[valid, "state"].values

        # EB states: all biological categories found in the data (not just RA core)
        eb_states = sorted(set(state_labels) - {"null"})
        if "null" in set(state_labels):
            eb_states.append("null")
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
