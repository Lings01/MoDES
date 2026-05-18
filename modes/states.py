"""Event state classification: grammar-based multi-modal scoring.

v2.0: All states defined via StateRule grammar. All applicable rules are scored
simultaneously and the best-matching state is selected by assignment_score.
State p-values are derived from the modalities that trigger each state.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

from modes._types import EventEvidence, ModalityEffect
from modes.utils import compute_quality_score


class EvidenceBuilder:
    """Construct evidence vectors D_e from effect estimates."""

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
        extra_cols = []
        for mod_name in self.extra_effects:
            for suffix in ["_z", "_fdr", "_pval", "_direction", "_coef"]:
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
                rna_after = ModalityEffect(coef=0, se=1e6, z_score=0, p_value=1.0)

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
                "atac_coef": atac_eff.coef,
                "rna_coef": rna_eff.coef,
            }
            # v2.0: add extra modality evidence columns
            for mod_name, eff_dict in self.extra_effects.items():
                spec = data.modality_specs.get(mod_name)
                mod_eff = None
                if spec and spec.assay == "PROTEIN":
                    # Protein-to-gene links must be explicitly provided.
                    links = getattr(data, 'protein_gene_links', None)
                    if links is not None and len(links) > 0:
                        matched = links[links["gene"] == gene]
                        if len(matched) == 0:
                            gene_short = str(gene).split(":")[0]
                            matched = links[links["gene"].astype(str) == gene_short]
                        for pid in matched["protein_id"].values:
                            mod_eff = eff_dict.get(str(pid))
                            if mod_eff:
                                break
                elif spec and spec.feature_type == "region":
                    feature = peak
                    mod_eff = eff_dict.get(feature)
                    if mod_eff is None:
                        from modes.utils import interval_overlap
                        best_ov = 0.0
                        for k, v in eff_dict.items():
                            ov = interval_overlap(str(feature), str(k))
                            if ov and ov["match"] and ov["min_reciprocal_overlap"] > best_ov:
                                mod_eff = v
                                best_ov = ov["min_reciprocal_overlap"]
                        # Legacy fallback: string match
                        if mod_eff is None:
                            for k, v in eff_dict.items():
                                base = k.split("|")[0] if "|" in str(k) else str(k)
                                if base == str(feature).split("|")[0]:
                                    mod_eff = v
                                    break
                else:
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
                    record[f"{mod_name}_coef"] = mod_eff.coef
                else:
                    record[f"{mod_name}_z"] = 0.0
                    record[f"{mod_name}_fdr"] = 1.0
                    record[f"{mod_name}_pval"] = 1.0
                    record[f"{mod_name}_direction"] = 0
                    record[f"{mod_name}_coef"] = np.nan
            records.append(record)

        return pd.DataFrame(records)


class StateClassifier:
    """Classify events by scoring all applicable StateRules simultaneously.

    Parameters
    ----------
    fdr_threshold : float
        FDR threshold for significance calls (default 0.1).
    quality_threshold : float
        Below this, artifact_risk becomes high (default 0.3).
    use_empirical_bayes : bool
        Apply EB smoothing. Default False — EB is experimental only.
    modality_specs : dict
        Modality specifications for dynamic rule selection.
    """

    def __init__(
        self,
        fdr_threshold: float = 0.1,
        quality_threshold: float = 0.3,
        use_empirical_bayes: bool = False,
        modality_specs: dict | None = None,
    ):
        self.fdr_threshold = fdr_threshold
        self.quality_threshold = quality_threshold
        self.use_empirical_bayes = use_empirical_bayes
        self.modality_specs = modality_specs or {}
        self._rules_cache: list | None = None

    @property
    def rules(self) -> list:
        """Lazy-load applicable StateRules from grammar module."""
        if self._rules_cache is not None:
            return self._rules_cache
        from modes.modalities.state_rules import get_all_rules

        has_epi = any(
            hasattr(s, 'is_epigenomic') and s.is_epigenomic()
            for s in self.modality_specs.values()
        )
        has_protein = any(
            s.assay == "PROTEIN" for s in self.modality_specs.values()
        )
        has_spatial = any(
            s.assay == "SPATIAL" for s in self.modality_specs.values()
        )
        self._rules_cache = get_all_rules(
            include_epigenomic=has_epi,
            include_protein=has_protein,
            include_spatial=has_spatial,
        )
        return self._rules_cache

    @classmethod
    def get_applicable_states(cls, modality_specs: dict | None = None) -> set:
        """Return all applicable biological states for the given modalities."""
        from modes.modalities.state_rules import get_all_rules

        has_epi = any(
            hasattr(s, 'is_epigenomic') and s.is_epigenomic()
            for s in (modality_specs or {}).values()
        )
        has_protein = any(
            s.assay == "PROTEIN" for s in (modality_specs or {}).values()
        )
        has_spatial = any(
            s.assay == "SPATIAL" for s in (modality_specs or {}).values()
        )
        rules = get_all_rules(
            include_epigenomic=has_epi,
            include_protein=has_protein,
            include_spatial=has_spatial,
        )
        return {r.name for r in rules} | {"mixed_evidence", "unresolved"}

    # ── Public API ────────────────────────────────────────────────────

    def classify(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """Classify events by scoring all applicable StateRules.

        Returns DataFrame with: event_id, state, state_assignment_score,
        state_support_pval, state_support_qval, supporting_modalities,
        neutral_modalities, conflicting_modalities, artifact_risk, artifact_reason.
        """
        cols = [
            "event_id", "state", "state_assignment_score",
            "state_support_pval", "state_support_qval",
            "supporting_modalities", "neutral_modalities",
            "conflicting_modalities", "artifact_risk", "artifact_reason",
        ]
        if evidence.empty:
            return pd.DataFrame(columns=cols)

        results = []
        for _, row in evidence.iterrows():
            result = self._classify_event(row)
            results.append(result)

        states_df = pd.DataFrame(results)

        # BH correction on state_support_pval
        from modes.utils import benjamini_hochberg
        pvals = states_df["state_support_pval"].values.astype(float)
        qvals = benjamini_hochberg(pvals)
        states_df["state_support_qval"] = qvals

        return states_df

    # ── Per-event classification ──────────────────────────────────────

    def _classify_event(self, row: pd.Series) -> dict:
        """Score all applicable rules for one event, select best state."""
        eid = row["event_id"]
        quality = float(row.get("quality_score", 0.5))

        # Invalid evidence → unresolved
        if quality <= 0 or np.isnan(quality):
            return self._make_result(eid, "unresolved", np.nan, 1.0,
                                     "", "", "", "low", "")

        # Score all rules
        scores = []
        for rule in self.rules:
            score = self._score_rule(rule, row, quality)
            scores.append(score)

        # Filter to valid scores
        valid_scores = [s for s in scores if s["score"] > 0]
        if not valid_scores:
            return self._make_result(eid, "unresolved", np.nan, 1.0,
                                     "", "", "", "low", "")

        # Sort by assignment_score descending
        valid_scores.sort(key=lambda s: s["score"], reverse=True)
        best = valid_scores[0]

        # Ambiguity check: second-best close to best AND different name.
        # Ties (equal scores, equal n_satisfied) are NOT ambiguous — first in rule order wins.
        if len(valid_scores) > 1:
            second = valid_scores[1]
            close = second["score"] >= best["score"] * 0.85
            better_specificity = second.get("n_satisfied", 0) > best.get("n_satisfied", 0)
            if close and second["name"] != best["name"] and better_specificity:
                return self._make_result(
                    eid, "mixed_evidence", best["score"],
                    best["support_pval"],
                    best["supporting"],
                    best["neutral"],
                    best["conflicting"],
                    "low", "",
                )

        # Compute artifact risk
        art_risk, art_reason = self._compute_artifact_risk(row)

        return self._make_result(
            eid, best["name"], best["score"],
            best["support_pval"],
            best["supporting"],
            best["neutral"],
            best["conflicting"],
            art_risk, art_reason,
        )

    def _score_rule(self, rule, row: pd.Series, quality: float) -> dict:
        """Score one StateRule against evidence for a single event.

        Returns dict with name, score, support_pval, supporting, neutral, conflicting.
        Score of 0 means the rule is invalid (required evidence not met).
        """
        from modes.modalities.state_rules import directed_pvalue

        # Collect modality evidence
        ev = self._modality_evidence_map(row)

        n_required = len(rule.required)
        n_satisfied = 0
        n_conflicts = 0
        supporting = []
        neutral = []
        conflicting = []
        req_pvals = []

        # Check required evidence
        for req in rule.required:
            mod_key = req.modality
            # Map rule modality name to evidence column name
            mod_ev = self._resolve_modality_evidence(ev, req)
            sig = mod_ev["fdr"] < self.fdr_threshold if mod_ev else False
            direction_ok = mod_ev["direction"] == req.direction if mod_ev else False
            coef = mod_ev.get("coef", 0.0) if mod_ev else 0.0
            pval = mod_ev.get("pval", 1.0) if mod_ev else 1.0

            if sig and direction_ok:
                n_satisfied += 1
                dp = directed_pvalue(pval, coef, req.direction)
                req_pvals.append(dp)
                supporting.append(req.modality)
            elif sig and not direction_ok:
                n_conflicts += 1
                conflicting.append(f"{req.modality}(dir_mismatch)")
            else:
                # Required but not significant
                pass

        # Check forbidden evidence
        for fb in rule.forbidden:
            mod_ev = self._resolve_modality_evidence(ev, fb)
            sig = mod_ev["fdr"] < self.fdr_threshold if mod_ev else False
            direction_ok = mod_ev["direction"] == fb.direction if mod_ev else False
            if sig and direction_ok:
                n_conflicts += 1
                conflicting.append(f"{fb.modality}(forbidden)")

        # Check neutral evidence
        for nt in rule.neutral:
            mod_ev = self._resolve_modality_evidence(ev, nt)
            sig = mod_ev["fdr"] < self.fdr_threshold if mod_ev else False
            if sig:
                neutral.append(nt.modality)

        # A rule with no required evidence always matches (null-like)
        if n_required == 0:
            n_satisfied = 1  # dummy to pass validity check

        # Rule is invalid if required evidence not satisfied
        if n_satisfied == 0:
            return {"name": rule.name, "score": 0, "support_pval": 1.0,
                    "supporting": "", "neutral": "", "conflicting": ""}

        # State support p-value: intersection test (max of directed p-values)
        state_support_pval = max(req_pvals) if req_pvals else 1.0

        # Evidence strength: sum of -log10(p) for ALL satisfied required modalities.
        # This naturally differentiates multi-evidence states (e.g., concordant with
        # ATAC↑ + RNA↑) from single-evidence states (chromatin_primed with ATAC↑ only).
        evidence_strength = sum(-np.log10(max(p, 1e-15)) for p in req_pvals) if req_pvals else 0.0
        evidence_strength = min(evidence_strength, 30.0)
        if n_required == 0:
            evidence_strength = max(evidence_strength, 0.1)

        # Penalties
        missing_required = n_required - n_satisfied
        missing_penalty = 0.5 ** missing_required if missing_required > 0 else 1.0
        conflict_penalty = 0.5 ** n_conflicts if n_conflicts > 0 else 1.0

        # Assignment score with specificity bonus for multi-evidence states
        specificity = 1.0 + 0.4 * max(n_satisfied - 1, 0)  # bonus per extra req
        assignment_score = evidence_strength * quality * conflict_penalty * missing_penalty * specificity

        return {
            "name": rule.name,
            "score": float(assignment_score),
            "support_pval": float(state_support_pval),
            "supporting": ";".join(supporting),
            "neutral": ";".join(neutral),
            "conflicting": ";".join(conflicting),
            "n_satisfied": n_satisfied,
        }

    def _resolve_modality_evidence(self, ev: dict, req) -> dict | None:
        """Map a rule modality name to actual evidence columns.

        The rule uses abstract modality names (e.g., 'atac', 'rna', 'protein',
        'cuttag_activating', 'cuttag_repressive', 'spatial'). This method
        resolves them to the actual column prefixes in the evidence DataFrame.
        """
        mod_name = req.modality

        # Direct match: atac, rna are always present
        if mod_name in ("atac", "rna"):
            return ev.get(mod_name)

        # Protein: match any modality with assay == "PROTEIN"
        if mod_name == "protein":
            for key, val in ev.items():
                if key not in ("atac", "rna"):
                    spec = self.modality_specs.get(key)
                    if spec and spec.assay == "PROTEIN":
                        return val
            return None

        # CUT&Tag activating: match modalities with expected_rna_direction == 1
        if mod_name == "cuttag_activating":
            for key, val in ev.items():
                if key not in ("atac", "rna", "protein"):
                    spec = self.modality_specs.get(key)
                    if spec and getattr(spec, 'expected_rna_direction', None) == 1:
                        return val
            return None

        # CUT&Tag repressive: match modalities with expected_rna_direction == -1
        if mod_name == "cuttag_repressive":
            for key, val in ev.items():
                if key not in ("atac", "rna", "protein"):
                    spec = self.modality_specs.get(key)
                    if spec and getattr(spec, 'expected_rna_direction', None) == -1:
                        return val
            return None

        # Spatial: match modalities with assay == "SPATIAL"
        if mod_name == "spatial":
            for key, val in ev.items():
                spec = self.modality_specs.get(key)
                if spec and spec.assay == "SPATIAL":
                    return val
            return None

        # Exact modality name match (fallback)
        return ev.get(mod_name)

    def _modality_evidence_map(self, row: pd.Series) -> dict:
        """Extract modality evidence from a row into a lookup dict.

        Returns dict like:
          {'atac': {'fdr': 0.01, 'direction': 1, 'pval': 0.005, 'coef': 1.2}, ...}
        """
        ev = {}
        # Always available
        ev["atac"] = {
            "fdr": float(row.get("atac_fdr", 1.0)),
            "direction": int(row.get("atac_direction", 0)),
            "pval": float(row.get("atac_pval", 1.0)),
            "coef": float(row.get("atac_coef", 0.0)),
        }
        ev["rna"] = {
            "fdr": float(row.get("rna_fdr", 1.0)),
            "direction": int(row.get("rna_direction", 0)),
            "pval": float(row.get("rna_pval", 1.0)),
            "coef": float(row.get("rna_coef", 0.0)),
        }
        # Extra modalities: detect from column suffixes
        extra_mods = set()
        for col in row.index:
            if col.endswith("_z") and col not in ("z_atac", "z_rna", "z_rna_given_atac"):
                extra_mods.add(col[:-2])  # strip _z
        for mod_name in extra_mods:
            ev[mod_name] = {
                "fdr": float(row.get(f"{mod_name}_fdr", 1.0)),
                "direction": int(row.get(f"{mod_name}_direction", 0)),
                "pval": float(row.get(f"{mod_name}_pval", 1.0)),
                "coef": float(row.get(f"{mod_name}_coef", 0.0)),
            }
        return ev

    def _compute_artifact_risk(self, row: pd.Series) -> tuple:
        """Compute artifact risk level and reason."""
        reasons = []
        qs = float(row.get("quality_score", 1.0))
        atac_sig = float(row.get("atac_fdr", 1.0)) < self.fdr_threshold
        rna_sig = float(row.get("rna_fdr", 1.0)) < self.fdr_threshold
        z_atac = abs(float(row.get("z_atac", 0.0)))
        z_rna = abs(float(row.get("z_rna", 0.0)))

        if qs < self.quality_threshold:
            reasons.append("low_quality_score")
        elif qs < self.quality_threshold * 2:
            reasons.append("borderline_quality")
        if z_atac < 0.5 and atac_sig:
            reasons.append("low_atac_depth")
        if z_rna < 0.5 and rna_sig:
            reasons.append("low_rna_depth")
        if qs < self.quality_threshold * 2 and (atac_sig ^ rna_sig):
            reasons.append("single_modality_low_quality")

        if not reasons:
            return "low", ""
        if "single_modality_low_quality" in reasons:
            return "high", ";".join(reasons)
        return "medium", ";".join(reasons)

    @staticmethod
    def _make_result(eid, state, score, support_pval,
                     supporting, neutral, conflicting,
                     art_risk, art_reason):
        return {
            "event_id": eid,
            "state": state,
            "state_assignment_score": float(score) if not np.isnan(score) else np.nan,
            "state_support_pval": float(support_pval),
            "state_support_qval": 1.0,  # filled later via BH
            "supporting_modalities": supporting,
            "neutral_modalities": neutral,
            "conflicting_modalities": conflicting,
            "artifact_risk": art_risk,
            "artifact_reason": art_reason,
        }


def summarize_states(states: pd.DataFrame) -> pd.DataFrame:
    """Return per-state count summary."""
    summary = states.groupby("state").size().reset_index(name="count")
    summary["fraction"] = summary["count"] / summary["count"].sum()
    return summary.sort_values("count", ascending=False)
