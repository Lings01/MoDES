"""MoDES orchestrator and result container.

v2.0: Fixed event_table schema + long-format event_modality_evidence.
All multi-modal effects live in the long-format table, not as dynamic columns.
State classification is grammar-driven with state_support_pval from required modalities.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from modes.data import MoDEData
from modes.decompose import ConditionalDecomposition
from modes.effects import EffectEstimator
from modes.events import EventCandidateBuilder
from modes.states import EvidenceBuilder, StateClassifier


def _event_result_columns():
    """Return fixed event_table columns (v2.0 frozen schema)."""
    return [
        "event_id", "tf_name", "peak_id", "gene", "context",
        "link_source", "link_score",
        "state_family", "state", "state_assignment_score", "raw_state_assignment_score",
        "state_support_score", "state_support_adjusted_score",
        "supporting_modalities", "absent_modalities", "conflicting_modalities",
        "missing_modalities",
        "atac_coef", "atac_se", "atac_pval", "atac_fdr", "atac_direction",
        "rna_coef", "rna_se", "rna_pval", "rna_fdr", "rna_direction",
        "rna_after_atac_coef", "rna_after_atac_se",
        "rna_after_atac_pval", "rna_after_atac_fdr",
        "artifact_risk", "artifact_reason",
        "quality_score",
        # Backward compat: kept for existing consumers
        "state_confidence_deprecated", "event_pval_deprecated", "event_fdr_deprecated",
    ]


def _read_optional_tsv(output_dir: str, filename: str):
    """Read a TSV file if it exists, else return None."""
    path = os.path.join(output_dir, filename)
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    return None


class MoDES:
    """Multi-Omics Discordance/Event State annotation engine.

    Main orchestrator for the MoDES pipeline.

    Parameters
    ----------
    data : MoDEData
        Input data container.
    condition_col : str
        Column in data.obs specifying contrast of interest.
    covariate_cols : list of str, optional
        Additional covariates.
    donor_col : str, optional
        Donor/replicate column.
    batch_col : str, optional
        Batch column.
    fdr_threshold : float
        FDR cutoff for significance. Default 0.1.
    genome_annotation : str, optional
        Path to GTF for TSS annotation.
    external_links : pd.DataFrame, optional
        Pre-computed peak-to-gene links.
    motif_annotation : pd.DataFrame, optional
        Peak-to-TF motif mapping.
    tss_map : dict, optional
        Manual gene -> (name, chr, tss_pos) mapping.
    contrast : tuple, optional
        (reference, target) levels.
    """

    def __init__(
        self,
        data: MoDEData,
        condition_col: str,
        covariate_cols: list[str] | None = None,
        donor_col: str | None = None,
        batch_col: str | None = None,
        fdr_threshold: float = 0.1,
        genome_annotation: str | None = None,
        external_links: pd.DataFrame | None = None,
        motif_annotation: pd.DataFrame | None = None,
        tss_map: dict | None = None,
        contrast: tuple | None = None,
        allow_poisson_fallback: bool = True,
        allow_simplified_fallback: bool = False,
        conditional_mode: str = "auto",
        cov_type: str = "nonrobust",
        min_nonzero_samples: int = 3,
        min_total_count: float = 10.0,
        batch_col_quality: str | None = None,
        use_empirical_bayes: bool = False,
    ):
        self.data = data
        self.condition_col = condition_col
        self.covariate_cols = covariate_cols or []
        self.donor_col = donor_col
        self.batch_col = batch_col
        self.fdr_threshold = fdr_threshold
        self.genome_annotation = genome_annotation
        self.external_links = external_links
        self.motif_annotation = motif_annotation
        self.tss_map = tss_map
        self.contrast = contrast
        self.allow_poisson_fallback = allow_poisson_fallback
        self.allow_simplified_fallback = allow_simplified_fallback
        self.conditional_mode = conditional_mode
        self.cov_type = cov_type
        self.min_nonzero_samples = min_nonzero_samples
        self.min_total_count = min_total_count
        self.batch_col_quality = batch_col_quality
        self.use_empirical_bayes = use_empirical_bayes

        # Pipeline state
        self.events: pd.DataFrame | None = None
        self.effects: dict[str, dict] = {}
        self.atac_effects: dict | None = None
        self.rna_effects: dict | None = None
        self.conditional_effects: pd.DataFrame | None = None
        self.evidence: pd.DataFrame | None = None
        self.states: pd.DataFrame | None = None
        self.results: MoDESResult | None = None
        self.modality_evidence: pd.DataFrame | None = None
        # v2.0: multi-model conditional decomposition results
        self.conditional_models: dict[str, pd.DataFrame] = {}

    def run(self) -> "MoDESResult":
        """Execute the full MoDES pipeline."""
        self.build_events()
        self.estimate_effects()
        self.decompose()
        self.build_evidence()
        self.classify_states()
        self._build_modality_evidence()
        self.results = self._assemble_results()
        return self.results

    def build_events(
        self,
        external_links: pd.DataFrame | None = None,
        motif_annotation: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Step 1: Build candidate regulatory events."""
        builder = EventCandidateBuilder(promoter_window=2000, distal_window=250000)
        ext_links = self.external_links if external_links is None else external_links
        motifs = self.motif_annotation if motif_annotation is None else motif_annotation
        self.events = builder.build(
            gene_names=list(self.data.gene_names),
            peak_names=list(self.data.peak_names),
            external_links=ext_links,
            motif_annotation=motifs,
            genome_annotation=self.genome_annotation,
            tss_map=self.tss_map,
        )
        if self.events is None or len(self.events) == 0:
            raise ValueError(
                "No candidate events were generated. "
                "Provide genome_annotation, tss_map, or external_links."
            )
        return self.events

    def estimate_effects(self) -> tuple[dict, dict]:
        """Step 2: Estimate condition effects for all modalities."""
        if self.events is None:
            raise RuntimeError("Call build_events() first")

        estimator = EffectEstimator(
            condition_col=self.condition_col,
            covariate_cols=self.covariate_cols,
            donor_col=self.donor_col,
            batch_col=self.batch_col,
            use_empirical_bayes=True,
            contrast=self.contrast,
            allow_poisson_fallback=self.allow_poisson_fallback,
            allow_simplified_fallback=self.allow_simplified_fallback,
            cov_type=self.cov_type,
            min_nonzero_samples=self.min_nonzero_samples,
            min_total_count=self.min_total_count,
        )

        peak_names = list(self.events["peak_id"].unique())
        gene_names = list(self.events["gene"].unique())

        self.atac_effects, self.rna_effects = estimator.estimate_effects(
            self.data, peak_names, gene_names
        )
        self.effects["atac"] = self.atac_effects
        self.effects["rna"] = self.rna_effects
        for mod_name in self.data.modalities:
            if mod_name in ("rna", "atac"):
                continue
            feat_names = list(self.data.modalities[mod_name].columns)
            self.effects[mod_name] = estimator.estimate_modality_effects(
                self.data, feat_names, modality_name=mod_name,
            )
        return self.atac_effects, self.rna_effects

    def decompose(self) -> pd.DataFrame:
        """Step 3: Conditional decomposition (RNA after ATAC, plus multi-modal)."""
        if self.atac_effects is None or self.rna_effects is None:
            raise RuntimeError("Call estimate_effects() first")

        decomposer = ConditionalDecomposition(
            condition_col=self.condition_col,
            covariate_cols=self.covariate_cols,
            donor_col=self.donor_col,
            batch_col=self.batch_col,
            contrast=self.contrast,
            conditional_mode=self.conditional_mode,
        )

        self.conditional_effects = decomposer.decompose(
            data=self.data,
            events=self.events,
            atac_effects=self.atac_effects,
            rna_effects=self.rna_effects,
        )
        self.conditional_models["rna_after_atac"] = self.conditional_effects

        # v2.0: Multi-modal conditional decomposition
        from modes._types import (
            ConditionalModelSpec, RNA_AFTER_H3K27AC, RNA_AFTER_ATAC_H3K27AC,
            PROTEIN_AFTER_RNA,
        )
        multi_specs = []
        # Add models for available modalities
        for mod_name, spec in self.data.modality_specs.items():
            if hasattr(spec, 'is_epigenomic') and spec.is_epigenomic():
                if spec.expected_rna_direction == 1:
                    multi_specs.append(RNA_AFTER_H3K27AC)
                    multi_specs.append(RNA_AFTER_ATAC_H3K27AC)
                break
        if any(s.assay == "PROTEIN" for s in self.data.modality_specs.values()):
            multi_specs.append(PROTEIN_AFTER_RNA)

        extra_eff = {k: v for k, v in self.effects.items() if k not in ("rna", "atac")}
        multi_df = decomposer.decompose_multi(
            self.data, self.events, multi_specs,
            self.atac_effects, self.rna_effects, extra_eff,
        )
        for model_name in multi_df["model_name"].unique() if not multi_df.empty else []:
            sub = multi_df[multi_df["model_name"] == model_name]
            self.conditional_models[model_name] = sub

        return self.conditional_effects

    def build_evidence(self) -> pd.DataFrame:
        """Step 4: Build evidence vectors."""
        if self.conditional_effects is None:
            raise RuntimeError("Call decompose() first")

        extra_effects = {k: v for k, v in self.effects.items()
                        if k not in ("rna", "atac")}
        builder = EvidenceBuilder(
            batch_col=self.batch_col_quality or self.batch_col,
            extra_modality_effects=extra_effects if extra_effects else None,
        )
        self.evidence = builder.build(
            events=self.events,
            atac_effects=self.atac_effects,
            rna_effects=self.rna_effects,
            conditional_effects=self.conditional_effects,
            data=self.data,
        )
        return self.evidence

    def classify_states(self) -> pd.DataFrame:
        """Step 5: Grammar-based multi-modal state classification."""
        if self.evidence is None:
            raise RuntimeError("Call build_evidence() first")

        classifier = StateClassifier(
            fdr_threshold=self.fdr_threshold,
            use_empirical_bayes=self.use_empirical_bayes,
            modality_specs=self.data.modality_specs,
        )
        self.states = classifier.classify(self.evidence)
        return self.states

    def _assemble_results(self) -> "MoDESResult":
        """Combine all pipeline outputs into MoDESResult."""
        if self.events is None or len(self.events) == 0:
            params = {
                "condition_col": self.condition_col,
                "fdr_threshold": self.fdr_threshold,
                "n_events": 0,
                "n_samples": self.data.n_samples,
            }
            return MoDESResult(
                event_table=pd.DataFrame(columns=_event_result_columns()),
                params=params,
            )

        # Pre-build lookup maps
        cond_map = {}
        if self.conditional_effects is not None:
            for _, cr in self.conditional_effects.iterrows():
                cond_map[cr["event_id"]] = cr
        state_map = {}
        if self.states is not None:
            for _, sr in self.states.iterrows():
                state_map[sr["event_id"]] = sr

        from modes.utils import benjamini_hochberg

        records = []
        for _, event in self.events.iterrows():
            eid = event["event_id"]
            peak = event["peak_id"]
            gene = event["gene"]

            atac = self.atac_effects.get(peak) if self.atac_effects else None
            rna = self.rna_effects.get(gene) if self.rna_effects else None

            # Conditional effect
            cr = cond_map.get(eid)
            rna_after_coef = cr["rna_after_atac_coef"] if cr is not None else np.nan
            rna_after_se = cr["rna_after_atac_se"] if cr is not None else np.nan
            rna_after_pval = cr["rna_after_atac_pval"] if cr is not None else 1.0
            rna_after_fdr = cr["rna_after_atac_fdr"] if cr is not None else 1.0

            # State
            sr = state_map.get(eid)
            if sr is not None:
                state_family = sr.get("state_family", sr["state"])
                state = sr["state"]
                assignment_score = sr["state_assignment_score"]
                support_score = sr.get("state_support_score", 0.0)
                support_adj = sr.get("state_support_adjusted_score", 0.0)
                supporting = sr.get("supporting_modalities", "")
                absent = sr.get("absent_modalities", "")
                conflicting = sr.get("conflicting_modalities", "")
                missing = sr.get("missing_modalities", "")
                artifact_risk = sr.get("artifact_risk", "low")
                artifact_reason = sr.get("artifact_reason", "")
            else:
                state_family = "unresolved"
                state = "unresolved"
                assignment_score = np.nan
                support_score = 0.0
                support_adj = 0.0
                supporting = ""
                absent = ""
                conflicting = ""
                missing = ""
                artifact_risk = "low"
                artifact_reason = ""

            # Quality
            quality = float(event.get("quality_score", 0.5)) if "quality_score" in event.index else 0.5

            # Link score: confidence in peak-gene link (0-1)
            link_source = str(event.get("link_source", ""))
            link_score = _compute_link_score(link_source, event.get("distance_to_tss", 0))

            # v2.0: assignment_score incorporates link uncertainty per reviewer feedback
            adjusted_assignment = (assignment_score * link_score) if not np.isnan(assignment_score) else np.nan

            rec = {
                "event_id": eid,
                "tf_name": event.get("tf_name"),
                "peak_id": peak,
                "gene": gene,
                "context": event.get("context", ""),
                "link_source": link_source,
                "link_score": link_score,
                "state_family": state_family,
                "state": state,
                "state_assignment_score": adjusted_assignment,
                "raw_state_assignment_score": assignment_score,
                "state_support_score": support_score,
                "state_support_adjusted_score": support_adj,
                "supporting_modalities": supporting,
                "absent_modalities": absent,
                "conflicting_modalities": conflicting,
                "missing_modalities": missing,
                "atac_coef": atac.coef if atac else np.nan,
                "atac_se": atac.se if atac else np.nan,
                "atac_coef": atac.coef if atac else np.nan,
                "atac_se": atac.se if atac else np.nan,
                "atac_pval": atac.p_value if atac else 1.0,
                "atac_fdr": atac.fdr if atac else 1.0,
                "atac_direction": atac.direction if atac else 0,
                "rna_coef": rna.coef if rna else np.nan,
                "rna_se": rna.se if rna else np.nan,
                "rna_pval": rna.p_value if rna else 1.0,
                "rna_fdr": rna.fdr if rna else 1.0,
                "rna_direction": rna.direction if rna else 0,
                "rna_after_atac_coef": rna_after_coef,
                "rna_after_atac_se": rna_after_se,
                "rna_after_atac_pval": rna_after_pval,
                "rna_after_atac_fdr": rna_after_fdr,
                "artifact_risk": artifact_risk,
                "artifact_reason": artifact_reason,
                "quality_score": quality,
                # Deprecated backward-compat aliases
                "state_confidence_deprecated": assignment_score,
                "event_pval_deprecated": np.nan,
                "event_fdr_deprecated": np.nan,
            }
            records.append(rec)

        event_table = pd.DataFrame(records, columns=_event_result_columns())

        # Params
        import sys
        import numpy as _np
        import pandas as _pd
        import statsmodels as _sm
        from modes import __version__ as _modes_ver

        params = {
            "modes_version": _modes_ver,
            "python_version": sys.version.split()[0],
            "numpy_version": _np.__version__,
            "pandas_version": _pd.__version__,
            "statsmodels_version": _sm.__version__,
            "condition_col": self.condition_col,
            "contrast": str(self.contrast) if self.contrast else "auto",
            "fdr_threshold": self.fdr_threshold,
            "use_empirical_bayes": self.use_empirical_bayes,
            "n_events": len(records),
            "n_samples": self.data.n_samples,
            "n_genes": self.data.n_genes,
            "n_peaks": self.data.n_peaks,
            "n_external_links": len(self.external_links) if self.external_links is not None else 0,
        }

        model_diag = self._build_model_diagnostics()

        # Combine multi-model conditional effects
        all_cond = []
        for name, df in self.conditional_models.items():
            if df is not None and not df.empty:
                all_cond.append(df)
        cond_combined = pd.concat(all_cond, ignore_index=True) if all_cond else pd.DataFrame()

        return MoDESResult(
            event_table=event_table,
            state_probabilities=self.states.copy() if self.states is not None else None,
            layer_effects=self._build_layer_effects_df(),
            evidence_vectors=self.evidence.copy() if self.evidence is not None else None,
            model_diagnostics=model_diag,
            modality_evidence=self.modality_evidence.copy() if self.modality_evidence is not None and len(self.modality_evidence) > 0 else None,
            conditional_effects=cond_combined,
            params=params,
        )

    def _build_layer_effects_df(self) -> pd.DataFrame:
        """Assemble per-event layer effect estimates."""
        records = []
        for _, event in self.events.iterrows():
            eid = event["event_id"]
            peak = event["peak_id"]
            gene = event["gene"]
            atac = self.atac_effects.get(peak) if self.atac_effects else None
            rna = self.rna_effects.get(gene) if self.rna_effects else None
            records.append({
                "event_id": eid, "peak_id": peak, "gene": gene,
                "atac_coef": atac.coef if atac else np.nan,
                "atac_se": atac.se if atac else np.nan,
                "atac_z": atac.z_score if atac else np.nan,
                "atac_fdr": atac.fdr if atac else 1.0,
                "rna_coef": rna.coef if rna else np.nan,
                "rna_se": rna.se if rna else np.nan,
                "rna_z": rna.z_score if rna else np.nan,
                "rna_fdr": rna.fdr if rna else 1.0,
            })
        return pd.DataFrame(records)

    def _build_modality_evidence(self) -> pd.DataFrame:
        """Build long-format event x modality evidence table (v2.0).

        Each row includes directed_score and per-modality quality components
        (detection_score, depth_score, batch_score).
        """
        from modes.utils import compute_quality_score, interval_overlap
        from modes.modalities.state_rules import directed_score

        rows = []
        bc = self.batch_col_quality or self.batch_col
        batch_labels = self.data.obs[bc].values if bc in self.data.obs.columns else None

        for _, event in self.events.iterrows():
            eid = event["event_id"]
            gene = event["gene"]
            peak = event["peak_id"]

            def _qc_scores(counts):
                qc = compute_quality_score(counts, batch_labels)
                return (qc.get("detection_score", 0.5),
                        qc.get("depth_score", 0.5),
                        qc.get("batch_score", 1.0),
                        qc.get("quality_score", 0.5))

            # RNA evidence
            rna_eff = self.rna_effects.get(gene) if self.rna_effects else None
            if rna_eff and gene in self.data.rna.columns:
                det, dep, bat, qs = _qc_scores(self.data.rna[gene].values)
                ds = directed_score(rna_eff.p_value, rna_eff.coef, rna_eff.direction)
                rows.append({
                    "event_id": eid, "modality": "rna", "assay": "RNA",
                    "target": None, "feature_id": gene, "role": "transcript_output",
                    "coef": rna_eff.coef, "se": rna_eff.se,
                    "pval": rna_eff.p_value, "fdr": rna_eff.fdr,
                    "direction": rna_eff.direction, "directed_score": ds,
                    "quality_score": qs, "detection_score": det,
                    "depth_score": dep, "batch_score": bat,
                    "model_used": _model_used(rna_eff),
                    "converged": rna_eff.convergence,
                })

            # ATAC evidence
            atac_eff = self.atac_effects.get(peak) if self.atac_effects else None
            if atac_eff and peak in self.data.atac.columns:
                det, dep, bat, qs = _qc_scores(self.data.atac[peak].values)
                ds = directed_score(atac_eff.p_value, atac_eff.coef, atac_eff.direction)
                rows.append({
                    "event_id": eid, "modality": "atac", "assay": "ATAC",
                    "target": None, "feature_id": peak, "role": "chromatin_accessibility",
                    "coef": atac_eff.coef, "se": atac_eff.se,
                    "pval": atac_eff.p_value, "fdr": atac_eff.fdr,
                    "direction": atac_eff.direction, "directed_score": ds,
                    "quality_score": qs, "detection_score": det,
                    "depth_score": dep, "batch_score": bat,
                    "model_used": _model_used(atac_eff),
                    "converged": atac_eff.convergence,
                })

            # Extra modalities
            for mod_name in self.data.modalities:
                if mod_name in ("rna", "atac"):
                    continue
                spec = self.data.modality_specs.get(mod_name)
                eff_dict = self.effects.get(mod_name, {})
                mod_eff = None
                feature = peak
                region_match = 1.0

                if spec and spec.assay == "PROTEIN":
                    links = getattr(self.data, 'protein_gene_links', None)
                    if links is not None and len(links) > 0:
                        matched = links[links["gene"].astype(str) == str(gene)]
                        if len(matched) == 0:
                            gene_short = str(gene).split(":")[0]
                            matched = links[links["gene"].astype(str) == gene_short]
                        for pid in matched["protein_id"].values:
                            mod_eff = eff_dict.get(str(pid))
                            if mod_eff:
                                feature = str(pid)
                                break
                elif spec and spec.feature_type == "region":
                    feature = peak
                    mod_eff = eff_dict.get(feature)
                    # Compute region_match_score via interval overlap
                    if mod_eff is None:
                        best_ov = 0.0
                        for k, v in eff_dict.items():
                            ov = interval_overlap(str(peak), str(k))
                            if ov and ov["min_reciprocal_overlap"] > best_ov:
                                mod_eff = v
                                region_match = ov["min_reciprocal_overlap"]
                                best_ov = ov["min_reciprocal_overlap"]
                    else:
                        ov = interval_overlap(str(peak), str(feature))
                        if ov:
                            region_match = ov["min_reciprocal_overlap"]
                else:
                    feature = gene
                    mod_eff = eff_dict.get(feature)

                if mod_eff and np.isfinite(mod_eff.coef):
                    mat = self.data.modalities.get(mod_name)
                    det, dep, bat, qs = (0.5, 0.5, 1.0, 0.5)
                    if mat is not None and feature in mat.columns:
                        det, dep, bat, qs = _qc_scores(mat[feature].values)
                    ds = directed_score(mod_eff.p_value, mod_eff.coef, mod_eff.direction)
                    rows.append({
                        "event_id": eid, "modality": mod_name,
                        "assay": spec.assay if spec else "unknown",
                        "target": spec.target if spec else None,
                        "feature_id": feature,
                        "role": spec.regulatory_role if spec else "unknown",
                        "coef": mod_eff.coef, "se": mod_eff.se,
                        "pval": mod_eff.p_value, "fdr": mod_eff.fdr,
                        "direction": mod_eff.direction, "directed_score": ds,
                        "quality_score": qs, "detection_score": det,
                        "depth_score": dep, "batch_score": bat,
                        "region_match_score": region_match,
                        "model_used": _model_used(mod_eff),
                        "converged": mod_eff.convergence,
                    })

        self.modality_evidence = pd.DataFrame(rows) if rows else pd.DataFrame()
        return self.modality_evidence

    def _build_model_diagnostics(self) -> pd.DataFrame:
        """Collect model diagnostics from all effect estimates."""
        def _add_rows(feature_effects, modality_label):
            out = []
            for feat_id, e in (feature_effects or {}).items():
                s = e.model_summary or {}
                out.append({
                    "feature_id": feat_id,
                    "modality": modality_label,
                    "model_used": s.get("model_used", "unknown"),
                    "family": s.get("family", "unknown"),
                    "alpha": s.get("alpha"),
                    "alpha_estimated": s.get("alpha_estimated", False),
                    "converged": s.get("converged", e.convergence),
                    "dropped_covariates": s.get("dropped_covariates", False),
                    "warning": s.get("warning", ""),
                })
            return out

        rows = []
        rows.extend(_add_rows(self.atac_effects, "ATAC"))
        rows.extend(_add_rows(self.rna_effects, "RNA"))
        for mod_name in self.effects:
            if mod_name in ("rna", "atac"):
                continue
            rows.extend(_add_rows(self.effects[mod_name], mod_name.upper()))
        return pd.DataFrame(rows)


def _compute_link_score(link_source: str, distance: int | float) -> float:
    """Compute peak-gene link confidence score.

    Returns 0-1 where 1 = highest confidence (promoter), lower = less certain.
    """
    if not link_source or link_source == "":
        return 0.5
    if "promoter" in link_source:
        return 1.0
    if "external" in link_source and "distal" not in link_source:
        return 0.85
    if "distal" in link_source:
        d = abs(float(distance))
        return max(0.1, 1.0 / (1.0 + d / 10000.0))
    return 0.5


def _model_used(eff) -> str:
    if isinstance(eff.model_summary, dict):
        return eff.model_summary.get("model_used", "unknown")
    return "unknown"


class MoDESResult:
    """Container for MoDES output."""

    def __init__(
        self,
        event_table: pd.DataFrame,
        state_probabilities: pd.DataFrame | None = None,
        layer_effects: pd.DataFrame | None = None,
        evidence_vectors: pd.DataFrame | None = None,
        model_diagnostics: pd.DataFrame | None = None,
        modality_evidence: pd.DataFrame | None = None,
        conditional_effects: pd.DataFrame | None = None,
        params: dict | None = None,
    ):
        self.event_table = event_table
        self.state_probabilities = state_probabilities
        self.layer_effects = layer_effects
        self.evidence_vectors = evidence_vectors
        self.model_diagnostics = model_diagnostics
        self.modality_evidence = modality_evidence
        self.conditional_effects = conditional_effects
        self.params = params or {}

    def summary(self) -> str:
        """Return a text summary of the results."""
        lines = []
        lines.append("=" * 60)
        lines.append("MoDES Results Summary")
        lines.append("=" * 60)
        lines.append(f"Total events:     {len(self.event_table)}")
        lines.append("")
        state_counts = self.event_table["state"].value_counts()
        lines.append("State distribution:")
        for state, count in state_counts.items():
            lines.append(f"  {state:25s}: {count:6d} "
                        f"({count / len(self.event_table) * 100:5.1f}%)")
        n_with_tf = self.event_table["tf_name"].notna().sum()
        if n_with_tf > 0:
            lines.append(f"\nEvents with TF annotation: {n_with_tf}")
        return "\n".join(lines)

    def filter(
        self,
        state: str | None = None,
        states: list[str] | None = None,
        min_assignment_score: float = 0.0,
        fdr_threshold: float | None = None,
        exclude_high_artifact: bool = False,
        max_event_fdr: float | None = None,
        max_state_support_adjusted_score: float | None = None,
        min_quality_score: float | None = None,
        genes: list[str] | None = None,
        peaks: list[str] | None = None,
        context: str | None = None,
        state_family: str | None = None,
    ) -> pd.DataFrame:
        """Filter the event table."""
        import warnings
        df = self.event_table.copy()

        if state_family is not None and "state_family" in df.columns:
            df = df[df["state_family"] == state_family]
        if state is not None:
            df = df[df["state"] == state]
        if states is not None:
            df = df[df["state"].isin(states)]
        if min_assignment_score > 0:
            col = "state_assignment_score"
            if col in df.columns:
                df = df[df[col] >= min_assignment_score]
        if max_event_fdr is not None:
            warnings.warn("max_event_fdr is deprecated; use max_state_support_adjusted_score",
                          DeprecationWarning)
            col = "state_support_adjusted_score" if "state_support_adjusted_score" in df.columns else "event_fdr_deprecated"
            df = df[df[col].fillna(1.0) <= max_event_fdr]
        if max_state_support_adjusted_score is not None and "state_support_adjusted_score" in df.columns:
            df = df[df["state_support_adjusted_score"] <= max_state_support_adjusted_score]
        if fdr_threshold is not None:
            df = df[(df["atac_fdr"] < fdr_threshold) | (df["rna_fdr"] < fdr_threshold)]
        if exclude_high_artifact and "artifact_risk" in df.columns:
            df = df[df["artifact_risk"] != "high"]
        if min_quality_score is not None and "quality_score" in df.columns:
            df = df[df["quality_score"] >= min_quality_score]
        if genes is not None:
            df = df[df["gene"].isin(genes)]
        if peaks is not None:
            df = df[df["peak_id"].isin(peaks)]
        if context is not None:
            df = df[df["context"] == context]
        return df.reset_index(drop=True)

    def to_tsv(self, output_dir: str) -> None:
        """Write TSV output files."""
        os.makedirs(output_dir, exist_ok=True)
        self.event_table.to_csv(
            os.path.join(output_dir, "event_table.tsv"), sep="\t", index=False)
        if self.state_probabilities is not None:
            self.state_probabilities.to_csv(
                os.path.join(output_dir, "event_state_confidence.tsv"),
                sep="\t", index=False)
        if self.layer_effects is not None:
            self.layer_effects.to_csv(
                os.path.join(output_dir, "event_layer_effects.tsv"),
                sep="\t", index=False)
        if self.evidence_vectors is not None:
            self.evidence_vectors.to_csv(
                os.path.join(output_dir, "event_evidence_vectors.tsv"),
                sep="\t", index=False)
        if self.modality_evidence is not None and len(self.modality_evidence) > 0:
            self.modality_evidence.to_csv(
                os.path.join(output_dir, "event_modality_evidence.tsv"),
                sep="\t", index=False)
        if self.model_diagnostics is not None:
            self.model_diagnostics.to_csv(
                os.path.join(output_dir, "model_diagnostics.tsv"),
                sep="\t", index=False)
        if self.conditional_effects is not None and len(self.conditional_effects) > 0:
            self.conditional_effects.to_csv(
                os.path.join(output_dir, "conditional_effects.tsv"),
                sep="\t", index=False)
        pd.DataFrame(list(self.params.items()), columns=["parameter", "value"]).to_csv(
            os.path.join(output_dir, "run_params.tsv"), sep="\t", index=False)

    def save(self, output_dir: str) -> None:
        """Save all result tables. Alias for to_tsv()."""
        self.to_tsv(output_dir)
        import json
        with open(os.path.join(output_dir, "run_params.json"), "w") as f:
            json.dump(self.params, f, indent=2, default=str)

    @classmethod
    def load(cls, output_dir: str) -> "MoDESResult":
        """Load results from a directory."""
        et = pd.read_csv(os.path.join(output_dir, "event_table.tsv"), sep="\t")
        for col in ["state", "artifact_risk", "artifact_reason", "tf_name", "context"]:
            if col in et.columns:
                et[col] = et[col].fillna("").astype(str).replace({"nan": "", "None": ""})
                if col in ("state", "artifact_risk"):
                    et[col] = et[col].replace({"": "null", "nan": "null"})
        sp = _read_optional_tsv(output_dir, "event_state_confidence.tsv")
        le = _read_optional_tsv(output_dir, "event_layer_effects.tsv")
        ev = _read_optional_tsv(output_dir, "event_evidence_vectors.tsv")
        md = _read_optional_tsv(output_dir, "model_diagnostics.tsv")
        me = _read_optional_tsv(output_dir, "event_modality_evidence.tsv")
        import json
        params = {}
        rp_json = os.path.join(output_dir, "run_params.json")
        if os.path.exists(rp_json):
            with open(rp_json) as f:
                params = json.load(f)
        elif os.path.exists(os.path.join(output_dir, "run_params.tsv")):
            rp = pd.read_csv(os.path.join(output_dir, "run_params.tsv"), sep="\t")
            params = dict(zip(rp["parameter"], rp["value"]))
        return cls(
            event_table=et, state_probabilities=sp, layer_effects=le,
            evidence_vectors=ev, model_diagnostics=md,
            modality_evidence=me, params=params)

    def to_graphml(self, path: str) -> None:
        """Write event network as GraphML."""
        try:
            import networkx as nx
        except ImportError:
            raise ImportError("networkx is required for to_graphml()")
        G = nx.Graph()
        for _, row in self.event_table.iterrows():
            gene_node = f"gene:{row['gene']}"
            peak_node = f"peak:{row['peak_id']}"
            if gene_node not in G:
                G.add_node(gene_node, type="gene")
            if peak_node not in G:
                G.add_node(peak_node, type="peak")
            G.add_edge(
                gene_node, peak_node, type="event",
                state=row["state"],
                state_assignment_score=float(row.get("state_assignment_score", np.nan)),
                artifact_risk=str(row.get("artifact_risk", "low")),
            )
            if pd.notna(row.get("tf_name")):
                tf_node = f"tf:{row['tf_name']}"
                if tf_node not in G:
                    G.add_node(tf_node, type="tf")
                G.add_edge(tf_node, peak_node, type="motif")
        nx.write_graphml(G, path)

    def to_report(self, path: str) -> None:
        """Generate HTML report."""
        from modes.report import generate_report
        generate_report(self, path)


def run_by_context(
    modes: MoDES,
    context_col: str = "context",
    min_samples_per_context: int = 4,
) -> pd.DataFrame:
    """Run MoDES separately for each context (e.g., cell type)."""
    all_events = []
    contexts = modes.data.obs[context_col].unique()
    for ctx in contexts:
        ctx_mask = modes.data.obs[context_col] == ctx
        n_ctx = ctx_mask.sum()
        if n_ctx < min_samples_per_context:
            continue
        ctx_data = MoDEData(
            rna=modes.data.rna.loc[ctx_mask],
            atac=modes.data.atac.loc[ctx_mask],
            obs=modes.data.obs.loc[ctx_mask],
        )
        ctx_modes = MoDES(
            data=ctx_data,
            condition_col=modes.condition_col,
            covariate_cols=modes.covariate_cols,
            donor_col=modes.donor_col,
            batch_col=modes.batch_col,
            fdr_threshold=modes.fdr_threshold,
            external_links=modes.external_links,
            tss_map=modes.tss_map,
            contrast=modes.contrast,
        )
        result = ctx_modes.run()
        result.event_table["context"] = ctx
        all_events.append(result.event_table)
    if all_events:
        return pd.concat(all_events, ignore_index=True)
    return pd.DataFrame(columns=_event_result_columns())
